using Distributed
@everywhere begin
    const filedir = @__DIR__
end

using Pkg;
Pkg.activate(filedir);
using Distributed
using DelimitedFiles

@everywhere begin
    using Pkg
    Pkg.activate(filedir)
end

using Comonicon
using VLBIImagingSummaryStats
using DataFrames
using CSV
using VLBISkyModels
@everywhere begin
    using VIDA
    using VLBISkyModels
    using VLBIImagingSummaryStats
end

function loaddir(file)
    return open(file) do io
        return readlines(io)
    end
end

"""
Extract summary statistics from a set of images.

# Arguments

- `imfiles`: The file containing the absolute paths to all the GRMHD images we will analyze
- `outname`: The base file name for the output files

# Options

- `-c, --code=<string>`: The file containing the code used to generate the images. If nothing will ignore.
- `-f, --fevals=<int>`: The number of evaluations allowed when extracting params
- `-s, --stride=<int>`: The checkpointing stride
- `-o, --order=<int>`: The order of the ring parameters
- `-b, --blur=<float>`: Gaussian FWHM blurring kernel to apply to the images in μas

# Flags

- `-r, --regrid`: Regrid the images before extracting
- `--restart`: Restart the extraction process from before
"""
@main function main(imfiles::String, outname::String; fevals::Int = 40_000, stride::Int = 2 * nworkers(), code::String = "", order::Int = 4, blur::Float64 = 0.0, regrid::Bool = false, restart::Bool = false)
    @info "Image files path: $(imfiles)"
    @info "Outputting results to $(outname)"
    @info "Using a $(order) order ring model"

    imfs = loaddir(imfiles)
    @info "Loaded $(length(imfs)) files"
    if code != ""
        @info "Code file: $(code)"
        cfs = loaddir(code)
    else
        cfs = fill("unknown", length(imfs))
    end

    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 50, 50)

    @info "Regridding image : $(regrid)"
    @info "Blurring kernel: $(blur) μas"

    @assert length(imfs) == length(cfs) "Length of imfiles ($(length(imfs))) and code ($(length(cfs))) do not match"

    @info "Starting to extract summary statistics"
    if restart
        df = CSV.read(outname, DataFrame)
    else
        df = DataFrame()
    end
    startindex = nrow(df) + 1
    indexpart = Iterators.partition(startindex:length(imfs), stride)
    for ii in indexpart
        @info "Extracting from $(ii[begin]) to $(ii[end])"
        res = pmap(imfs[ii]) do f
            img = center_image(load_image(f; polarization = true))

            if blur > 0.0
                img = smooth(img, μas2rad(blur) / (2 * sqrt(2 * log(2))))
            end

            if regrid
                gim = g
            else
                gim = nothing
            end
            rimg = img
            

            stats = summary_ringparams(rimg; maxiters = fevals, order, divergence = NxCorr, grid=gim)
            return stats
        end
        dftmp = DataFrame(res)
        dftmp.code = cfs[ii]
        dftmp.files = imfs[ii]
        df = vcat(df, dftmp)
        CSV.write(outname, df)
    end
end
