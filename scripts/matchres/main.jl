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
- `bimg`: The base image to match the resolution
- `outdir`: The output directory for the results. Note the image names will be identical to those in `imfiles`

# Options

- `-s, --stride=<int>`: The checkpointing stride

# Flags

- `--restart`: Restart the extraction process from before
- `--regrid`: Regrid the images before extracting
"""
@main function main(imfiles::String, bimg::String, outdir::String; stride::Int=2 * nworkers(), regrid::Bool=false, restart::Bool=false)
    @info "Image files path: $(imfiles)"
    @info "Base image to match res $(bimg)"
    @info "Outputting results to $(outdir)"
    !isfile(bimg) && throw(ArgumentError("Base image to match resolution $(bimg) does not exist"))

    imfs = loaddir(imfiles)
    @info "Loaded $(length(imfs)) files"

    mkpath(outdir)

    g = imagepixels(μas2rad(150.0), μas2rad(150.0), 50, 50)

    # Flip this because by default we want to regrid
    regrid = !regrid

    if regrid
        bcimg = VLBISkyModels.regrid(load_image(bimg; polarization=true), g)
    else
        bcimg = load_image(bimg; polarization=true)
    end

    @info "Regridding image : $(regrid)"

    @info "Starting to extract summary statistics"
    if restart
        df = CSV.read(joinpath(outdir, "parameters.csv"), DataFrame)
    else
        df = DataFrame()
    end
    startindex = nrow(df) + 1
    indexpart = Iterators.partition(startindex:length(imfs), stride)
    for ii in indexpart
        @info "Extracting from $(ii[begin]) to $(ii[end])"
        res = pmap(imfs[ii]) do f
            img = center_image(load_image(f; polarization=true))

            if regrid
                rimg = VLBISkyModels.regrid(img, g)
            else
                rimg = img
            end

            cimg, xopt = match_center_and_res(stokes(bcimg, :I), img; maxiters=15_000, grid=g, divergence=NxCorr)
            # now regrid after blurring to prevent interpolation error
            rcimg = VLBISkyModels.regrid(cimg, g)
            nxc = nxcorr(rcimg, bcimg)
            xopt = merge(xopt, (; nxI=nxc.I, nxP=nxc.LP, nxV=nxc.V))
            save_fits(joinpath(outdir, replace(basename(f), ".fits" => "_matchres.fits")), rcimg)
            return xopt
        end

        dftmp = DataFrame(res)
        dftmp.files = imfs[ii]
        df = vcat(df, dftmp)
        CSV.write(joinpath(outdir, "parameters.csv"), df)
    end
end
