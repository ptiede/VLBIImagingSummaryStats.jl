@info "Initiating Setup"

using Distributed

if nworkers() > 1
    @error "Please run setup with only a single process"
end

using Pkg; Pkg.activate(@__DIR__)
using Distributed
using DelimitedFiles

try
    Pkg.rm("Comrade")
catch e

end

try 
    Pkg.rm("VLBIImagingSummaryStats")
catch e

end

if isfile(joinpath(@__DIR__, "Manifest.toml"))
    @warn "Manifest file detected deleting and creating new one for setup"
    rm(joinpath(@__DIR__, "Manifest.toml"))
end

# Pkg.rm("VLBIImagingSummaryStats") #Remove the package if it is in the project otherwise you will get an error
# Pkg.rm("Comrade")
Pkg.develop(name="VLBIImagingSummaryStats", url="https://github.com/ptiede/VLBIImagingSummaryStats.jl")
Pkg.instantiate()
Pkg.precompile()
Pkg.add(;name="Comrade", rev="ptiede-enzymeswitch")

@info "Testing whether we can import `VLBIImagingSummaryStats`"
using VLBIImagingSummaryStats
@info "Successfully imported `VLBIImagingSummaryStats`"

@info "Finished setup have ready to run"
