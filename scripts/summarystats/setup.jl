@info "Initiating Setup"

using Distributed

if nworkers() > 1
    @error "Please run setup with only a single process"
end

using Pkg; Pkg.activate(@__DIR__)
using Distributed
using DelimitedFiles


Pkg.develop(name="VLBIImagingSummaryStats")
Pkg.instantiate()
Pkg.precompile()
Pkg.add(;name="Comrade", rev="ptiede-enzymeswitch")

@info "Testing whether we can import `VLBIImagingSummaryStats`"
using VLBIImagingSummaryStats
@info "Successfully imported `VLBIImagingSummaryStats`"

@info "Finished setup have ready to run"
