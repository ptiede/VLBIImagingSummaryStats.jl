"""
    summary_ringparams(img::IntensityMap;
                        lpmode=(2,), cpmode=(1,),
                        order=1, maxiters=1000)

Analyes the image `img` and produces summary statistics about the it.
For stokes I this includes all the ring parameters from [`center_ring`](@ref).
If the image is polarized this also includes the linear and circular polarization modes
beta modes given by `lpmode` and `cpmode` respectively, as well as `mnet` and `vnet`.
"""
function summary_ringparams(img::IntensityMap{<:StokesParams};
                            lpmode=(2,), cpmode=(1,),
                            order=1, maxiters=1000)
    xopt = summary_ringparams(stokes(img, :I); order, maxiters)
    simg = shifted(img, -xopt.x0, -xopt.y0)
    m_net = mnet(simg)
    v_net = vnet(simg)

    βlp = lpmodes(simg, lpmode)
    βcp = cpmodes(simg, cpmode)

    n_amplp = Tuple((Symbol("amp_betalp", "_$n") for n in lpmode))
    amp_betalp = NamedTuple{n_amplp}(map(abs, βlp))
    n_arglp = Tuple((Symbol("arg_betalp", "_$n") for n in lpmode))
    arg_betalp = NamedTuple{n_arglp}(map(angle, βlp))

    n_ampcp = Tuple((Symbol("amp_betacp", "_$n") for n in cpmode))
    amp_betacp = NamedTuple{n_ampcp}(map(abs, βcp))
    n_argcp = Tuple((Symbol("arg_betacp", "_$n") for n in cpmode))
    arg_betacp = NamedTuple{n_argcp}(map(angle, βcp))

    return merge(xopt, (;m_net), amp_betalp, arg_betalp, (;v_net), amp_betacp, arg_betacp)
end

function summary_ringparams(img::IntensityMap{<:Real};
                            lpmodes=(2,), cpmodes=(1,),
                            order=1, maxiters=1000)
    rimg = regrid(img, imagepixels(μas2rad(120.0), μas2rad(120.0), 32, 32))
    _, xopt, _ = center_ring(rimg; order, maxiters)
    return _flatten_tuple(xopt)
end

function _flatten_tuple(nt::NamedTuple)
    names = keys(nt)
    vals = values(nt)
    mapreduce(merge, zip(names, vals)) do (name, value)
        return _flatten_tuple(name, value)
    end
end

function _flatten_tuple(name::Symbol, a)
    return NamedTuple{(name,)}((a,))
end

function _flatten_tuple(name::Symbol, t::NTuple{N}) where {N}
    names = Tuple((Symbol("$(name)_$i") for i in 1:N))
    return NamedTuple{names}(t)
end
