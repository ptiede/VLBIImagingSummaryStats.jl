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

    m_avg = mavg(simg)
    v_avg = vavg(simg)

    βlp = lpmodes(simg, lpmode)
    βcp = cpmodes(simg, cpmode)

    n_relp = Tuple((Symbol("re_betalp", "_$n") for n in lpmode))
    real_betalp = NamedTuple{n_relp}(map(real, βlp))
    n_imlp = Tuple((Symbol("im_betalp", "_$n") for n in lpmode))
    imag_betalp = NamedTuple{n_imlp}(map(imag, βlp))

    n_recp = Tuple((Symbol("re_betacp", "_$n") for n in cpmode))
    real_betacp = NamedTuple{n_recp}(map(real, βcp))
    n_imcp = Tuple((Symbol("im_betacp", "_$n") for n in cpmode))
    imag_betacp = NamedTuple{n_imcp}(map(imag, βcp))

    return merge(xopt, (;m_net, m_avg), real_betalp, imag_betalp, (;v_net, v_avg), real_betacp, imag_betacp)
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
