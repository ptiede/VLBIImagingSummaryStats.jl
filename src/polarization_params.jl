"""
    lpmodes(img::IntensityMap{<:StokesParams}, modes::NTuple{N, Int}; rmin=0, rmax=μas2rad(60.0))

Computes the linear polarization modes of an image.

## Arguments:
 - `img::IntensityMap{<:StokesParams}`: The image to analyze.
 - `modes::NTuple{N, Int}`: The modes to compute.

## Keyword Arguments:
    - `rmin::Real=0`: The minimum radius to consider.
    - `rmax::Real=μas2rad(60.0)`: The maximum radius to consider.

## Returns:
    - `betas::NTuple{N,Complex}`: The coefficients of the modes.

## Warning:

These modes implicitly assume that the image is centered on the ring.
If the image is not centered, the results will be incorrect. If you
want to center the image, use `center_template` before calling this function.

## Example:

```julia-repl
julia> img = load_image("image.fits")
julia> cimg = center_template(img)
julia> βs = lpmodes(cimg, (1, 2)) # compute the first 2 modes
```
"""
function lpmodes(img::IntensityMap{<:StokesParams}, modes::NTuple{<:Any, Int}; rmin=0, rmax=μas2rad(60.0))
    I = stokes(img, :I)
    lp = linearpol.(img)
    g = domainpoints(img)
    r = splat(hypot).(values.(g))
    ang(p) = atan(-p.X, p.Y)
    θ = ang.(g)

    mask = (r .<= rmax).&(r .>= rmin)
    flux = sum(I[mask])
    betas = map(modes) do m
        pbasis = cis.(m.*θ)
        prod = lp.*pbasis
        coef = sum(prod[mask])
        return coef/flux
    end
    return betas
end

"""
    cpmodes(img::IntensityMap{<:StokesParams}, modes::NTuple{<:Any, Int}; rmin=0, rmax=μas2rad(60.0))

Computes the circular polarization modes of an image.

## Arguments:
 - `img::IntensityMap{<:StokesParams}`: The image to analyze.
 - `modes::NTuple{<:Any, Int}`: The modes to compute.

## Keyword Arguments:
    - `rmin::Real=0`: The minimum radius to consider.
    - `rmax::Real=μas2rad(60.0)`: The maximum radius to consider.

## Returns:
    - `betas::NTuple{Complex}`: The coefficients of the modes.

## Warning:

These modes implicitly assume that the image is centered on the ring.
If the image is not centered, the results will be incorrect. If you
want to center the image, use `center_template` before calling this function.

## Example:

```julia-repl
julia> img = load_image("image.fits")
julia> cimg = center_template(img, MRing{2})
julia> βs = cpmodes(cimg, (1, 2)) # compute the first 2 modes
```
"""
function cpmodes(img::IntensityMap{<:StokesParams}, modes::NTuple{<:Any, Int}; rmin=0, rmax=μas2rad(60.0))
    I = stokes(img, :I)
    lp = stokes(img, :V)
    g = domainpoints(img)
    r = splat(hypot).(values.(g))
    ang(p) = atan(-p.X, p.Y)
    θ = ang.(g)

    mask = (r .<= rmax).&(r .>= rmin)
    flux = sum(I[mask])
    betas = map(modes) do m
        pbasis = cis.(m.*θ)
        prod = lp.*pbasis
        coef = sum(prod[mask])
        return coef/flux
    end
    return betas
end

"""
    mnet(img::IntensityMap{<:StokesParams})

Compute the net fractional linear polarization of an image.
"""
mnet(img::AbstractArray{<:StokesParams}) = abs(sum(linearpol, img))/sum(stokes(img, :I))

"""
    mnet(img::IntensityMap{<:StokesParams})

Compute the net fractional linear polarization of an image.
"""
vnet(img::AbstractArray{<:StokesParams}) = sum(stokes(img, :V))/sum(stokes(img, :I))

"""
    mavg(img::IntensityMap{<:StokesParams})

Compute the average fractional linear polarization of an image.

## Warning:
This is a biased estimator that depends on the effective resolution of the image and the
field of view. To make comparissons with data you should first debias this result based
on synthetic data tests.
"""
mavg(img::AbstractArray{<:StokesParams}) = sum(abs.(linearpol.(img)))/sum(stokes(img, :I))

"""
    vavg(img::IntensityMap{<:StokesParams})

Compute the average fractional circular polarization of an image.

## Warning:
This is a biased estimator that depends on the effective resolution of the image and the
field of view. To make comparissons with data you should first debias this result based
on synthetic data tests.
"""
vavg(img::AbstractArray{<:StokesParams}) = sum(abs.(stokes(img, :V)))/sum(stokes(img, :I))


"""
    netevpa(img::IntensityMap{<:StokesParams})

Computes the image integrated EVPA, which we define as the `evpa(mean(img))`
"""
netevpa(img::AbstractArray{<:StokesParams}) = evpa(mean(img))


"""
    polnxcorr(x::IntensityMap{<:StokesParams}, y::IntensityMap{<:StokesParams})

Computes the polarization normalized cross-correlation of two images. This returns a named tuple
with fields `LP` and `V` which are the normalized cross-correlation of the linear polarization
and circular polarization.
"""
function polnxcorr(x::IntensityMap{<:StokesParams{T}}, y::IntensityMap{<:StokesParams}) where {T}
    sp = zero(Complex{T})
    sv = zero(T)
    nx = zero(T)
    ny = zero(T)
    nvx= zero(T)
    nvy= zero(T)
    for I in eachindex(x, y)
        xx = x[I]
        yy = y[I]
        lpx = linearpol(xx)
        lpy = linearpol(yy)
        sp += lpx*conj(lpy)
        nx += abs2(lpx)
        ny += abs2(lpy)
        sv += xx.V*yy.V
        nvx += abs2(xx.V)
        nvy += abs2(yy.V)
    end
    return (LP = real(sp)*inv(sqrt(nx*ny)),  V =sv*inv(sqrt(nvx*nvy)))
end

"""
    nxcorr(x::IntensityMap{<:StokesParams}, y::IntensityMap{<:StokesParams})

Computes the polarized normalized cross-correlation of two images. This returns a named tuple
with fields `I`, `LP` and `V` which are the normalized cross-correlation of the total intensity,
linear polarization and circular polarization.

The linear polarization NxCorr is defined as 

    LP(P1, P2) = real(<P1, P2>) / sqrt(<P1, P1> <P2, P2>)

where <P1, P2> is the dot product of the complex linear polarization vectors.

The circular polarization nxcorr of two circular polarization maps V1, V2 is defined as 

    V(V1, V2) =  <V1,V2>/ sqrt(<V1,V1> <V2,V2>)

where <V1V2> is the dot product of the circular polarization vectors.

"""
function VIDA.nxcorr(x::IntensityMap{<:StokesParams}, y::IntensityMap{<:StokesParams}) 
    nI = nxcorr(stokes(x, :I), stokes(y, :I))
    nP = polnxcorr(x, y)
    return merge((I=nI,), nP)
end
