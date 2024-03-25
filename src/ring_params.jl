@doc raw"""
    center_ring(img::IntensityMap; order=1, maxiters=10_000)

Recenter an image based on a ring model. This will attempt to find the center of an image
whose dominant component is a ring, and then translate the image so that the ring center
is at the origin.

To center the ring we use `VIDA` and a ring model with a radial Gaussian and an azimuthal
cosine expansions. The model is given by:

```math
I(x, y) = \sum_{n=0}^{order} s_n \cos(n\theta + \xi_n) \exp\left(-\frac{r^2}{2\sigma^2}\right) + f_0 + f_g \exp\left(-\frac{r^2}{2\sigma_g^2}\right)
```

where `r = \sqrt{x^2 + y^2}` and `θ = \arctan(-x, y)`. The parameters are:

 - `r0`: The radius of the ring.
 - `σ`: The width of the ring.
 - `s`: The tuple of coefficients of the azimuthal expansion.
 - `ξ`: The tuple of phases of the azimuthal expansion.
 - `x0, y0`: The center of the ring.
 - `σg`: The width of the Gaussian background.
 - `xg, yg`: The center of the Gaussian background.
 - `fg`: The amplitude of the Gaussian background.
 - `f0`: The constant background.

## Arguments:
 - `img::IntensityMap`: The image to recenter.
 - `order::Int=1`: The order of the ring model to use.
 - `maxiters::Int=10_000`: The maximum number of iterations to use when optimizing for the ring center.

## Returns:
    The tuple `(θopt, xopt, img)` where:
    - `θopt`: The optimal ring model
    - `xopt`: The optimal ring model parameters as a NamedTuple whose elements match the model above.
    - `img`: The recentered image.
"""
function center_ring(img::IntensityMap{<:Real}; order=1, maxiters=10_000)
    xopt, θopt = find_ring_center(img; order, maxiters)
    return θopt, xopt, shifted(img, -xopt.x0, -xopt.y0)
end

function center_ring(img::IntensityMap{<:StokesParams}; order=1, maxiters=10_000)
    xopt, θopt = find_ring_center(stokes(img, :I); order, maxiters)
    return θopt, xopt, shifted(img, -xopt.x0, -xopt.y0)
end

function find_ring_center(img::IntensityMap{<:Real}; order=1, maxiters=10_000)
    bh = LeastSquares(max.(img, 0.0))

    x0, y0 = centroid(img)

    ring(x) = modify(RingTemplate(RadialGaussian(x.σ/x.r0), AzimuthalCosine(x.s, x.ξ)),
                     Stretch(x.r0, x.r0), Shift(x.x0, x.y0))+
              modify(Gaussian(), Stretch(x.σg), Shift(x.xg, x.yg), Renormalize(x.fg)) +
              x.f0*Constant(fieldofview(img).X)
    lower = (r0 = μas2rad(15.0), σ = μas2rad(1.0),
             s = ntuple(_->0.001, order),
             ξ = ntuple(_->0.0, order),
             x0 = -μas2rad(15.0), y0 = -μas2rad(15.0),
             σg = μas2rad(30.0),
             xg = -fieldofview(img).X/4,
             yg = -fieldofview(img).Y/4,
             fg = 1e-6, f0 = 1e-6
             )
    upper = (r0 = μas2rad(30.0), σ = μas2rad(8.0),
             s = ntuple(_->0.999, order),
             ξ = ntuple(_->2π, order),
             x0 = μas2rad(15.0), y0 = μas2rad(15.0),
             σg = fieldofview(img).X/2,
             xg = fieldofview(img).X/4,
             yg = fieldofview(img).Y/4,
             fg = 20.0, f0=10.0
             )
    p0 = (r0 = μas2rad(20.0), σ = μas2rad(4.0),
             s = ntuple(_->0.01, order),
             ξ = ntuple(_->1π, order),
            x0 = x0, y0 = y0,
          σg = μas2rad(40.0), xg = 0.0, yg = 0.0, fg = 0.2, f0=1e-3)
    prob = VIDAProblem(bh, ring, lower, upper)
    xopt, θopt, _ = vida(prob, ECA(;options=Options(f_calls_limit = maxiters, f_tol = 1e-5)); init_params=p0)
    return xopt, θopt
end

"""
    match_center_and_res(target::IntensityMap, input::IntensityMap)

Blur and shift the center of an image `input` to match the center and resolution
of the image `target`. This is useful when wanting to compare images from different
codes that may be different intrinsic resolutions. The function will use `VIDA` to
find the optimal shift and blurring to match the two images based on the NXCORR metric.

## Arguments:
 - `target::IntensityMap`: The target image to match.
 - `input::IntensityMap`: The image to shift and blur.

## Returns:
    The tuple `(img, xopt)` where:
    - `img`: The shifted and blurred image.
    - `xopt`: The optimal shift and blurring parameters as a NamedTuple with fields `x`, `y`, and `σ`.
"""
function match_center_and_res(target::IntensityMap{<:StokesParams}, input::IntensityMap{<:StokesParams})
    target_I = stokes(target, :I)
    input_I = stokes(input, :I)
    _, xopt = match_center_and_res(target_I, input_I)
    return Comrade.smooth(shifted(input, xopt.x, xopt.y), xopt.σ), xopt
end

function match_center_and_res(target::IntensityMap, input::IntensityMap)
    cache = create_cache(FFTAlg(), input)
    f = let cimg=ContinuousImage(input, cache)
        function f(x)
            return smoothed(shifted(cimg, x.x, x.y), x.σ)
        end
    end
    div = NxCorr(target)
    lower = (x=-μas2rad(20.0), y=-μas2rad(20.0), σ=μas2rad(1.0))
    upper = (x=μas2rad(20.0), y=μas2rad(20.0), σ=μas2rad(30.0))
    p0 = (x=0.0, y=0.0, σ=μas2rad(10.0))

    prob = VIDAProblem(div, f, lower, upper)
    xopt, θopt, divmin = vida(prob, ECA(;options=Options(f_calls_limit = 2000, f_tol = 1e-5)); init_params=p0, maxiters=800)
    return Comrade.smooth(shifted(input, xopt.x, xopt.y), xopt.σ), xopt
end
