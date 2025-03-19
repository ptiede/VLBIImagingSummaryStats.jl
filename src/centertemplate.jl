@doc raw"""
    center_template(img::IntensityMap template::Type; 
                    grid=axisdims(img), 
                    div=LeastSquares,
                    maxiters=10_000)

Recenter an image based on some template model `template`. The set of implemented templates are:
  - `Disk`: A Gaussian disk with params (r0, σ, τ, ξτ, x0, y0, f0)
  - `Gaussian`: A Gaussian with params (σ1, x0, y0, σ2, x1, y1, f2, f0)
  - `MRing{N}`: A ring with N azimuthal terms and a radial Gaussian with params (r0, σ, s, ξ, τ, ξτ, x0, y0, f0)

Parameter meaning:
    - `r0`: The radius of the template.
    - `σ`: The std dev. of the Gaussian or its radial profile.
    - `s`: The tuple of coefficients of the azimuthal expansion.
    - `ξ`: The tuple of phases of the azimuthal expansion.
    - `τ`: The stretch of the template it is always semi-minor axis / semi-major axis.
    - `ξτ`: The rotation axis of the asymmetry measured East of North
    - `xi, yi`: The center of template component.
    - `f0`: The constant background flux.

## Arguments:
 - `img::IntensityMap`: The image to recenter.
 - `template::Type`: The template model to use.
 - `grid`: The grid to use for the template matching. Default is the image grid, but smaller grids can be used to speed up the optimization.
 - `div::Function=LeastSquares`: The divergence to use for the optimization function.
 - `maxiters::Int=10_000`: The maximum number of iterations to use when optimizing for the ring center.

## Returns:
    The tuple `(img, xopt, θopt)` where:
    - `img`: The recentered image.
    - `θopt`: The optimal template
    - `xopt`: The optimal template parameters as a NamedTuple whose elements match the model above.
"""
function center_template(img, template::Type;
                         grid = axisdims(img), 
                         div = LeastSquares, 
                         maxiters = 10_000)
    rimg = regrid(img, grid)
    xopt, θopt = _center_template(rimg, template, div, maxiters)
    return shifted(img, -xopt.x0, -xopt.y0), xopt, θopt
end

function center_template(img::IntensityMap{<:StokesParams}, template::Type;
                         grid = axisdims(img), 
                         div = LeastSquares, 
                         maxiters = 10_000)
    _, xopt, θopt = center_template(stokes(img, :I), template; grid, div, maxiters)
    return shifted(img, -xopt.x0, -xopt.y0), xopt, θopt
end

function _optimize(div, func, lower, upper, p0, maxiters = 8_000)
    prob = VIDAProblem(div, func, lower, upper)
    xopt, θopt, dmin = vida(prob, ECA(;options=Options(f_calls_limit = maxiters, f_tol = 1e-5)); init_params=p0)
    # @info dmin
    return xopt, θopt
end

function _center_template(img::IntensityMap, ::Type{<:Disk}, div, maxiters)
    bh = div(max.(img, 0.0))

    x0, y0 = centroid(img)

    temp(x) = modify(VLBISkyModels.GaussDisk(x.σ/x.r0), Stretch(x.r0, x.r0*(1+x.τ)), Rotate(x.ξτ), Shift(x.x0, x.y0))  +
              x.f0*VLBISkyModels.Constant(fieldofview(img).X)
    lower = (r0 = μas2rad(10.0), σ = μas2rad(1.0),
             τ  = 0.001, ξτ = 0.0,
             x0 = -μas2rad(15.0), y0 = -μas2rad(15.0),
             f0 = 1e-6
             )
    upper = (r0 = μas2rad(40.0), σ = μas2rad(40.0),
             τ = 1.0, ξτ = 1π,
             x0 = μas2rad(15.0), y0 = μas2rad(15.0),
             f0=10.0
             )
    p0 = (r0 = μas2rad(20.0), σ = μas2rad(4.0),
          τ = 0.01, ξτ = 0.1,
          x0 = x0, y0 = y0,
          f0=1e-3)
    return _optimize(bh, temp, lower, upper, p0, maxiters)
end

@doc raw"""
    center_ring(img::IntensityMap; order=1, maxiters=10_000)

!!! note
    This function just calls center_template. We recommend using that function directly.

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
function center_ring(img::IntensityMap; order=1, g=axisdims(img), maxiters=10_000)
    return center_template(img, MRing{order}, grid=g, maxiters=maxiters)
end


function _center_template(img::IntensityMap, ::Type{<:Gaussian}, div, maxiters)
    bh = div(max.(img, 0.0))

    x0, y0 = centroid(img)
    Gauss(σ, x, y) = modify(Gaussian(), Stretch(σ), Shift(x, y))
    fovx, fovy = fieldofview(img)

    temp(x) = Gauss(x.σ1, x.x0, x.y0) +
              x.f2*Gauss(x.σ2, x.x1, x.y1) +
              x.f0*VLBISkyModels.Constant(fieldofview(img).X)
    lower = (
             σ1 = 2*min(pixelsizes(img)...), 
             x0 = -fovx/2, y0 = -fovy/2,
             σ2 = 2*min(pixelsizes(img)...), 
             x1 = -fovx/2, y1 = -fovy/2,
             f2 = 0.0,
            f0 = 1e-6
             )
    upper = (
             σ1 = max(fovx, fovy)/3, 
             x0 = fovx/2, y0 = fovy/2,
             σ2 = max(fovx, fovy)/3, 
             x1 = fovx/2, y1 = fovy/2,
             f2 = 6.0,
             f0=10.0
             )
    p0 = (
            σ1 = (upper.σ1+lower.σ1)/2, 
            x0 = x0, y0 = y0,
            σ2 = (upper.σ2+lower.σ2)/2, 
            x1 = x0, y1 = y0,
            f2 = 0.5,
            f0=1e-3
        )
    xopt, θopt = _optimize(bh, temp, lower, upper, p0, maxiters)
    flip = xopt.x0 < xopt.x1
    if flip
        return (; σ1 = xopt.σ2, x0 = xopt.x1, y0 = xopt.y1,
                  σ2 = xopt.σ1, x1 = xopt.x0, y1 = xopt.y0,
                  f2 = inv(xopt.f2), f0=xopt.f0), θopt
    end
    return xopt, θopt
end

function _center_template(img::IntensityMap{<:Real}, ::Type{<:MRing{N}}, div, maxiters) where {N}
    bh = div(max.(img, 0.0))
    x0, y0 = centroid(img)

    temp(x) = modify(RingTemplate(RadialGaussian(x.σ/x.r0), AzimuthalCosine(x.s, x.ξ .- x.ξτ)),
                     Stretch(x.r0, x.r0*(1+x.τ)), Rotate(x.ξτ), Shift(x.x0, x.y0))+
            #   modify(Gaussian(), Stretch(x.σg), Shift(x.xg, x.yg), Renormalize(x.fg)) +
              x.f0*VLBISkyModels.Constant(fieldofview(img).X)
    lower = (r0 = μas2rad(10.0), σ = μas2rad(0.5),
             s = ntuple(_->0.001, N),
             ξ = ntuple(_->0.0, N),
             τ = 0.0,
             ξτ = 0.0,
             x0 = -μas2rad(20.0), y0 = -μas2rad(20.0),
            #  σg = μas2rad(30.0),
            #  xg = -fieldofview(img).X/4,
            #  yg = -fieldofview(img).Y/4,
            #  fg = 1e-6,
             f0 = 1e-6
             )
    upper = (r0 = μas2rad(30.0), σ = μas2rad(15.0),
             s = ntuple(_->0.999, N),
             ξ = ntuple(_->2π, N),
             τ = 1.0,
             ξτ = 1π,
             x0 = μas2rad(20.0), y0 = μas2rad(20.0),
            #  σg = fieldofview(img).X/2,
            #  xg = fieldofview(img).X/4,
            #  yg = fieldofview(img).Y/4,
            #  fg = 20.0,
             f0=10.0
             )
    p0 = (r0 = μas2rad(16.0), σ = μas2rad(4.0),
          s = ntuple(_->0.2, N),
          ξ = ntuple(_->1π, N),
          τ = 0.01,
          ξτ = 0.5π,
          x0 = x0, y0 = y0,
        #   σg = μas2rad(40.0), xg = 0.0, yg = 0.0, fg = 0.2,
          f0=0.1)
    return _optimize(bh, temp, lower, upper, p0, maxiters)
end
