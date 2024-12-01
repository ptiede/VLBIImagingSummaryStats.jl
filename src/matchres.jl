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
function match_center_and_res(target::IntensityMap, input::IntensityMap; divergence=NxCorr, maxiters=8000)

    if eltype(target) <: StokesParams
        timg = stokes(target, :I)
    else
        timg = target
    end

    if eltype(input) <: StokesParams
        iimg = stokes(input, :I)
    else
        iimg = input
    end


    f = let cimg=VLBISkyModels.InterpolatedImage(iimg)
        function f(x)
            return smoothed(shifted(cimg, x.x, x.y), x.σ)
        end
    end
    div = divergence(timg)
    lower = (x=-μas2rad(20.0), y=-μas2rad(20.0), σ=μas2rad(1e-2))
    upper = (x=μas2rad(20.0), y=μas2rad(20.0), σ=μas2rad(30.0))
    p0 = (x=0.0, y=0.0, σ=μas2rad(10.0))

    prob = VIDAProblem(div, f, lower, upper)
    xopt, θopt, divmin = vida(prob, ECA(;options=Options(f_calls_limit = 8000, f_tol = 1e-5)); 
                              init_params=p0, maxiters
                              )
    return smooth(shifted(input, xopt.x, xopt.y), xopt.σ), xopt
end
