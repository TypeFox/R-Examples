## plot x, y, and z,
## with plane fit and display of squared residuals
regr2.plot <- function(x, y, z,
                       main.in="put a useful title here",
                       resid.plot=FALSE,
                       plot.base.plane=TRUE,
                       plot.back.planes=TRUE,
                       plot.base.points=FALSE,
                       eye=NULL,                   ## S-Plus
                       theta=0, phi=15, r=sqrt(3), ticktype="detailed", ## R
                       ...) {

  if.R(r=perspp <- trans3d, s={})

  ## least squares fit
  zxy.lm <- lm(z  ~ x +  y)
  fit.points <- predict(zxy.lm)

  ## fitted plane
  fit.plane <- predict(zxy.lm,
                       expand.grid(x=pretty(x),
                                   y=pretty(y)))
  fit.plane <- matrix(fit.plane,
                      length(pretty(x)),
                      length(pretty(y)))

  ## 3D display box and axes
  if.R(r={},
       s=persp.setup(col=c(0,0,0)) # set to background color
       )
  persp.out <- if.R(s={
    if (is.null(eye))
      persp(pretty(x),
            pretty(y),
            fit.plane,
            ...)
    else
      persp(pretty(x),
            pretty(y),
            fit.plane,
            eye=eye,
            ...)
  },r={
    persp.out <- persp(pretty(x),
                       pretty(y),
                       fit.plane,
                       theta=theta, phi=phi, r=r, ticktype=ticktype,
                       ...)
  })
  title(main=main.in)
  if.R(r={},
       s=persp.setup(restore=TRUE)    # restore default values
       )
  
  ## 3D observed points
  persp.points <- if.R(s=
                       perspp(x, y, z, persp.out)
                       ,r=
                       trans3d(x, y, z, persp.out)
                       )
  points(persp.points, pch=16, col=1)

  ## 3D fitted points
  persp.fit.points <- if.R(s=
                           perspp(x, y, fit.points, persp.out)
                           ,r=
                           trans3d(x, y, fit.points, persp.out)
                           )
  points(persp.fit.points, pch=9, col=2)

  ## 3D fitted plane
  perspPlane(pretty(x), pretty(y), fit.plane, persp.out,
              col=2, lty=3)

  if (plot.base.points) {
    ## 3D base plane points
    persp.base.points <-
      perspp(x, y, fit.points*0, persp.out)
    points(persp.base.points, pch="o")

    ## 3D lines from base to fitted
    segments(persp.base.points$x, persp.base.points$y,
             persp.fit.points$x , persp.fit.points$y,
             col=2)
  }

  ## 3D residual lines from fitted to observed
  segments(persp.fit.points$x, persp.fit.points$y,
           persp.points$x,     persp.points$y,     col=2)


  ## 3D base plane and back walls
  if (plot.base.plane)
    perspFloor(pretty(x), pretty(y), fit.plane*0,
                persp.out,
                col=1, lty=2)
  if (plot.back.planes) {
    perspBack.wall.x(pretty(x), pretty(y), fit.plane,
                      persp.out,
                      col=1, lty=2)
    perspBack.wall.y(pretty(x), pretty(y), fit.plane,
                      persp.out,
                      col=1, lty=2)
  }


  ## 3D squared residuals in visual inches
  if (resid.plot != FALSE)
    resid.squares(persp.points$x, persp.points$y,
                  persp.fit.points$y, resid.plot)
  invisible(NULL)
}
