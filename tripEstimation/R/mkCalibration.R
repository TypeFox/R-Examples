"mkCalibration" <-
function (x, known = NULL, elim = c(-36, 12), choose = TRUE)
{
    ## removed by addition of namespace MDS 2011-10-06
    ## require(mgcv) || stop("required package mgcv not available")

    if (is.null(known))
        stop("known location must be provided")
    elevation.gmt <- mkElevationSeg(1, x$gmt)
    x$elevation <- elevation.gmt(1, known[1], known[2])
    keep <- x$elevation >= elim[1] & x$elevation <= elim[2]
    x$segment <- factor(cumsum(c(0, abs(diff(keep)))))
    x <- x[keep, ]

    if (choose)  {

      ## removed by addition of namespace MDS 2011-10-06
      ## require(lattice) || stop("required package lattice not available")

      print(xyplot(light ~ elevation | segment, data = x))

      cat("\n")
      print("Enter those segments (without quotes) which look OK, finish by entering nothing")
      cat("\n")
      keepsegs <- scan("", "numeric")
      x <- x[x$segment %in% keepsegs, ]
      x$segment <- factor(x$segment)
    }

    ## HUH?  elevation.gmt <- mkElevationSeg(1, x$gmt)

    if (nlevels(x$segment) == 1) {
    	 fit <- gam(light ~ s(elevation, k = 20), family = gaussian(link = "identity"),
        	data = x)
        	warning("calibration made from a single segment, ok but let it be noted\n\n")

        } else {
	    fit <- gam(light ~ s(elevation, k = 20) + segment, family = gaussian(link = "identity"),
	        data = x)
	        }
    elev.cal <- seq(min(x$elevation), max(x$elevation),
        length = 500)
    new.data <- data.frame(elevation = elev.cal, segment = rep(x$segment[1],
        length(elev.cal)))
    light.cal <- predict(fit, new.data, type = "response")
    approxfun(elev.cal, light.cal, rule = 2)
}

