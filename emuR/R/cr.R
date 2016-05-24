##' Plot digital sinuoids.
##' 
##' The function plots and/or sums digital sinusoids for different parameter
##' settings.
##' 
##' 
##' @param A A vector of amplitude values. Defaults to A = 1
##' @param k A vector of cycles (repetitions). Defauls to k = 1
##' @param p A vector of phase values between -pi/2 and pi/2. Defaults to 0.
##' @param N The number of points in the signal. Defaults to 16.
##' @param samfreq If NULL, then a sinusoid is plotted with a frequency of k
##' cycles per N points. Otherwise, if samfreq is an numeric, then the argument
##' to k is interpreted as the frequency in Hz and the sinusoid at that
##' frequency is plotted for however many points are specified by N. For
##' example, if samfreq is 40 (Hz), and if N is 40 and k = 1, then 1 cycle of a
##' 1 Hz sinusoid will be plotted.
##' @param duration Specify the duration in ms. If NULL, the default, then the
##' duration of the sinusoid is in points (N), otherwise if a numeric value is
##' supplied, then in ms. For example, 1/2 second of a 1 cycle sinusoid at a
##' sampling frequency of 40 Hz: duration = 500, k = 1, samfreq=40. A ms value
##' can be supplied only if the sampling frequency is also specified.
##' @param const A single numeric vector for shifting the entire sinusoid up or
##' down the y-axis. For example, when const is 5, then 5 + s, where s is the
##' sinusoid is plotted. Defaults to 0 (zero).
##' @param expon A numeric vector. If supplied, then what is plotted is
##' expon[j]\eqn{\mbox{\textasciicircum}}{^}(c(0:(N - 1) * A cos (2 * pi * k/N
##' * (0:(N-1))). For example, a decaying sinusoid is produced with
##' cr(expon=-0.9). Defaults to NULL (i.e. to expon = 1).
##' @param plotf A single-valued logical vector. If T (default), the sinusoid
##' is plotted.
##' @param ylim A two-valued numeric vector for specifying the y-axis range.
##' @param xlim A two-valued numeric vector for specifying the y-axis range.
##' @param values If T, then the values of the sinusoid are listed. Defaults to
##' F.
##' @param xlab A character vector for plotting the x-axis title.
##' @param ylab A character vector for plotting the y-axis title.
##' @param type A character vector for specifying the line type (see par)
##' @param bw A numeric vector for specifying the bandwidth, if the sampling
##' frequency is supplied. The bandwidth is converted to an exponential (see
##' expon using exp( - rad(bw/2, samfreq = samfreq).
##' @param dopoints this is now redundant.
##' @param \dots Option for supplying further graphical parameters - see par.
##' @author Jonathan Harrington
##' @seealso \code{\link{crplot}}
##' @keywords dplot
##' @examples
##' 
##' # cosine wave
##' cr()
##' 
##' # doubling the frequency, 1/3 amplitude, phase = pi/4, 50 points
##' cr(A=1/3, k=2, p=pi/4, N=50)
##' 
##' # sum 3 sinusoids of different frequencies)
##' cr(k=c(1, 3, 4))
##' 
##' # sum 2 sinusoids of different parameters
##' cr(c(1, 2), c(2, 10), c(0, -pi/3), N=200, type="l")
##' 
##' 
##' # store the above to a vector and overlay with noise
##' v = cr(c(1, 2), c(2, 10), c(0, -pi/3), N=200, type="l", values=TRUE)
##' r = runif(200, -3, 3)
##' v = v+r
##' plot(0:199, v, type="l")
##' 
##' 
##' # 100 points of a 50 Hz sinusoid with a 4 Hz bandwidth 
##' # at a sampling frequency of 200 Hz
##' cr(k=50, bw=4, samfreq=2000, N=100)
##' 
##' # the same but shift the y-axis by +4 (d.c. offset=+4)
##' cr(const=4, k=50, bw=4, samfreq=2000, N=100)
##' 
##' # sinusoid multiplied by a decaying exponential (same effect as bandwidth)
##' cr(expon=-0.95,  N=200, type="l")
##' 
##' 
##' @export cr
"cr" <- function(A = 1, k = 1, p = 0, N = 16, 
                 samfreq = NULL, duration = NULL, 
                 const = NULL, expon = NULL, plotf = TRUE, 
                 ylim = NULL,  xlim=NULL, values = FALSE, 
                 xlab = "Time (number of points)", ylab = "Amplitude", 
                 type = "b", bw = NULL, dopoints = FALSE, ...)
{
  ## A: amplitude (arbitrary units); defaults to 1
  ## to shift by a quarter wave, p = 2 * pi * .25
  ## k: Number of sine wave k; if k is n, then
  ## n cycles will be plotted
  ## samfreq: (sampling frequency in Hz)
  ## if NULL, then a sinusoid is plotted
  ## with a frequency of k cycles per N points.
  ## otherwise, if samfreq is an integer value, 
  ## then the argument to k is interpreted
  ## as the frequency in Hz and the sinusoid at
  ## that frequency is plotted for however many
  ## points are specified by N
  ## duration: the desired duration in ms
  ## of the sinusoid. If not supplied, the
  ## duration is N, the number of points
  ## bw: supply a bandwidth in Hz; only 
  ## works if a sampling frequency is also specified
  ## dopoints: plot superimposed large data points on a line
  # 
  vec.n <- c(length(A), length(k),  length(p))
  if(!all(diff(vec.n)==0))
  {
    which <- vec.n == max(vec.n)
    if(!which[1])
      A <- rep(A, max(vec.n) )
    if(!which[2])
      k <- rep(k, max(vec.n) )
    if(!which[3])
      p <- rep(p, max(vec.n) )
  }
  
  if(dopoints) type <- "l"
  if(!is.null(duration)) {
    if(is.null(samfreq)) {
      print("must supply a sampling frequency if duration supplied"
      )
      stop()
    }
    N <- round((samfreq * duration)/1000)
  }
  if(!is.null(bw)) {
    if(is.null(samfreq)) {
      print("must supply a sampling frequency if bandwidth supplied"
      )
      stop()
    }
    expon <- exp( - rad(bw/2, samfreq = samfreq))
  }
  if(!is.null(samfreq) & is.numeric(samfreq))
    k <- (k * N)/samfreq
  if(max(k) > N) {
    print("number of k must be less than N")
    stop()
  }
  if(any((p > pi) | (p <  - pi))) {
    print("p must be within plus or minus pi")
    stop()
  }
  
  mat <- rep(0, N)
  t <- 1:N
  k <- k + 1
  if(is.null(const))
    const <- rep(0, length(A))
  for(j in 1:length(A)) {
    if(is.null(expon))
      expon <- rep(1, N)
    else expon <- expon[j]^(c(0:(N - 1)))
    theta <- (2 * pi * (k[j] - 1) * (t - 1))/N
    wavef <- const[j] + A[j] * expon * cos(theta + p[j])
    mat <- mat + wavef
  }
  if(plotf) {
    if(is.null(ylim))
      ylim <- range(mat)
    if(is.null(xlim))
      xlim <- c(0, (N-1))
    graphics::plot(c(0:(N - 1)), mat, xlab = xlab, ylab = ylab, type = type, 
         ylim = ylim, xlim = xlim, ...)
    if(dopoints)
      graphics::points(c(0:(N - 1)), mat, pch = 16, mkh = 0.05)
    
  }
  if(values)
    mat
  else invisible()
}











##' Function to plot a digital sinusoid and the circle from which it is
##' derived.
##' 
##' A digital sinusoid is derived the movement of a point around a circle.  The
##' function shows the relationship between the two for various parameter
##' settings.
##' 
##' 
##' @param A Amplitude of the circle/sinusoid.
##' @param k Frequency of the sinusoid
##' @param p Phase of the sinusoid
##' @param N Number of points per cycle or revolution.
##' @param const A constant corresponding to k + A*cos(2*pi*k+p)
##' @param figsize Set the figure size as pin <- c(figsize, figsize/2).
##' Defaults to figsize = 8.
##' @param npoints The number of points used in plotting the circle. Defaults
##' to 500
##' @param col An integer for the color in plotting the sinusoid and points
##' around the circle
##' @param cplot Now redundant
##' @param splot Now redundant
##' @param numplot Logical. If T (defaults), the digital points around the
##' circle are numbered
##' @param axes Logical. If T, plot axes.
##' @param incircle Logical. If T, plot an the angle between digital points in
##' the circle.
##' @param arrow Logical. If T, plot an arrow on incircle showing the direction
##' of movement.
##' @param linetype Specify a linetype. Same function as lty in plot
##' @param textplot A list containing \$radius, \$textin, \$pivals for plotting
##' text at specified angles and radii on the circle. \$radius: a vector of
##' amplitudes of the radii at which the text is to be plotted; \$textin: a
##' vector of chacacter labels to be plotted; \$pivals: the angle, in radians
##' relative to zero radians (top of the circle) at which the text is to be
##' plotted. Defaults to NULL
##' @param lineplot Plot lines from the centre of the circle to the
##' circumference. lineplot is a vector specifying the angle in radians (zero
##' corresponds to the top of the circle)
##' @param ylab Specify a y-axis label.
##' @param super Superimpose a part solid circle and corresponding sinusoid.
##' This needs to be a list containing \$first and \$last, which are values
##' between 0 and 2*pi defining the beginning and ending of the part circle
##' which is to be superimposed
##' @param xaxlab Now redundant
##' @param xlab Specify an x-axis label.
##' @param type Specify a type.
##' @param fconst A single elment numeric vector for the aspect ratio in a
##' postscript plot. Defaults to 3.5/3.1 which is appropriate for a postscript
##' setting of setps(h=4, w=4)
##' @param pointconst The radius for plotting the numbers around the circle.
##' Defaults to 1.2 * A
##' @author Jonathan Harrington
##' @seealso \code{\link{cr}}
##' @references Harrington, J, & Cassidy, S. 1999. Techniques in Speech
##' Acoustics. Kluwer
##' @keywords dplot
##' @examples
##' 
##' crplot()
##' # sine wave
##' crplot(p=-pi/2)
##' 
##' crplot(k=3)
##' 
##' # aliasing
##' crplot(k=15)
##' 
##' @export crplot
"crplot" <- function(A = 1, k = 1, p = 0, N = 16, const = NULL, figsize = 8, 
                     npoints = 500, col = 1, cplot = TRUE, splot = TRUE, numplot = TRUE, axes = TRUE, 
                     incircle = TRUE, arrow = TRUE, linetype = 1, textplot = NULL, lineplot = NULL,
                     ylab = "Amplitude", super = NULL, xaxlab = NULL, type = "b", 
                     xlab = "Time (number of points)", fconst = 3.5/3.1, pointconst = 1.2)
{
  
  if(A < -2 | A > 2)
    stop("choose A to be between plus or minus 2")
  if(N > 50)
    numplot = FALSE
  
  "cr.lines" <- function(pivals, Ax = 1, Ay = 1, 
                         plotshift, col = 1, lty = 2)
  {
    # called from within crplot() only
    c.x <-  - plotshift
    c.y <- 0	# add lines to a circle plot extending from the centre of the
    # circle, to the circumference at pivals, where pivals is in radians
    for(j in 1:length(pivals)) {
      x.th1 <- c(c.x, Ax *  - sin(pivals[j]) - plotshift)
      y.th1 <- c(c.y, Ay * cos(pivals[j]))
      graphics::lines(x.th1, y.th1, col = col, lty = lty)
    }
  }
  
  "cr.super" <- function(first = 0, last = pi/3, Ax = 1, Ay = 1, 
                         k = 1, p = 0, const = NULL, figsize = 8, 
                         npoints = 500, N = npoints, col = 1, 
                         axes = TRUE, ylab = "amplitude", lwd = 2)
  {
    first <- round(1 + ((first * (npoints - 1))/(2 * pi)))
    last <- round(1 + ((last * (npoints - 1))/(2 * pi)))
    A <- (A * figsize)/8
    pin <- c(figsize, figsize/2)
    plotshift <- 0.5 * pin[1] * 0.625
    nums <- seq(0 + p, (2 * pi) + p, length = npoints)
    cosv <- Ay * cos(nums)
    sinv <- Ax *  - sin(nums)
    graphics::plot(sinv[first:last] - plotshift, cosv[first:last], type = "l", xlim
         = c( - pin[1]/2, pin[1]/2), ylim = c( - pin[2]/2, pin[2]/2), 
         axes = FALSE, ylab = "", xlab = "", lwd = lwd)
    if(k == 0)
      theta <- rep(0 + p, N)
    else theta <- seq(0 + p, k * 2 * pi + p, by = (k * 2 * 
                                                     pi)/N)[1:N]
    sinpoints <- Ax *  - sin(theta) - plotshift
    cospoints <- Ay * cos(theta)
    xvals <- seq(0, figsize/2, length = N)
    graphics::par(new = TRUE)
    graphics::plot(xvals[first:last], cospoints[first:last], type = "l", xlim = c( - 
                                                                           pin[1]/2, pin[1]/2), ylim = c( - pin[2]/2, pin[2]/2), axes = FALSE, 
         ylab = "", xlab = "", col = col, lwd = lwd)
  }
  
  "cr.text" <- function(Ax, Ay, radius, textin, pivals, plotshift, col = 1)
  {
    ## called by Cr(); plots text at specified points around
    ## the circle. See Cr() for documentation
    for(j in 1:length(textin)) {
      numsin <- radius[j] * Ax *  - sin(pivals[j]) - plotshift
      numcos <- radius[j] * Ay * cos(pivals[j])
      graphics::text(numsin, numcos, textin[j], col = col)
    }
  }
  
  
  ## draw a circle and the corresponding sinusoid
  ## A, amplitude i.e. the radius
  ## k: the frequency
  ## p: a p of zero corresponds to the top of the circle
  ## N: number of datapoints
  ## const: a constant corresponding to k + A*cos(2*pi*k+p)
  ## figsize in area. 8 corresponds to 4 (width) x 2 (height) inches
  ## npoints: no. of points used to plot the circle
  ## col: colour. defaults to 1
  ## cplot: do you want to plot both the circle ? default is T
  ## numplot: plot the numbers on the circle? default is T
  ## splot: do you want to plot the sinusoid? default is T
  ## plot the axes? defaults to T
  ## incircle: plot the inner circle showing the angle between points?
  ## arrow: plot an arrow on the part inner circle? defaults to T
  ## linetype: linetype for plotting the sinusoid. Defaults to a solid line
  ## textplot: a list containing $radius, $textin, $pivals
  ## for plotting text at specified angles and radii on
  ## the circle. $radius: a vector of  amplitudes of the radii at
  ## which the text is to be plotted; $textin: a vector
  ## of chacacter labels to be plotted; $pivals: the angle, in radians
  ## relative to  zero radians (top of the circle) at which
  ## the text is to be plotted. Defaults to NULL
  ## lineplot: plot lines from the centre of the circle
  ## to the circumference. lineplot should be a vector specifying
  ## the angle in radians (zero corresponds to the top of the circle)
  ## super: superimpose a part solid circle and corresponding
  ## sinusoid. This needs to be a list containing $first and
  ## $last, which are values between 0 and 2*pi defining
  ## the beginning and ending of the part circle which is
  ## to be superimposed
  ## xaxlab: a character vector. Add
  ## your own axis labels
  ## the sinusoid plot. The characters are placed at equal
  ## intervals from the beginning to the end of the sinusoid.
  ## fconst: this is to get the aspect ratio correct for 
  ## postscript using setps(h=4, w=4)
  ## pointconst: the radius of numbers around the circle
  if(k > N) {
    print("number of k must be less than N")
    stop()
  }
  if(p > pi | p <  - pi) {
    print("p must be within plus or minus pi")
    stop()
  }
  ## different scale factors for x and y to frob aspect ratio
  Ay <- (A * figsize)/8
  Ax <- fconst * Ay# par(pin = c(figsize, figsize/2))
  pin <- c(figsize, figsize/2)
  plotshift <- 0.5 * pin[1] * 0.625## .
  ## plot the circle
  nums <- seq(0 + p, (2 * pi) + p, length = npoints)
  cosv <- Ay * cos(nums)
  sinv <- Ax *  - sin(nums)
  if(cplot) {
    graphics::plot(sinv - plotshift, cosv, type = "l", xlim = c( - pin[1]/2, 
                                                       pin[1]/2), ylim = c( - pin[2]/2, pin[2]/2), axes = FALSE, 
         ylab = ylab, xlab = "", col = col, lty = linetype)
  }
  ## plot the points on the circle
  if(k == 0)
    theta <- rep(0 + p, N)
  else theta <- seq(0 + p, k * 2 * pi + p, by = (k * 2 * 
                                                   pi)/N)[1:N]
  sinpoints <- (Ax *  - sin(theta)) - plotshift
  cospoints <- Ay * cos(theta)
  if(cplot) {
    if(numplot) {
      graphics::points(cbind(sinpoints, cospoints), pch = 16, mkh = 
               0.05)
      ## plot the numbers around the circle with an extended radius
      numvals <- round(theta %% (2 * pi), 2)
      for(j in unique(numvals)) {
        temp <- numvals == j
        numsin <- pointconst * Ax *  - sin(theta[temp][
          1]) - plotshift
        numcos <- pointconst * Ay * cos(theta[temp][1])
        graphics::text(numsin, numcos, paste(c(0:(N - 1))[temp], 
                                   collapse = " "), cex = 1, col = col)
      }
    }
    if(incircle) {
      ## show the angle, theta, as radii from the first to the second point
      if(k != 0 & k != N) {
        x.th1 <- c( - plotshift, Ax *  - sin(theta[1]) - 
                      plotshift)
        y.th1 <- c(0, Ay * cos(theta[1]))
        x.th2 <- c( - plotshift, A *  - sin(theta[2]) - 
                      plotshift)
        y.th2 <- c(0, Ay * cos(theta[2]))
        graphics::lines(x.th1, y.th1, col = col)
        graphics::lines(x.th2, y.th2, col = col)## .
        ## draw a part circle between these lines
        cir.thet <- seq(theta[1], theta[2], length = 
                          round((npoints * (theta[2] - theta[1]))/(2 * 
                                                                     pi)))
        graphics::par(new = TRUE)
        cos.in <- Ay * cos(cir.thet) * 0.5
        sin.in <- Ax *  - sin(cir.thet) * 0.5
        graphics::plot(sin.in - plotshift, cos.in, type = "l", 
             xlim = c( - pin[1]/2, pin[1]/2), ylim = c( - 
                                                          pin[2]/2, pin[2]/2), axes = FALSE, ylab = "", 
             xlab = "", col = col)## .
        len.in <- length(cos.in)
        if(arrow)
          graphics::arrows(sin.in[len.in - 1] - plotshift, cos.in[
            len.in - 1], sin.in[len.in] - plotshift, 
            cos.in[len.in], col = col, code=2, length=.1)
      }
    }
  }
  ## .
  ## plot the cosine wave
  if(splot) {
    xvals <- seq(0, figsize/2, length = N) * fconst - 0.2
    graphics::par(new = TRUE)
    graphics::plot(xvals, cospoints, type = "l", xlim = c( - pin[1]/2, pin[1]/
                                                   2), ylim = c( - pin[2]/2, pin[2]/2), axes = FALSE, ylab = 
           "", xlab = "", col = col, lty = linetype)
    graphics::points(xvals, cospoints, pch = 16, mkh = 0.05)
  }
  ## add any text if specified
  if(!is.null(textplot)) cr.text(Ax, Ay, textplot$radius, textplot$textin,
                                 textplot$pivals, plotshift, col = col)
  ## add any lines from the centre to the circumference
  if(!is.null(lineplot))
    cr.lines(lineplot, Ax, Ay, plotshift, col = col)
  if(axes) {
    if(is.null(xaxlab))
      graphics::axis(side = 1, line =  - A * 1.1, at = xvals, labels = 
             c(0:(N - 1)))
    else graphics::axis(side = 1, line =  - A * 1.1, at = seq(0, figsize/2, 
                                                    length = length(xaxlab)) - 0.2, labels = xaxlab
    )
    graphics::mtext(xlab, at = figsize/4, line = A * 2, side = 1)
    if(is.null(const))
      graphics::axis(side = 2)
    else graphics::axis(side = 2, at = seq( - figsize/4, figsize/4, length = 
                                    5), labels = seq( - figsize/4, figsize/4, 
                                                      length = 5) + const)
    graphics::abline(h = 0, lty = 2)
  }
  ## superimpose a part circle
  if(!is.null(super)) {
    graphics::par(new = TRUE)
    cr.super(first = super$first, last = super$last, Ax = Ax, Ay = 
               Ay, k = k, p = p, figsize = figsize, 
             const = const, npoints = npoints)
  }
}
