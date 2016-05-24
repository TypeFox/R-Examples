#' Plot function for an RLum.Results S4 class object
#'
#' The function provides a standardised plot output for data of an RLum.Results
#' S4 class object
#'
#' The function produces a multiple plot output.  A file output is recommended
#' (e.g., \code{\link{pdf}}).
#'
#' @param object \code{\linkS4class{RLum.Results}} (\bold{required}): S4 object
#' of class \code{RLum.Results}
#'
#' @param single \code{\link{logical}} (with default): single plot output
#' (\code{TRUE/FALSE}) to allow for plotting the results in as few plot windows
#' as possible.
#'
#' @param \dots further arguments and graphical parameters will be passed to
#' the \code{plot} function.
#'
#' @return Returns multiple plots.
#'
#' @note Not all arguments available for \code{\link{plot}} will be passed!
#' Only plotting of \code{RLum.Results} objects are supported.
#'
#' @section Function version: 0.2.1
#'
#' @author Christoph Burow, University of Cologne (Germany), Sebastian Kreutzer, IRAMAT-CRP2A,
#' Universite Bordeaux Montaigne (France)
#'
#' @seealso \code{\link{plot}}, \code{\link{plot_RLum}},
#'
#' @references #
#'
#' @keywords aplot
#'
#' @examples
#'
#'
#' ###load data
#' data(ExampleData.DeValues, envir = environment())
#'
#' # apply the un-logged minimum age model
#' mam <- calc_MinDose(data = ExampleData.DeValues$CA1, sigmab = 0.2, log = TRUE, plot = FALSE)
#'
#' ##plot
#' plot_RLum.Results(mam)
#'
#' # estimate the number of grains on an aliquot
#' grains<- calc_AliquotSize(grain.size = c(100,150), sample.diameter = 1, plot = FALSE)
#'
#' ##plot
#' plot_RLum.Results(grains)
#'
#'
#' @export
plot_RLum.Results<- function(
  object,
  single = TRUE,
  ...
){

  ##============================================================================##
  ## CONSISTENCY CHECK OF INPUT DATA
  ##============================================================================##

  ##check if object is of class RLum.Data.Curve
  if(!is(object,"RLum.Results")){
    stop("[plot_RLum.Results()] Input object is not of type 'RLum.Results'")
  }

  ##============================================================================##
  ## SAFE AND RESTORE PLOT PARAMETERS ON EXIT
  ##============================================================================##
  par.old <- par(no.readonly = TRUE)
  on.exit(par(par.old))

  ##============================================================================##
  ## ... ARGUMENTS
  ##============================================================================##

  ##deal with addition arguments
  extraArgs <- list(...)

  ##main
  main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
  {""}
  ##mtext
  mtext <- if("mtext" %in% names(extraArgs)) {extraArgs$mtext} else
  {""}
  ##log
  log <- if("log" %in% names(extraArgs)) {extraArgs$log} else
  {""}
  ##lwd
  lwd <- if("lwd" %in% names(extraArgs)) {extraArgs$lwd} else
  {1}
  ##lty
  lty <- if("lty" %in% names(extraArgs)) {extraArgs$lty} else
  {1}
  ##type
  type <- if("type" %in% names(extraArgs)) {extraArgs$type} else
  {"l"}
  ##pch
  pch <- if("pch" %in% names(extraArgs)) {extraArgs$pch} else
  {1}
  ##col
  col <- if("col" %in% names(extraArgs)) {extraArgs$col} else
  {"black"}

  ##============================================================================##
  ## PLOTTING
  ##============================================================================##

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## CASE 1: Minimum Age Model / Maximum Age Model
  if(object@originator=="calc_MinDose" || object@originator=="calc_MaxDose") {

    ## single MAM estimate
    # plot profile log likelhood
    tryCatch({
      suppressWarnings(
        bbmle::plot(object@data$profile, show.points=FALSE, plot.confstr=TRUE, onepage = single, ask = FALSE)
      )
    }, error = function(e) {
      if (single)
        par(mfrow=c(2, 2))
      param <- c("gamma", "sigma", "p0", "mu")
      for (i in param) {
        if (object@data$summary$par == 3 && i == "mu")
          break
        tryCatch({
          bbmle::plot(object@data$profile, which = i)
        }, error = function(e)  {
          message(paste("Unable to plot the Likelihood profile for:", i))
        })
      }
      par(mfrow=c(1,1))
    })

    ## bootstrap MAM estimates
    if(object@data$args$bootstrap==TRUE) {

      # save previous plot parameter and set new ones
      .pardefault<- par(no.readonly = TRUE)

      # get De-llik pairs
      pairs<- object@data$bootstrap$pairs$gamma

      # get polynomial fit objects
      poly.lines<- list(poly.three=object@data$bootstrap$poly.fits$poly.three,
                        poly.four=object@data$bootstrap$poly.fits$poly.four,
                        poly.five=object@data$bootstrap$poly.fits$poly.five,
                        poly.six=object@data$bootstrap$poly.fits$poly.six)

      # define polynomial curve functions for plotting
      poly.curves<- list(poly.three.curve=function(x) { poly.lines$poly.three$coefficient[4]*x^3 + poly.lines$poly.three$coefficient[3]*x^2 + poly.lines$poly.three$coefficient[2]*x + poly.lines$poly.three$coefficient[1] },
                         poly.four.curve=function(x) { poly.lines$poly.four$coefficient[5]*x^4 + poly.lines$poly.four$coefficient[4]*x^3 + poly.lines$poly.four$coefficient[3]*x^2 + poly.lines$poly.four$coefficient[2]*x + poly.lines$poly.four$coefficient[1] },
                         poly.five.curve=function(x) { poly.lines$poly.five$coefficient[6]*x^5 + poly.lines$poly.five$coefficient[5]*x^4 + poly.lines$poly.five$coefficient[4]*x^3 + poly.lines$poly.five$coefficient[3]*x^2 + poly.lines$poly.five$coefficient[2]*x + poly.lines$poly.five$coefficient[1] },
                         poly.six.curve=function(x) { poly.lines$poly.six$coefficient[7]*x^6 + poly.lines$poly.six$coefficient[6]*x^5 + poly.lines$poly.six$coefficient[5]*x^4 + poly.lines$poly.six$coefficient[4]*x^3 + poly.lines$poly.six$coefficient[3]*x^2 + poly.lines$poly.six$coefficient[2]*x + poly.lines$poly.six$coefficient[1] })

      ## --------- PLOT "RECYCLE" BOOTSTRAP RESULTS ------------ ##

      if(single==TRUE) {
        layout(cbind(c(1,1,2, 5,5,6), c(3,3,4, 7,7,8)))
        par(cex = 0.6)
      } else {
        layout(matrix(c(1,1,2)),2,1)
        par(cex = 0.8)
      }

      for(i in 1:4) {
        ## ----- LIKELIHOODS

        # set margins (bottom, left, top, right)
        par(mar=c(0,5,5,3))

        # sort De and likelihoods by De (increasing)
        pairs<- pairs[order(pairs[,1]),]

        # remove invalid NA values
        pairs<- na.omit(pairs)

        plot(x=pairs[,1],
             y=pairs[,2],
             xlab="Equivalent Dose [Gy]",
             ylab="Likelihood",
             xlim=range(pretty(pairs[,1])),
             ylim=range(pretty(c(0, as.double(quantile(pairs[,2],probs=0.98))))),
             xaxt = "n",
             xaxs = "i",
             yaxs = "i",
             bty = "l",
             main="Recycled bootstrap MAM-3")

        axis(side = 1, labels = FALSE, tick = FALSE)

        # add subtitle
        mtext(as.expression(bquote(italic(M) == .(object@data$args$bs.M) ~ "|" ~
                                     italic(N) == .(object@data$args$bs.N) ~ "|" ~
                                     italic(sigma[b])  == .(object@data$args$sigmab) ~
                                     "\u00B1" ~ .(object@data$args$sigmab.sd) ~ "|" ~
                                     italic(h) == .(round(object@data$args$bs.h,1))
        )
        ),
        side = 3, line = 0.3, adj = 0.5,
        cex = if(single){0.5}else{0.8})

        # add points
        points(x=pairs[,1], y=pairs[,2], pch=1, col = "grey80")

        # get polynomial function
        poly.curve<- poly.curves[[i]]

        # add curve to plot
        curve(poly.curve, from = min(pairs[,1]), to = (max(pairs[,1])),
              col = "black", add = TRUE, type = "l")

        # add legend
        legend<- c("Third degree", "Fourth degree", "Fifth degree", "Sixth degree")
        legend("topright",  xjust = 0,
               legend = legend[i],
               y.intersp = 1.2,
               bty = "n",
               title = "Polynomial Fit",
               lty = 1,
               lwd= 1.5)

        ## ----- RESIDUALS

        # set margins (bottom, left, top, right)
        par(mar=c(5,5,0,3))

        plot(x = pairs[,1],
             y = residuals(poly.lines[[i]]),
             ylim = c(min(residuals(poly.lines[[i]]))*1.2,
                      as.double(quantile(residuals(poly.lines[[i]]),probs=0.99))),
             xlim=range(pretty(pairs[,1])),
             xaxt = "n",
             bty = "l",
             xaxs = "i",
             col = "grey80",
             ylab = "Fit residual",
             xlab = "Equivalent dose [Gy]")

        axis(side = 1, labels = TRUE, tick = TRUE)

        # add horizontal line
        abline(h = 0, lty=2)

        # calculate residual sum of squares (RSS) and add to plot
        rss<- sum(residuals(poly.lines[[i]])^2)
        mtext(text = paste("RSS =",round(rss,3)), adj = 1,
              side = 3, line = -2,
              cex = if(single){0.6}else{0.8})

        ## ----- PROPORTIONS

      }##EndOf::Plot_loop

      # restore previous plot parameters
      par(.pardefault)

      ### TODO: plotting of the LOESS fit needs to be fleshed out
      ### possibly integrate this in the prior polynomial plot loop

      ### LOESS PLOT
      pairs<- object@data$bootstrap$pairs$gamma
      pred<- predict(object@data$bootstrap$loess.fit)
      loess<- cbind(pairs[,1], pred)
      loess<- loess[order(loess[,1]),]

      # plot gamma-llik pairs
      plot(pairs,
           ylim = c(0, as.double(quantile( pairs[,2],probs=0.99))),
           ylab = "Likelihood",
           xlab = "Equivalent dose [Gy]",
           col = "gray80")

      # add LOESS line
      lines(loess, type = "l", col = "black")

      ### ------ PLOT BOOTSTRAP LIKELIHOOD FIT

      par(mar=c(5,4,4,4))

      xlim<- range(pretty(object@data$data[,1]))
      xlim[1]<- xlim[1]-object@data$data[which.min(object@data$data[,1]),2]
      xlim[2]<- xlim[2]+object@data$data[which.max(object@data$data[,1]),2]
      xlim<- range(pretty(xlim))

      # empty plot
      plot(NA,NA,
           xlim=xlim,
           ylim=c(0,2),
           xlab="Equivalent dose [Gy]",
           ylab="",
           bty="l",
           axes=FALSE,
           xaxs="i",
           yaxs="i",
           yaxt="n")

      axis(side = 1)
      axis(side = 2, at = c(0,0.5,1))

      mtext(text = "Normalised likelihood / density", side = 2, line = 2.5, adj = 0)

      # set the polynomial to plot
      poly.curve<- poly.curves[[1]] # three degree poly

      # plot a nice grey polygon as in the publication
      step<- 0.1
      x<- seq(min(pairs[,1]), max(pairs[,1]), step)
      y<- poly.curve(x)
      # normalise y-values
      y<- y/max(y)

      x<- c(min(pairs[,1]), x, max(pairs[,1]))
      y<- c(0, y, 0)

      # cutoff negative y values
      neg<- which(y<0)
      y<- y[-neg]
      x<- x[-neg]

      # add bootstrap likelihood polygon to plot
      polygon(x, y, col = "grey80", border = NA)

      ### ----- PLOT MAM SINGLE ESTIMATE

      # symmetric errors, might not be appropriate
      mean<- object@data$summary$de
      sd<- object@data$summary$de_err

      x<- seq(mean-5*sd, mean+5*sd, 0.001)
      y<- dnorm(seq(mean-5*sd, mean+5*sd, 0.001), mean, sd)
      # normalise y-values
      y<- y/max(y)

      points(x, y,
             type="l",
             col="red")

      ## asymmetric errors
      x<- unlist(object@data$profile@profile$gamma$par.vals[,1])
      y<- abs(unlist(object@data$profile@profile$gamma$z))

      if(object@data$args$log == TRUE) {
        x<- exp(x)
      }

      # now invert the data by shifting
      y<- -y
      y<- y-min(y)
      y<- y/max(y)

      # fit a smoothing spline
      l<- spline(x = x, y = y, method = "n", n = 1000)

      # make the endpoints zero
      l$y[1]<- l$y[length(l$y)]<- 0

      # add profile log likelihood curve to plot
      lines(l, col="blue", lwd=1)

      # add vertical lines of the mean values
      #points(x = 80, y = 100,type = "l")

      #### ------ PLOT DE
      par(new = TRUE)

      # sort the data in ascending order
      dat<- object@data$data[order(object@data$data[,1]),]

      x<- dat[,1]
      y<- 1:length(object@data$data[,1])

      plot(x = x, y = y,
           xlim=xlim,
           ylim=c(0, max(y)+1),
           axes = FALSE,
           pch = 16,
           xlab = "",
           ylab="",
           xaxs="i",
           yaxs="i")

      axis(side = 4)
      mtext(text = "# Grain / aliquot", side = 4, line = 2.5)

      # get sorted errors
      err<- object@data$data[order(object@data$data[,1]),2]

      # fancy error bars
      arrows(x0 = x-err, y0 = y,
             x1 =  x+err, y1 = y,
             code = 3, angle = 90, length = 0.05)

      ### ---- AUXILLARY

      # add legend
      legend("bottomright",
             bty = "n",
             col = c("grey80", "red", "blue", "black"),
             pch = c(NA,NA,NA,16),
             lty = c(1,1,1,1),
             lwd=c(10,2,2,2),
             legend = c("Bootstrap likelihood", "Profile likelihood (gaussian fit)","Profile likelihood", "Grain / aliquot"),
      )

    }##EndOf::Bootstrap_plotting
  }#EndOf::CASE1_MinimumAgeModel-3


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## CASE 2: Central Age Model
  if(object@originator=="calc_CentralDose") {

    # get profile log likelihood data
    sig<- object@data$profile$sig*100
    llik<- object@data$profile$llik

    # save previous plot parameter and set new ones
    .pardefault<- par(no.readonly = TRUE)

    # plot the profile log likeihood
    par(oma=c(2,1,2,1),las=1,cex.axis=1.2, cex.lab=1.2)
    plot(sig,llik,type="l",xlab=as.expression(bquote(sigma[OD]~"[%]")),ylab="Log likelihood",lwd=1.5)
    abline(h=0,lty=3)
    abline(h=-1.92,lty=3)
    title(as.expression(bquote("Profile log likelihood for" ~ sigma[OD])))

    # find upper and lower confidence limits for sigma
    sigmax<- sig[which.max(llik)]
    tf<- abs(llik+1.92) < 0.05
    sig95<- sig[tf]
    ntf<- length(sig95)
    sigL<- sig95[1]
    sigU<- sig95[ntf]

    # put them on the graph
    abline(v=sigL)
    abline(v=sigmax)
    abline(v=sigU)
    dx<- 0.006
    dy<- 0.2
    ytext<- min(llik) + dy
    res<- c(sigL,sigmax,sigU)
    text(res+dx,rep(ytext,3),round(res,2),adj=0)

    # restore previous plot parameters
    par(.pardefault)
    rm(.pardefault)
  }##EndOf::Case 2 - calc_CentralDose()


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## CASE 3: Fuchs & Lang 2001
  if(object@originator=="calc_FuchsLang2001") {

    ##deal with addition arguments
    extraArgs <- list(...)

    main <- if("main" %in% names(extraArgs)) {extraArgs$main} else {"Fuchs & Lang (2001)"}
    xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else {expression(paste(D[e]," [s]"))}
    ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else {"# Aliquots"}
    sub <-  if("sub" %in% names(extraArgs)) {extraArgs$sub} else {""}
    cex <- if("cex" %in% names(extraArgs)) {extraArgs$cex} else {1}
    lwd <- if("lwd" %in% names(extraArgs)) {extraArgs$lwd} else {1}
    pch <- if("pch" %in% names(extraArgs)) {extraArgs$pch} else {19}
    ylim <- if("ylim" %in% names(extraArgs)) {extraArgs$ylim} else {c(1,length(object@data$data[,1])+3)}
    xlim <- if("xlim" %in% names(extraArgs)) {extraArgs$xlim} else {c(min(object@data$data[,1])-max(object@data$data[,2]), max(object@data$data[,1])+max(object@data$data[,2]))}
    mtext <- if("mtext" %in% names(extraArgs)) {extraArgs$mtext} else {"unknown sample"}

    # extract relevant plotting parameters
    o<- order(object@data$data[1])
    data_ordered<- object@data$data[o,]
    usedDeValues<- object@data$usedDeValues
    n.usedDeValues<- object@data$summary$n.usedDeValues

    par(cex = cex, mfrow=c(1,1))

    ##PLOT
    counter<-seq(1,max(o))

    plot(NA,NA,
         ylim = ylim,
         xlim = xlim,
         xlab = xlab,
         ylab = ylab,
         main = main,
         sub = sub)

    ##SEGMENTS
    segments(data_ordered[,1]-data_ordered[,2],1:length(data_ordered[,1]),
             data_ordered[,1]+data_ordered[,2],1:length(data_ordered[,1]),
             col="gray")


    ##POINTS
    points(data_ordered[,1], counter,pch=pch)

    ##LINES
    ##BOUNDARY INFORMATION
    ##lower boundary
    lines(c(
      usedDeValues[length(usedDeValues[,1])-n.usedDeValues+1,1], #boundary_counter for incorporate skipped values
      usedDeValues[length(usedDeValues[,1])-n.usedDeValues+1,1]),
      c(min(o)-0.5,max(o)+0.5),
      col="red",
      lty="dashed", lwd = lwd)


    #upper boundary
    lines(c(max(usedDeValues[,1]),max(usedDeValues[,1])),c(min(o)-0.5,max(o)+0.5),
          col="red",lty="dashed", lwd = lwd)

    #plot some further informations into the grafik
    arrows(
      usedDeValues[length(usedDeValues[,1])-n.usedDeValues+1,1]+usedDeValues[length(usedDeValues[,1])-n.usedDeValues+1,1]*0.02, #x1
      max(o)+0.5, #y1
      max(usedDeValues[,1]-usedDeValues[,1]*0.02), #x2
      max(o)+0.5, #y2,
      code=3,
      length=0.03)

    text(
      c(
        usedDeValues[length(usedDeValues[,1])-n.usedDeValues+1,1],
        usedDeValues[length(usedDeValues[,1])-n.usedDeValues+1,1]),
      c(max(o)+2,max(o)+2),
      labels=paste("used values = ", n.usedDeValues),
      cex=0.6*cex,
      adj=0)

    ##MTEXT
    mtext(side=3,mtext,cex=cex)
  }


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## CASE 4: Finite Mixture Model

  if(object@originator == "calc_FiniteMixture") {
    if(length(object@data$args$n.components) > 1L) {

      ##deal with addition arguments
      extraArgs <- list(...)

      main <- if("main" %in% names(extraArgs)) {extraArgs$main} else {"Finite Mixture Model"}
      plot.proportions<- if("plot.proportions" %in% names(extraArgs)) {extraArgs$plot.proportions} else {TRUE}
      pdf.colors<- if("pdf.colors" %in% names(extraArgs)) {extraArgs$pdf.colors} else {"gray"}
      pdf.weight<- if("pdf.weight" %in% names(extraArgs)) {extraArgs$pdf.weight} else {TRUE}
      pdf.sigma<- if("pdf.sigma" %in% names(extraArgs)) {extraArgs$pdf.sigma} else {"sigmab"}

      # extract relevant data from object
      n.components<- object@data$args$n.components
      comp.n<- object@data$components
      sigmab<- object@data$args$sigmab
      BIC.n<- object@data$BIC$BIC
      LLIK.n<- object@data$llik$llik

      # save previous plot parameter and set new ones
      .pardefault<- par(no.readonly = TRUE)

      ## DEVICE AND PLOT LAYOUT
      n.plots<- length(n.components) #number of PDF plots in plotarea #1
      seq.vertical.plots<- seq(from = 1, to = n.plots, by = 1) #indices
      ID.plot.two<- n.plots+if(plot.proportions==TRUE){1}else{0} #ID of second plot area
      ID.plot.three<- n.plots+if(plot.proportions==TRUE){2}else{1} #ID of third plot area

      #empty vector for plot indices
      seq.matrix<- vector(mode="integer", length=4*n.plots)

      #fill vector with plot indices in correct order
      cnt<- 1
      seq<- seq(1,length(seq.matrix),4)
      for(i in seq) {
        seq.matrix[i]<- cnt
        seq.matrix[i+1]<- cnt
        seq.matrix[i+2]<- if(plot.proportions==TRUE){ID.plot.two}else{cnt}
        seq.matrix[i+3]<- ID.plot.three

        cnt<- cnt+1
      }

      # create device layout
      layout(matrix(c(seq.matrix),4,n.plots))

      # outer margins (bottom, left, top, right)
      par(oma=c(2.5,5,3,5))

      # general plot parameters (global scaling, allow overplotting)
      par(cex = 0.8, xpd = NA)

      # define color palette for prettier output
      if(pdf.colors == "colors") {
        col.n<- c("red3", "slateblue3", "seagreen", "tan3", "yellow3",
                  "burlywood4", "magenta4", "mediumpurple3", "brown4","grey",
                  "aquamarine")
        poly.border<- FALSE
      }
      if(pdf.colors == "gray" || pdf.colors == "grey") {
        col.n<- gray.colors(length(n.components)*2)
        poly.border<- FALSE
      }
      if(pdf.colors == "none") {
        col.n<- NULL
        poly.border<- TRUE
      }

      ##--------------------------------------------------------------------------
      ## PLOT 1: EQUIVALENT DOSES OF COMPONENTS

      ## create empty plot without x-axis
      for(i in 1:n.plots) {

        pos.n<- seq(from = 1, to = n.components[i]*3, by = 3)

        # set margins (bottom, left, top, right)
        par(mar=c(1,0,2,0))

        # empty plot area
        plot(NA, NA,
             xlim=c(min(n.components)-0.2, max(n.components)+0.2),
             ylim=c(min(comp.n[pos.n,]-comp.n[pos.n+1,], na.rm = TRUE),
                    max((comp.n[pos.n,]+comp.n[pos.n+1,])*1.1, na.rm = TRUE)),
             ylab="",
             xaxt="n",
             yaxt="n",
             xlab="")

        # add text in upper part of the plot ("k = 1,2..n")
        mtext(bquote(italic(k) == .(n.components[i])),
              side = 3, line = -2, cex=0.8)

        # add y-axis label (only for the first plot)
        if(i==1) {
          mtext(expression(paste("D"[e]," [Gy]")), side=2,line=2.7, cex=1)
        }

        # empty list to store normal distribution densities
        sapply.storage<- list()

        ## NORMAL DISTR. OF EACH COMPONENT
        options(warn=-1) #supress warnings for NA values

        # LOOP - iterate over number of components
        for(j in 1:max(n.components)) {

          # draw random values of the ND to check for NA values
          comp.nd.n<- sort(rnorm(n = length(object@data$data[,1]),
                                 mean = comp.n[pos.n[j],i],
                                 sd = comp.n[pos.n[j]+1,i]))

          # proceed if no NA values occured
          if(length(comp.nd.n)!=0) {

            # weight - proportion of the component
            wi<- comp.n[pos.n[j]+2,i]

            # calculate density values with(out) weights
            fooX<- function(x) {
              dnorm(x, mean = comp.n[pos.n[j],i],
                    sd = if(pdf.sigma=="se"){comp.n[pos.n[j]+1,i]}
                    else{if(pdf.sigma=="sigmab"){comp.n[pos.n[j],i]*sigmab}}
              )*
                if(pdf.weight==TRUE){wi}else{1}
            }

            # x-axis scaling - determine highest dose in first cycle
            if(i==1 && j==1){
              max.dose<- max(object@data$data[,1])+sd(object@data$data[,1])/2
              min.dose<- min(object@data$data[,1])-sd(object@data$data[,1])/2

              # density function to determine y-scaling if no weights are used
              fooY<- function(x) {
                dnorm(x, mean = na.exclude(comp.n[pos.n,]),
                      sd = na.exclude(comp.n[pos.n+1,]))
              }
              # set y-axis scaling
              dens.max<-max(sapply(0:max.dose, fooY))
            }##EndOfIf::first cycle settings


            # override y-axis scaling if weights are used
            if(pdf.weight==TRUE){
              sapply.temp<- list()
              for(b in 1:max(n.components)){

                # draw random values of the ND to check for NA values
                comp.nd.n<- sort(rnorm(n = length(object@data$data[,1]),
                                       mean = comp.n[pos.n[b],i],
                                       sd = comp.n[pos.n[b]+1,i]))

                # proceed if no NA values occured
                if(length(comp.nd.n)!=0) {

                  # weight - proportion of the component
                  wi.temp<- comp.n[pos.n[b]+2,i]

                  fooT<- function(x) {
                    dnorm(x, mean = comp.n[pos.n[b],i],
                          sd = if(pdf.sigma=="se"){comp.n[pos.n[b]+1,i]}
                          else{if(pdf.sigma=="sigmab"){comp.n[pos.n[b],i]*sigmab}}
                    )*wi.temp
                  }
                  sapply.temp[[b]]<- sapply(0:max.dose, fooT)
                }
              }
              dens.max<- max(Reduce('+', sapply.temp))
            }

            # calculate density values for 0 to maximum dose
            sapply<- sapply(0:max.dose, fooX)

            # save density values in list for sum curve of gaussians
            sapply.storage[[j]]<- sapply

            ## determine axis scaling
            # x-axis (dose)
            if("dose.scale" %in% names(extraArgs)) {
              y.scale<- extraArgs$dose.scale
            } else {
              y.scale<- c(min.dose,max.dose)
            }
            # y-axis (density)
            if("pdf.scale" %in% names(extraArgs)) {
              x.scale<- extraArgs$pdf.scale
            } else {
              x.scale<- dens.max*1.1
            }

            ## PLOT Normal Distributions
            par(new=TRUE)
            plot(sapply, 1:length(sapply)-1,
                 type="l", yaxt="n", xaxt="n", col=col.n[j], lwd=1,
                 ylim=y.scale,
                 xlim=c(0,x.scale),
                 xaxs="i", yaxs="i",
                 ann=FALSE, xpd = FALSE)

            # draw colored polygons under curve
            polygon(x=c(min(sapply), sapply,  min(sapply)),
                    y=c(0, 0:max.dose,  0),
                    col = adjustcolor(col.n[j], alpha.f = 0.66),
                    yaxt="n", border=poly.border, xpd = FALSE, lty = 2, lwd = 1.5)

          }
        }##EndOf::Component loop

        #turn warnings on again
        options(warn=0)

        # Add sum of gaussians curve
        par(new = TRUE)

        plot(Reduce('+', sapply.storage),1:length(sapply)-1,
             type="l", yaxt="n", xaxt="n", col="black",
             lwd=1.5, lty = 1,
             ylim=y.scale,
             xlim=c(0,x.scale),
             xaxs="i", yaxs="i", ann=FALSE, xpd = FALSE)

        # draw additional info during first k-cycle
        if(i == 1) {

          # plot title
          mtext("Normal distributions",
                side = 3, font = 2, line = 0, adj = 0, cex = 0.8)

          # main title
          mtext(main,
                side = 3, font = 2, line = 3.5, adj = 0.5,
                at = grconvertX(0.5, from = "ndc", to = "user"))

          # subtitle
          mtext(as.expression(bquote(italic(sigma[b]) == .(sigmab) ~
                                       "|" ~ n == .(length(object@data$data[,1])))),
                side = 3, font = 1, line = 2.2, adj = 0.5,
                at = grconvertX(0.5, from = "ndc", to = "user"), cex = 0.9)

          # x-axis label
          mtext("Density [a.u.]",
                side = 1, line = 0.5, adj = 0.5,
                at = grconvertX(0.5, from = "ndc", to = "user"))

          # draw y-axis with proper labels
          axis(side=2, labels = TRUE)
        }

        if(pdf.colors == "colors") {
          # create legend labels
          dose.lab.legend<- paste("c", 1:n.components[length(n.components)], sep="")

          if(max(n.components)>8) {
            ncol.temp<- 8
            yadj<- 1.025
          } else {
            ncol.temp<- max(n.components)
            yadj<- 0.93
          }

          # add legend
          if(i==n.plots) {
            legend(grconvertX(0.55, from = "ndc", to = "user"),
                   grconvertY(yadj, from = "ndc", to = "user"),
                   legend = dose.lab.legend,
                   col = col.n[1:max(n.components)],
                   pch = 15, adj = c(0,0.2), pt.cex=1.4,
                   bty = "n", ncol=ncol.temp, x.intersp=0.4)

            mtext("Components: ", cex = 0.8,
                  at = grconvertX(0.5, from = "ndc", to = "user"))
          }
        }

      }##EndOf::k-loop and Plot 1

      ##--------------------------------------------------------------------------
      ## PLOT 2: PROPORTION OF COMPONENTS
      if(plot.proportions==TRUE) {
        # margins for second plot
        par(mar=c(2,0,2,0))

        # create matrix with proportions from a subset of the summary matrix
        prop.matrix<- comp.n[pos.n+2,]*100

        # stacked barplot of proportions without x-axis
        barplot(prop.matrix,
                width=1,
                xlim=c(0.2, length(n.components)-0.2),
                ylim=c(0,100),
                axes=TRUE,
                space=0,
                col=col.n,
                xpd=FALSE,
                xaxt="n")

        # y-axis label
        mtext("Proportion [%]",
              side=2,line=3, cex=1)


        # add x-axis with corrected tick positions
        axis(side = 1, labels = n.components, at = n.components+0.5-n.components[1])

        # draw a box (not possible with barplot())
        box(lty=1, col="black")

        # add subtitle
        mtext("Proportion of components",
              side = 3, font = 2, line = 0, adj = 0, cex = 0.8)

      }
      ##--------------------------------------------------------------------------
      ## PLOT 3: BIC & LLIK

      # margins for third plot
      par(mar=c(2,0,2,0))

      # prepare scaling for both y-axes
      BIC.scale<- c(min(BIC.n)*if(min(BIC.n)<0){1.2}else{0.8},
                    max(BIC.n)*if(max(BIC.n)<0){0.8}else{1.2})
      LLIK.scale<- c(min(LLIK.n)*if(min(LLIK.n)<0){1.2}else{0.8},
                     max(LLIK.n)*if(max(LLIK.n)<0){0.8}else{1.2})

      # plot BIC scores
      plot(n.components, BIC.n,
           main= "",
           type="b",
           pch=22,
           cex=1.5,
           xlim=c(min(n.components)-0.2, max(n.components)+0.2),
           ylim=BIC.scale,
           xaxp=c(min(n.components), max(n.components), length(n.components)-1),
           xlab=expression(paste(italic(k), " Components")),
           ylab=expression(paste("BIC")),
           cex.lab=1.25)

      # following plot should be added to previous
      par(new = TRUE)

      # plot LLIK estimates
      plot(n.components, LLIK.n,
           xlim=c(min(n.components)-0.2, max(n.components)+0.2),
           xaxp=c(min(n.components), max(n.components), length(n.components)-1),
           ylim=LLIK.scale,
           yaxt="n", type="b", pch=16, xlab="", ylab="", lty=2, cex = 1.5)

      # subtitle
      mtext("Statistical criteria",
            side = 3, font = 2, line = 0, adj = 0, cex = 0.8)

      # second y-axis with proper scaling
      axis(side = 4, ylim=c(0,100))

      # LLIK axis label
      mtext(bquote(italic(L)[max]),
            side=4,line=3, cex=1.3)

      # legend
      legend(grconvertX(0.75, from = "nfc", to = "user"),
             grconvertY(0.96, from = "nfc", to = "user"),
             legend = c("BIC", as.expression(bquote(italic(L)[max]))),
             pch = c(22,16), pt.bg=c("white","black"),
             adj = 0, pt.cex=1.3, lty=c(1,2),
             bty = "n", horiz = TRUE, x.intersp=0.5)


      ## restore previous plot parameters
      par(.pardefault)
    }
  }##EndOf::Case 4 - Finite Mixture Model

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## CASE 5: Aliquot Size
  if(object@originator=="calc_AliquotSize") {
    if(object@data$args$MC == TRUE) {

      extraArgs <- list(...)

      main<- if("main" %in% names(extraArgs)) { extraArgs$main } else { "Monte Carlo Simulation"  }
      xlab<- if("xlab" %in% names(extraArgs)) { extraArgs$xlab } else { "Amount of grains on aliquot" }

      # extract relevant data
      MC.n<- object@data$MC$estimates
      MC.n.kde<- object@data$MC$kde
      MC.stats<- object@data$MC$statistics
      MC.q<- object@data$MC$quantile
      MC.iter<- object@data$args$MC.iter

      # set layout of plotting device
      layout(matrix(c(1,1,2)),2,1)
      par(cex = 0.8)

      ## plot MC estimate distribution

      # set margins (bottom, left, top, right)
      par(mar=c(2,5,5,3))

      # plot histogram
      hist(MC.n, freq=FALSE, col = "gray90",
           main="", xlab=xlab,
           xlim = c(min(MC.n)*0.95, max(MC.n)*1.05),
           ylim = c(0, max(MC.n.kde$y)*1.1))

      # add rugs to histogram
      rug(MC.n)

      # add KDE curve
      lines(MC.n.kde, col = "black", lwd = 1)

      # add mean, median and quantils (0.05,0.95)
      abline(v=c(MC.stats$mean, MC.stats$median, MC.q),
             lty=c(2, 4, 3,3), lwd = 1)

      # add main- and subtitle
      mtext(main, side = 3, adj = 0.5,
            line = 3, cex = 1)
      mtext(as.expression(bquote(italic(n) == .(MC.iter) ~ "|" ~
                                   italic(hat(mu)) == .(round(MC.stats$mean)) ~ "|" ~
                                   italic(hat(sigma))  == .(round(MC.stats$sd.abs)) ~ "|" ~
                                   italic(frac(hat(sigma),sqrt(n))) == .(round(MC.stats$se.abs))  ~ "|" ~
                                   italic(v) == .(round(MC.stats$skewness, 2))
      )
      ),
      side = 3, line = 0.3, adj = 0.5,
      cex = 0.9)

      # add legend
      legend("topright", legend = c("mean","median", "0.05 / 0.95 quantile"),
             lty = c(2, 4, 3), bg = "white", box.col = "white", cex = 0.9)

      ## BOXPLOT
      # set margins (bottom, left, top, right)
      par(mar=c(5,5,0,3))

      plot(NA, type="n", xlim=c(min(MC.n)*0.95, max(MC.n)*1.05),
           xlab=xlab,  ylim=c(0.5,1.5),
           xaxt="n", yaxt="n", ylab="")
      par(bty="n")
      boxplot(MC.n, horizontal = TRUE, add = TRUE, bty="n")
    }
  }#EndOf::Case 5 - calc_AliqoutSize()

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## CASE 6: calc_SourceDoseRate()
  if(object@originator=="calc_SourceDoseRate") {

    ##prepare data
    ##get data
    df <- get_RLum(object = object, data.object = "dose.rate")

    ##reduce the size for plotting, more than 100 points makes no sense
    if(nrow(df)>100) {
      df <- df[seq(1,nrow(df), length = 100),]

    }


    ##plot settings
      plot.settings <- list(
        main = "Source Dose Rate Prediction",
        xlab = "Date",
        ylab = paste0(
          "Dose rate/(",get_RLum(object = object, data.object = "parameters")$dose.rate.unit,")"),
        log = "",
        cex = 1,
        xlim = NULL,
        ylim = c(min(df[,1]) - max(df[,2]), max(df[,1]) + max(df[,2])),
        pch = 1,
        mtext = paste0(
          "source type: ", get_RLum(object = object, data.object = "parameters")$source.type,
          " | ",
          "half-life: ", get_RLum(object = object, data.object = "parameters")$halflife,
          " a"
        ),
        grid = expression(nx = 10, ny = 10),
        col = 1,
        type = "b",
        lty = 1,
        lwd = 1,
        segments = ""
      )

      ##modify list if something was set
      plot.settings <- modifyList(plot.settings, list(...))


    ##plot
      plot(
        df[,3], df[,1],
        main = plot.settings$main,
        xlab = plot.settings$xlab,
        ylab = plot.settings$ylab,
        xlim = plot.settings$xlim,
        ylim = plot.settings$ylim,
        log = plot.settings$log,
        pch = plot.settings$pch,
        col = plot.settings$pch,
        type = plot.settings$type,
        lty = plot.settings$lty,
        lwd = plot.settings$lwd
      )

      if(!is.null(plot.settings$segments)){
        segments(
          x0 = df[,3], y0 = df[,1] + df[,2],
          x1 = df[,3], y1 = df[,1] - df[,2]
        )
      }

      mtext(side = 3, plot.settings$mtext)

      if(!is.null(plot.settings$grid)){
        grid(eval(plot.settings$grid))

      }

  }#EndOf::Case 6 - calc_SourceDoseRate()

}
