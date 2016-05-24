
"plot.SimpleTable" <- function(x, estimand=c("ATE", "ATT", "ATC",
                                    "RR", "RRT", "RRC",
                                    "logRR", "logRRT", "logRRC"),
                               percent=95,
                               plot.bounds=TRUE, plot.pf=TRUE,
                               plot.sens=TRUE,  plot.prior=FALSE,
                               color.bounds="cyan",
                               color1.pf="lawngreen",
                               color2.pf="green",
                               color1.sens="magenta3",
                               color2.sens="purple4",
                               color.prior="lightgray", ymax=NULL,
                               ...){
  S <- x
  
  estimand <- match.arg(estimand)
  
  S.estimand.pf <- paste("S$", estimand, ".pf", sep="")
  S.estimand.sens <- paste("S$", estimand, ".sens", sep="")
  S.estimand.prior <- paste("S$", estimand, ".prior", sep="")

  stuff.pf <- eval(parse(text=S.estimand.pf))
  stuff.sens <- eval(parse(text=S.estimand.sens))
  stuff.prior <- eval(parse(text=S.estimand.prior))
  
  S.estimand.min <- paste("S$", estimand, ".min", sep="")
  S.estimand.max <- paste("S$", estimand, ".max", sep="")
  
  stuff.min <- eval(parse(text=S.estimand.min))
  stuff.max <- eval(parse(text=S.estimand.max))
  stuff.bounds <- paste("[", round(stuff.min, 5), ", ",
                        round(stuff.max, 5), "]", sep="")

  if (estimand %in% c("ATE", "ATT", "ATC")){
    
    x <- seq(from=-1+1e-6, to=1-1e-6, by=.0005)
    
    fit.pf <- y.pf <- hdr.pf <- fit.sens <- y.sens <- hdr.sens <- fit.prior <-
      y.prior <- hdr.prior <- NULL
    if (plot.pf){
      fit.pf <- locfit(~lp(stuff.pf, nn=.3, deg=3), renorm=TRUE,
                       ev=rbox(cut=.3), xlim=c(-1,1), maxk=500)
      y.pf <- predict(fit.pf, newdata=x)
      hdr.pf <- hdr(den=list(x=x, y=y.pf), all.modes=FALSE, prob=percent)      
    }
    if (plot.sens){
      fit.sens <- locfit(~lp(stuff.sens, nn=.3, deg=3),
                         renorm=TRUE,
                         ev=rbox(cut=.3), xlim=c(-1,1), maxk=500)
      y.sens <- predict(fit.sens, newdata=x)
      hdr.sens  <- hdr(den=list(x=x, y=y.sens), all.modes=FALSE, prob=percent)
    }
    if (plot.prior){
      fit.prior <- locfit(~lp(stuff.prior, nn=.3, deg=3), renorm=TRUE,
                          ev=rbox(cut=.3), xlim=c(-1,1), maxk=500)
      y.prior <- predict(fit.prior, newdata=x)
      hdr.prior  <- hdr(den=list(x=x, y=y.prior), all.modes=FALSE, prob=percent)
    }
    
        
    if (is.null(ymax)){
      ymax <- max(c(y.pf, y.sens, y.prior))
    }

    makeplot(x=x, y.pf=y.pf, y.sens=y.sens, y.prior=y.prior,
             hdr.pf=hdr.pf, hdr.sens=hdr.sens, hdr.prior=hdr.prior,
             stuff.min=stuff.min, stuff.max=stuff.max, estimand=estimand,
             plot.bounds=plot.bounds, plot.pf=plot.pf,
             plot.sens=plot.sens,  plot.prior=plot.prior,
             color.bounds=color.bounds,
             color1.pf=color1.pf,
             color2.pf=color2.pf,
             color1.sens=color1.sens,
             color2.sens=color2.sens,
             color.prior=color.prior, ymax=ymax, ...)
    
  }
  if (estimand %in% c("RR", "RRT", "RRC")){
    fit.pf <- y.pf <- hdr.pf <- fit.sens <- y.sens <- hdr.sens <- fit.prior <-
      y.prior <- hdr.prior <- NULL

    fit.pf <- locfit(~lp(stuff.pf, nn=.3, deg=3), renorm=TRUE,
                     ev=rbox(cut=.3), xlim=c(0, 1e10), maxk=500)
    
    fit.sens <- locfit(~lp(stuff.sens, nn=.3, deg=3),
                       renorm=TRUE,
                       ev=rbox(cut=.3), xlim=c(0, 1e10), maxk=500)

    x <- seq(from=min(c(fit.pf$box, fit.sens$box)),
             to=max(c(quantile(stuff.pf, .999),
               quantile(stuff.sens, .999))),
             length.out=4001)
    
    
    if (plot.pf){
      y.pf <- predict(fit.pf, newdata=x)
      hdr.pf <- hdr(den=list(x=x, y=y.pf), all.modes=FALSE, prob=percent)
    }
    if (plot.sens){
      y.sens <- predict(fit.sens, newdata=x)
      hdr.sens  <- hdr(den=list(x=x, y=y.sens), all.modes=FALSE, prob=percent)
    }
    if (plot.prior){
      warning("Plotting of prior not allowed for estimand RR, RRT, or RRC\n")
      plot.prior <- FALSE
    }


    if (is.null(ymax)){
      ymax <- max(c(y.pf, y.sens, y.prior))
    }

    makeplot(x=x, y.pf=y.pf, y.sens=y.sens, y.prior=y.prior,
             hdr.pf=hdr.pf, hdr.sens=hdr.sens, hdr.prior=hdr.prior,
             stuff.min=stuff.min, stuff.max=stuff.max, estimand=estimand,
             plot.bounds=plot.bounds, plot.pf=plot.pf,
             plot.sens=plot.sens,  plot.prior=plot.prior,
             color.bounds=color.bounds,
             color1.pf=color1.pf,
             color2.pf=color2.pf,
             color1.sens=color1.sens,
             color2.sens=color2.sens,
             color.prior=color.prior, ymax=ymax, ...)
    
  }
  if (estimand %in% c("logRR", "logRRT", "logRRC")){

    fit.pf <- y.pf <- hdr.pf <- fit.sens <- y.sens <- hdr.sens <- fit.prior <-
      y.prior <- hdr.prior <- NULL

    fit.pf <- locfit(~lp(stuff.pf, nn=.3, deg=3), renorm=TRUE,
                     ev=rbox(cut=.3), maxk=500)
    
    fit.sens <- locfit(~lp(stuff.sens, nn=.3, deg=3),
                       renorm=TRUE,
                       ev=rbox(cut=.3), maxk=500)
    
    x <- seq(from=min(c(fit.pf$box, fit.sens$box)),
             to=max(c(quantile(stuff.pf, .999),
               quantile(stuff.sens, .999))),
             length.out=4001)
    
    
    if (plot.pf){
      y.pf <- predict(fit.pf, newdata=x)
      hdr.pf <- hdr(den=list(x=x, y=y.pf), all.modes=FALSE, prob=percent)
    }
    if (plot.sens){
      y.sens <- predict(fit.sens, newdata=x)
      hdr.sens  <- hdr(den=list(x=x, y=y.sens), all.modes=FALSE, prob=percent)
    }
    if (plot.prior){
      warning("Plotting of prior not allowed for estimand logRR, logRRT, or logRRC\n")
      plot.prior <- FALSE
    }


    if (is.null(ymax)){
      ymax <- max(c(y.pf, y.sens, y.prior))
    }

    makeplot(x=x, y.pf=y.pf, y.sens=y.sens, y.prior=y.prior,
             hdr.pf=hdr.pf, hdr.sens=hdr.sens, hdr.prior=hdr.prior,
             stuff.min=stuff.min, stuff.max=stuff.max, estimand=estimand,
             plot.bounds=plot.bounds, plot.pf=plot.pf,
             plot.sens=plot.sens,  plot.prior=plot.prior,
             color.bounds=color.bounds,
             color1.pf=color1.pf,
             color2.pf=color2.pf,
             color1.sens=color1.sens,
             color2.sens=color2.sens,
             color.prior=color.prior, ymax=ymax, ...)
  }


  
} ## end plotSimpleTable





## function to actually do the plotting
"makeplot" <- function(x, y.pf, y.sens, y.prior, hdr.pf, hdr.sens, hdr.prior,
                       stuff.min, stuff.max, estimand,
                       plot.bounds, plot.pf, plot.sens,  plot.prior,
                       color.bounds,
                       color1.pf,
                       color2.pf,
                       color1.sens,
                       color2.sens,
                       color.prior, ymax, ...){
  
  if (estimand == "ATE"){
    supplement = "Average Treatment Effect"
  }
  if (estimand == "ATT"){
    supplement = "Average Treatment Effect Within the Treated"
  }
  if (estimand == "ATC"){
    supplement = "Average Treatment Effect Within Controls"
  }
  if (estimand == "RR"){
    supplement = "Relative Risk"
  }
  if (estimand == "RRT"){
    supplement = "Relative Risk Within the Treated"
  }
  if (estimand == "RRC"){
    supplement = "Relative Risk Within Controls"
  }
  if (estimand == "logRR"){
    supplement = "Log Relative Risk"
  }
  if (estimand == "logRRT"){
    supplement = "Log Relative Risk Within the Treated"
  }
  if (estimand == "logRRC"){
    supplement = "Log Relative Risk Within Controls"
  }

  if (stuff.min == -Inf){
    stuff.min <- min(x)
  }
  if (stuff.max == Inf){
    stuff.max <- max(x)
  }
  

  if (estimand %in% c("ATE", "ATT", "ATC")){
    plot(0:1, 0:1, type="n", xlim=c(-1, 1),
         ylim=c(0, ymax), ylab="Density",
         xlab=supplement, xaxt="n", ...)
  
    axis(side=1, at=c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1), ...)
  }
  else{
    plot(0:1, 0:1, type="n", xlim=c(min(x, stuff.min), max(x, stuff.max)),
         ylim=c(0, ymax), ylab="Density",
         xlab=supplement, ...)
  }
    
  
  if (plot.sens){
    hdrCheckInterval(hdr.sens, x, y.sens)
    n.hdr.sens <- ncol(hdr.sens$hdr)/2
    polygon(c(x, rev(x)), c(y.sens, rep(0, length(y.sens))),
            border=NA, col=color1.sens)
    
    for (i in 1:n.hdr.sens){
      indic <- x >= hdr.sens$hdr[,i*2-1] & x <= hdr.sens$hdr[,i*2]
      polygon(c(x[indic], rev(x[indic])), c(y.sens[indic], rep(0, sum(indic))),
              col=color2.sens,
              border=NA)
    }
  }
  
  if (plot.prior){
    hdrCheckInterval(hdr.prior, x, y.prior)   
    n.hdr.prior <- ncol(hdr.prior$hdr)/2
    lines(x, y.prior, lwd=3, col=color.prior)
  }
  
  if (plot.pf){
    hdrCheckInterval(hdr.pf, x, y.pf)
    n.hdr.pf <- ncol(hdr.pf$hdr)/2
    lines(x, y.pf, lwd=3, lty=2, col=color1.pf)
    
    for (i in 1:n.hdr.pf){
      indic <- x >= hdr.pf$hdr[,i*2-1] & x <= hdr.pf$hdr[,i*2]
      lines(x[indic], y.pf[indic], col=color2.pf, lwd=3)
    }
  }
  
  if (plot.bounds){
    if (estimand %in% c("RRT", "logRRT")){
      par(lend=2, ljoin=1)
      arrows(stuff.min, -0.02*ymax, stuff.max, -0.02*ymax, col=color.bounds,
             lwd=5)
    }
    else if (estimand == "logRRC"){
      par(lend=2, ljoin=1)
      arrows(stuff.max, -0.02*ymax, stuff.min, -0.02*ymax, col=color.bounds,
             lwd=5)
    }
    else{
      polygon(c(stuff.min, stuff.min, stuff.max, stuff.max),
              c(-.01*ymax, -.0275*ymax, -.0275*ymax, -.01*ymax),
              col=color.bounds, border=NA)
    }
  }
  
} ## end makeplot





