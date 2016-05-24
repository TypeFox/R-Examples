"vis.gamma" <-
function(){

  if(!exists('slider.env')) slider.env<<-new.env()
  if(!requireNamespace('tcltk', quietly = TRUE)) {
      stop('This function needs the tcltk package')
  }

  shape <- 1; assign('shape',tcltk::tclVar(shape),envir=slider.env)
  rate  <- 1; assign('rate',tcltk::tclVar(rate),envir=slider.env)
  scale <- 1; assign('scale',tcltk::tclVar(scale),envir=slider.env)
  mean  <- 1; assign('mean', tcltk::tclVar(mean), envir=slider.env)
  sd    <- 1; assign('sd',tcltk::tclVar(sd), envir=slider.env)

  se <- 0; assign('se', tcltk::tclVar(se), envir=slider.env)
  sc2 <- 0; assign('sc2', tcltk::tclVar(sc2), envir=slider.env)
  sg <- 1; assign('sg', tcltk::tclVar(sg), envir=slider.env)

  xmin <- 0; assign('xmin',tcltk::tclVar(xmin),envir=slider.env)
  xmax <- 5; assign('xmax',tcltk::tclVar(xmax),envir=slider.env)
  ymin <- 0; assign('ymin',tcltk::tclVar(ymin),envir=slider.env)
  ymax <- 1; assign('ymax',tcltk::tclVar(ymax),envir=slider.env)

  old.shape <- shape
  old.rate  <- rate
  old.scale <- scale
  old.mean  <- mean
  old.sd    <- sd


  gamma.refresh <- function(...){

    shape <- as.numeric(evalq(tcltk::tclvalue(shape), envir=slider.env))
    rate  <- as.numeric(evalq(tcltk::tclvalue(rate), envir=slider.env))
    scale <- as.numeric(evalq(tcltk::tclvalue(scale), envir=slider.env))
    mean  <- as.numeric(evalq(tcltk::tclvalue(mean), envir=slider.env))
    sd    <- as.numeric(evalq(tcltk::tclvalue(sd), envir=slider.env))

    if ( shape != old.shape ) {

      mean <- shape * scale
      sd <- round( sqrt(shape)*scale, 6 );

      try(eval(parse(text=paste("tcltk::tclvalue(mean)<-",
                       mean,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(sd)<-",
                       sd,sep="")),envir=slider.env));
      old.shape <<- shape; old.mean <<- mean; old.sd <<- sd

    }

    if ( rate != old.rate ) {

      scale <- round(1/rate, 6)
      mean <- shape * scale
      sd <- round( sqrt(shape)*scale, 6 );

      try(eval(parse(text=paste("tcltk::tclvalue(scale)<-",
                       scale,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(mean)<-",
                       mean,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(sd)<-",
                       sd,sep="")),envir=slider.env));
      old.rate <<- rate; old.scale <<- scale;
      old.mean <<- mean; old.sd <<- sd

    }

    if ( scale != old.scale ) {

      rate <- round(1/scale, 6)
      mean <- shape * scale
      sd <- round( sqrt(shape)*scale, 6 );

      try(eval(parse(text=paste("tcltk::tclvalue(rate)<-",
                       rate,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(mean)<-",
                       mean,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(sd)<-",
                       sd,sep="")),envir=slider.env));
      old.rate <<- rate; old.scale <<- scale;
      old.mean <<- mean; old.sd <<- sd

    }

    if ( mean != old.mean ) {

      shape <- round( (mean/sd)^2, 6 )
      scale <- round( mean/shape, 6 )
      rate <- round(1/scale, 6)

      try(eval(parse(text=paste("tcltk::tclvalue(rate)<-",
                       rate,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(shape)<-",
                       shape,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(scale)<-",
                       scale,sep="")),envir=slider.env));
      old.shape <<- shape; old.rate <<- rate; old.scale <<- scale;
      old.mean <<- mean; old.sd <<- sd

    }


    if ( sd != old.sd ) {

      shape <- round( (mean/sd)^2, 6 )
      scale <- round( mean/shape, 6 )
      rate <- round(1/scale, 6)

      try(eval(parse(text=paste("tcltk::tclvalue(rate)<-",
                       rate,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(shape)<-",
                       shape,sep="")),envir=slider.env));
      try(eval(parse(text=paste("tcltk::tclvalue(scale)<-",
                       scale,sep="")),envir=slider.env));
      old.shape <<- shape; old.rate <<- rate; old.scale <<- scale;
      old.mean <<- mean; old.sd <<- sd

    }


    se <- as.numeric(evalq(tcltk::tclvalue(se), envir=slider.env))
    sc2 <- as.numeric(evalq(tcltk::tclvalue(sc2), envir=slider.env))
    sg <- as.numeric(evalq(tcltk::tclvalue(sg), envir=slider.env))

    xmin <- as.numeric(evalq(tcltk::tclvalue(xmin), envir=slider.env))
    xmax <- as.numeric(evalq(tcltk::tclvalue(xmax), envir=slider.env))
    ymin <- as.numeric(evalq(tcltk::tclvalue(ymin), envir=slider.env))
    ymax <- as.numeric(evalq(tcltk::tclvalue(ymax), envir=slider.env))

    xx <- seq(xmin,xmax, length=500)

    plot(xx,xx, xlim=c(xmin,xmax),ylim=c(ymin,ymax),
         xlab='x', ylab='y',type='n')

    if(se) {
      yye <- dexp(xx,1/mean)
      lines(xx,yye, lwd=3, col='green')
      lines(c(mean,mean),c(ymin,dexp(mean,1/mean)), lty=2, col='green')
      lines(c(mean,mean*2), dexp(mean*2, 1/mean)*c(1,1), lty=2, col='green')
    }

    if(sc2) {
      yyc <- dchisq(xx,mean)
      lines(xx,yyc, lwd=3, col='blue')
      lines(c(mean,mean),c(ymin,dchisq(mean,mean)), lty=2, col='blue')
      lines(c(mean,mean+sqrt(2*mean)), dchisq(mean+sqrt(2*mean), mean)*c(1,1),
            lty=2, col='blue')
    }

    if(sg) {
      yyg <- dgamma(xx,shape,rate)
      lines(xx,yyg, lwd=2)
      lines(c(mean,mean),c(ymin,dgamma(mean,shape,rate)), lty=2)
      lines(c(mean,mean+sd), dgamma(mean+sd, shape, rate)*c(1,1),
            lty=2)
    }

  }


  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m,'Visualizing the Gamma Distribution')
  tcltk::tkwm.geometry(m,'+0+0')

  # shape
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Shape', width='10'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=gamma.refresh, from=0.1, to=10,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=shape),envir=slider.env)

  # rate
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Rate', width='10'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=gamma.refresh, from=0.1, to=10,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=rate),envir=slider.env)

  # scale
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Scale', width='10'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=gamma.refresh, from=0.1, to=10,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=scale),envir=slider.env)

  # mean
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Mean', width='10'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=gamma.refresh, from=0.1, to=100,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=mean),envir=slider.env)

  # sd
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='S.D.', width='10'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=gamma.refresh, from=0.1, to=40,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sd),envir=slider.env)


  # show exponential
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(sc <- tcltk::tkcheckbutton(fr, command=gamma.refresh),
         side='left')
  tcltk::tkpack(tcltk::tklabel(fr, text='Show Exponential Distribution', width='25'),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=se),envir=slider.env)

  # show chisquared
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(sc <- tcltk::tkcheckbutton(fr, command=gamma.refresh),
         side='left')
  tcltk::tkpack(tcltk::tklabel(fr, text='Show Chi-squared Distribution', width='25'),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sc2),envir=slider.env)

  # show gamma
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(sc <- tcltk::tkcheckbutton(fr, command=gamma.refresh),
         side='left')
  tcltk::tkpack(tcltk::tklabel(fr, text='Show Gamma Distribution', width='25'),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sg),envir=slider.env)


  # xmin
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Xmin:', width=6), side='left')
  tcltk::tkpack(e <- tcltk::tkentry(fr,width=8), side='left')
  assign('e',e,envir=slider.env)
  evalq(tcltk::tkconfigure(e, textvariable=xmin), envir=slider.env)

  # xmax
  tcltk::tkpack(tcltk::tklabel(fr, text='Xmax:', width=6), side='left')
  tcltk::tkpack(e <- tcltk::tkentry(fr,width=8), side='left')
  assign('e',e,envir=slider.env)
  evalq(tcltk::tkconfigure(e, textvariable=xmax), envir=slider.env)

  # ymin
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Ymin:', width=6), side='left')
  tcltk::tkpack(e <- tcltk::tkentry(fr,width=8), side='left')
  assign('e',e,envir=slider.env)
  evalq(tcltk::tkconfigure(e, textvariable=ymin), envir=slider.env)

  # ymax
  tcltk::tkpack(tcltk::tklabel(fr, text='Ymax:', width=6), side='left')
  tcltk::tkpack(e <- tcltk::tkentry(fr,width=8), side='left')
  assign('e',e,envir=slider.env)
  evalq(tcltk::tkconfigure(e, textvariable=ymax), envir=slider.env)


  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=gamma.refresh),side='left')

  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')

}

