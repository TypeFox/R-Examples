"vis.normal" <-
function(){

    if( !requireNamespace('tcltk', quietly = TRUE) ) stop('This function depends on the tcltk package')


  if(!exists('slider.env')) slider.env<<-new.env()
  #library(tcltk)

  mu <- 0; assign('mu',tcltk::tclVar(mu),envir=slider.env)
  sd <- 1; assign('sd',tcltk::tclVar(sd),envir=slider.env)
  s2 <- 1; assign('s2',tcltk::tclVar(s2),envir=slider.env)
  xmin <- -5; assign('xmin',tcltk::tclVar(xmin),envir=slider.env)
  xmax <- 5; assign('xmax',tcltk::tclVar(xmax),envir=slider.env)
  ymin <- 0; assign('ymin',tcltk::tclVar(ymin),envir=slider.env)
  ymax <- round(dnorm(0,0,.5),2); assign('ymax',tcltk::tclVar(ymax),envir=slider.env)

  sd.old <- sd
  s2.old <- s2

  norm.refresh <- function(...){

    mu <- as.numeric(evalq(tcltk::tclvalue(mu), envir=slider.env))
    sd <- as.numeric(evalq(tcltk::tclvalue(sd), envir=slider.env))
    s2 <- as.numeric(evalq(tcltk::tclvalue(s2), envir=slider.env))

    if(sd != sd.old) {
      s2 <- round(sd^2,5); # assign('s2',tclVar(s2),envir=slider.env)
      try(eval(parse(text=paste("tcltk::tclvalue(s2)<-",
                       s2,sep="")),envir=slider.env));
      sd.old <<- sd; s2.old <<- s2
    }

    if(s2 != s2.old) {
      s2 <- as.numeric(evalq(tcltk::tclvalue(s2), envir=slider.env))
      sd <- round(sqrt(s2),5); # assign('sd',tclVar('sd'), envir=slider.env)
      try(eval(parse(text=paste("tcltk::tclvalue(sd)<-",
                       sd,sep="")),envir=slider.env));
      sd.old <<- sd; s2.old <<- s2
    }

    xmin <- as.numeric(evalq(tcltk::tclvalue(xmin), envir=slider.env))
    xmax <- as.numeric(evalq(tcltk::tclvalue(xmax), envir=slider.env))
    ymin <- as.numeric(evalq(tcltk::tclvalue(ymin), envir=slider.env))
    ymax <- as.numeric(evalq(tcltk::tclvalue(ymax), envir=slider.env))

    xx <- seq(xmin,xmax, length=500)
    yy <- dnorm(xx,mu,sd)
    plot(xx,yy,type='l', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
         ylab='',xlab='x')
    lines(c(mu,mu),c(par('usr')[3],dnorm(0,0,sd)), lty=2, col='blue')
    lines(c(mu,mu+sd), dnorm(sd,0,sd)*c(1,1), lty=2, col='blue')

  }



  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m,'Visualizing the Normal Distribution')
  tcltk::tkwm.geometry(m,'+0+0')

  # mean
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Mean', width='20'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=norm.refresh, from=-3, to=3,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=mu),envir=slider.env)

  # sd
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Standard Deviation', width='20'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=norm.refresh, from=.5, to=3,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sd),envir=slider.env)

  # variance
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Variance', width='20'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=norm.refresh, from=.25, to=9,
                       orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=s2),envir=slider.env)


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


  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=norm.refresh),side='left')

  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')

}

