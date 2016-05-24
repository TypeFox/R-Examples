"vis.t" <-
function(){

    if( !requireNamespace('tcltk', quietly = TRUE) ) stop('This function depends on the tcltk package')

    if(!exists('slider.env')) slider.env<<-new.env()

  df <- 1; assign('df',tcltk::tclVar(df),envir=slider.env)
  sn <- 0; assign('sn',tcltk::tclVar(sn),envir=slider.env)

  xmin <- -5; assign('xmin',tcltk::tclVar(xmin),envir=slider.env)
  xmax <- 5; assign('xmax',tcltk::tclVar(xmax),envir=slider.env)
  ymin <- 0; assign('ymin',tcltk::tclVar(ymin),envir=slider.env)
  ymax <- round(dnorm(0,0,1),2); assign('ymax',tcltk::tclVar(ymax),envir=slider.env)



  t.refresh <- function(...){

    df <- as.numeric(evalq(tcltk::tclvalue(df), envir=slider.env))
    sn <- as.numeric(evalq(tcltk::tclvalue(sn), envir=slider.env))

    xmin <- as.numeric(evalq(tcltk::tclvalue(xmin), envir=slider.env))
    xmax <- as.numeric(evalq(tcltk::tclvalue(xmax), envir=slider.env))
    ymin <- as.numeric(evalq(tcltk::tclvalue(ymin), envir=slider.env))
    ymax <- as.numeric(evalq(tcltk::tclvalue(ymax), envir=slider.env))

    xx <- seq(xmin,xmax, length=500)
    yyt <- dt(xx,df)

    if(sn){
      yyn <- dnorm(xx)
      plot(xx,yyn, lwd=3, col='skyblue', type='l',
           xlim=c(xmin,xmax), ylim=c(ymin,ymax),
           xlab='x', ylab='')
      lines(xx,yyt,lwd=2)
    } else {
      plot(xx,yyt,type='l', xlim=c(xmin,xmax), ylim=c(ymin,ymax),
           ylab='',xlab='x',lwd=2)
    }
  }


  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m,'Visualizing the t-Distribution')
  tcltk::tkwm.geometry(m,'+0+0')

  # df
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='d.f.', width='5'),side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=t.refresh, from=1, to=50,
                       orient='horiz',
                       resolution=1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=df),envir=slider.env)

  # show normal
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='Show Normal Distribution', width='25'),side='right')
  tcltk::tkpack(sc <- tcltk::tkcheckbutton(fr, command=t.refresh),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sn),envir=slider.env)


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


  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=t.refresh),side='left')

  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')

}

