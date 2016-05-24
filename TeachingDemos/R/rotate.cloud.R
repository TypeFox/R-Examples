"rotate.cloud" <-
function(x, ...){

  if(!requireNamespace('tcltk', quietly=TRUE)){stop('The tcltk package is needed')}
  if(!exists('slider.env')) slider.env <<-new.env()

  lab1 <- 'z'; assign('lab1', tcltk::tclVar(lab1), envir=slider.env)
  lab2 <- 'y'; assign('lab2', tcltk::tclVar(lab2), envir=slider.env)
  lab3 <- 'x'; assign('lab3', tcltk::tclVar(lab3), envir=slider.env)

  val1 <-  40; assign('val1', tcltk::tclVar(val1), envir=slider.env)
  val2 <-   0; assign('val2', tcltk::tclVar(val2), envir=slider.env)
  val3 <- -60; assign('val3', tcltk::tclVar(val3), envir=slider.env)

  cloud.options <- list(...)

  cloud.refresh <- function(...){

    lab1 <- evalq(tcltk::tclvalue(lab1), envir=slider.env)
    lab2 <- evalq(tcltk::tclvalue(lab2), envir=slider.env)
    lab3 <- evalq(tcltk::tclvalue(lab3), envir=slider.env)

    val1 <- as.numeric(evalq(tcltk::tclvalue(val1), envir=slider.env))
    val2 <- as.numeric(evalq(tcltk::tclvalue(val2), envir=slider.env))
    val3 <- as.numeric(evalq(tcltk::tclvalue(val3), envir=slider.env))


    sl <- list(val1,val2,val3)
    names(sl) <- c(lab1,lab2,lab3)

    cloud.options$x <- x
    cloud.options$screen <- sl

    print( do.call('cloud',cloud.options) )


  }

  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m,'Rotate Cloud plot')
  tcltk::tkwm.geometry(m,'+0+0')

  # one
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(e <- tcltk::tkentry(fr, width=2), side='left')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=cloud.refresh, from=-180, to=180,
                       orient='horiz',
                       resolution=1, showvalue=T),
         side='left')

  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=val1), envir=slider.env)
  assign('e',e,envir=slider.env)
  evalq(tcltk::tkconfigure(e,textvariable=lab1), envir=slider.env)

  # two
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(e <- tcltk::tkentry(fr, width=2), side='left')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=cloud.refresh, from=-180, to=180,
                       orient='horiz',
                       resolution=1, showvalue=T),
         side='left')

  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=val2), envir=slider.env)
  assign('e',e,envir=slider.env)
  evalq(tcltk::tkconfigure(e,textvariable=lab2), envir=slider.env)

  # three
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(e <- tcltk::tkentry(fr, width=2), side='left')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=cloud.refresh, from=-180, to=180,
                       orient='horiz',
                       resolution=1, showvalue=T),
         side='left')

  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=val3), envir=slider.env)
  assign('e',e,envir=slider.env)
  evalq(tcltk::tkconfigure(e,textvariable=lab3), envir=slider.env)

  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=cloud.refresh),side='left')

  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')


}

