run.power.examp.old <-
function(){
  if(!requireNamespace('tcltk', quietly = TRUE)){stop('The tcltk package is needed')}
  slider( power.refresh,
         c('Sample Size','Standard Deviation','True Difference',
           'Alpha level'),
         c(1,0.25,-1,0.01),
         c(100,5,3,0.99),
         c(1,0.25,0.1,0.01),
         c(1,1,1,0.05),
         title="Power Demo") }

run.power.examp <- function(hscale=1.5, vscale=1.5, wait=FALSE) {

    if( !requireNamespace('tkrplot', quietly = TRUE) ) stop('This function depends on the tkrplot package being available')

    n <- tcltk::tclVar()
    stdev <- tcltk::tclVar()
    diff <- tcltk::tclVar()
    alpha <- tcltk::tclVar()
    xmin <- tcltk::tclVar()
    xmax <- tcltk::tclVar()

    tcltk::tclvalue(n) <- 1
    tcltk::tclvalue(stdev) <- 1
    tcltk::tclvalue(diff) <- 1
    tcltk::tclvalue(alpha) <- 0.05
    tcltk::tclvalue(xmin) <- -2
    tcltk::tclvalue(xmax) <- 4

    hsc <- tcltk::tclVar()
    tcltk::tclvalue(hsc) <- hscale
    vsc <- tcltk::tclVar()
    tcltk::tclvalue(vsc) <- vscale

    out <- numeric(1)

    replot <- function(...) {
        out <<- power.examp( as.numeric(tcltk::tclvalue(n)),
                            as.numeric(tcltk::tclvalue(stdev)),
                            as.numeric(tcltk::tclvalue(diff)),
                            as.numeric(tcltk::tclvalue(alpha)),
                            as.numeric(tcltk::tclvalue(xmin)),
                            as.numeric(tcltk::tclvalue(xmax)) )
    }

    tt <- tcltk::tktoplevel()
    tcltk::tkwm.title(tt, "Power Demo")

    img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)
    tcltk::tkpack(img, side='left')

    tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top', fill='x')
    tcltk::tkpack(tcltk::tklabel(fr, text="n: "), side='left')
    tcltk::tkpack(tdspinner(fr, values=c(1,2,3,4,5,10,20,30,40,50),
                      width=5, textvariable=n,
                      command=function(...) tkrplot::tkrreplot(img,
                        hscale=as.numeric(tcltk::tclvalue(hsc)),
                        vscale=as.numeric(tcltk::tclvalue(vsc)) )
                      ), side='left')

    tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top',fill='x')
    tcltk::tkpack(tcltk::tklabel(fr, text="Standard Deviation: "), side='left')
    tcltk::tkpack(tcltk::tkscale(fr, variable=stdev, orient='horizontal',
                   command=function(...) tkrplot::tkrreplot(img,
                     hscale=as.numeric(tcltk::tclvalue(hsc)),
                     vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                   from=0.1, to=4, resolution=.05), side='right')

    tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top',fill='x')
    tcltk::tkpack(tcltk::tklabel(fr, text="True Difference: "), side='left')
    tcltk::tkpack(tcltk::tkscale(fr, variable=diff, orient='horizontal',
                   command=function(...) tkrplot::tkrreplot(img,
                     hscale=as.numeric(tcltk::tclvalue(hsc)),
                     vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                   from=0, to=4, resolution=.05), side='right')

    tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top',fill='x')
    tcltk::tkpack(tcltk::tklabel(fr, text="alpha: "), side='left')
    tcltk::tkpack(tcltk::tkscale(fr, variable=alpha, orient='horizontal',
                   command=function(...) tkrplot::tkrreplot(img,
                     hscale=as.numeric(tcltk::tclvalue(hsc)),
                     vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                   from=0.001, to=0.2, resolution=0.001), side='right')



    tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='top', fill='x')
    tcltk::tkpack(tcltk::tklabel(tfr,text="x min: "), side='left')
    tcltk::tkpack(tcltk::tkentry(tfr,textvariable=xmin,width=6), side='left')
    tcltk::tkpack(tcltk::tklabel(tfr,text="      x max: "), side='left')
    tcltk::tkpack(tcltk::tkentry(tfr,textvariable=xmax,width=6), side='left')



    tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
    tcltk::tkpack(tcltk::tklabel(tfr,text="Hscale: "), side='left')
    tcltk::tkpack(tcltk::tkentry(tfr,textvariable=hsc,width=6), side='left')
    tcltk::tkpack(tcltk::tklabel(tfr,text="      Vscale: "), side='left')
    tcltk::tkpack(tcltk::tkentry(tfr,textvariable=vsc,width=6), side='left')

    tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
    tcltk::tkpack(tcltk::tkbutton(tfr, text="Refresh", command=function() tkrplot::tkrreplot(img,
                                         hscale=as.numeric(tcltk::tclvalue(hsc)),
                                         vscale=as.numeric(tcltk::tclvalue(vsc)) ) ),
           side='left',anchor='s')

    tcltk::tkpack(tcltk::tkbutton(tfr, text="Exit", command=function()tcltk::tkdestroy(tt)),
           side='right',anchor='s')

    if(wait) {
        tcltk::tkwait.window(tt)
        return(list( n=as.numeric(tcltk::tclvalue(n)),
                     stdev=as.numeric(tcltk::tclvalue(stdev)),
                     diff=as.numeric(tcltk::tclvalue(diff)),
                     alpha=as.numeric(tcltk::tclvalue(alpha)),
                     power=out
                    ))
    } else {
        return(invisible(NULL))
    }
}
