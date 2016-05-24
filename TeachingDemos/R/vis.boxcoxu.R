vis.boxcoxu.old <-
function(lambda = sample( c(-1, -0.5, 0, 1/3, 1/2, 1, 2), 1)) {

  if( !requireNamespace('tkrplot', quietly = TRUE) ) stop('This function depends on the tkrplot package being available')


  y <- rnorm(1000, 7, 2)
  if( min(y) <= 0 ) y <- y - min(y)+0.05
  if (lambda==0) {
    y <- exp(y)
  } else {
    y <- y^(1/lambda)
  }

    if(!exists('slider.env')) slider.env <<-new.env()
  #library(tcltk)

  lam <- 1 ; assign('lam',tcltk::tclVar(lam), envir=slider.env)

  bc.refresh <- function(...){
    lam <- as.numeric(evalq(tcltk::tclvalue(lam), envir=slider.env))

    old.par <- par(mfcol=c(2,2))
    on.exit(par(old.par))

    ty <- bct(y,lam)
    hist(y, prob=T, xlab='x', main='Histogram of x')
    xx <- seq(min(y),max(y), length=250)
    lines(xx, dnorm( xx, mean(y), sqrt(var(y)) ))
    qqnorm(y, xlab='x')
    qqline(y)

    hist(ty, prob=T, xlab='Transformed x', main = 'Histogram of Transformed x')
    xx <- seq(min(ty),max(ty), length=250)
    lines(xx,dnorm(xx, mean(ty), sqrt(var(ty)) ) )

    qqnorm(ty, xlab='Transformed x')
    qqline(ty)
  }

  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m, 'Box Cox Transform')
  tcltk::tkwm.geometry(m,'+0+0')

  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='lambda', width='10'), side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=bc.refresh, from=-2, to=3, orient='horiz',
                       resolution=0.1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=lam), envir=slider.env)

  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=bc.refresh), side='left')
  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')

}




vis.boxcoxu <- function(lambda = sample( c(-1,-0.5,0,1/3,1/2,1,2), 1),
                        y, xlab=deparse(substitute(y)),
                        hscale=1.5, vscale=1.5, wait=FALSE) {

  if( missing(y) ) {
    if(missing(xlab)) xlab <- 'y'
    y <- rnorm(1000, 7, 2)
    if( min(y) <= 0 ) y <- y - min(y) + 0.05
    if(lambda==0) {
      y <- exp(y)
    } else {
      y <- y^(1/lambda)
    }
  }


  lam <- tcltk::tclVar()
  tcltk::tclvalue(lam) <- 1
  hsc <- tcltk::tclVar()
  tcltk::tclvalue(hsc) <- hscale
  vsc <- tcltk::tclVar()
  tcltk::tclvalue(vsc) <- hscale

  replot <- function(...) {
    tmp.l <- as.numeric(tcltk::tclvalue(lam))

    par(mfcol=c(2,2))

    ty <- bct(y,tmp.l)
    hist(y, prob=TRUE, xlab=xlab, main = paste('Histogram of',xlab))
    xx <- seq(min(y),max(y), length=250)
    lines(xx, dnorm(xx, mean(y), sd(y)) )
    qqnorm(y, xlab=xlab)
    qqline(y)

    hist(ty, prob=TRUE, xlab=paste("Transformed",xlab),
         main=paste("Histogram of Transformed",xlab))
    xx <- seq(min(ty),max(ty), length=250)
    lines(xx,dnorm(xx, mean(ty), sd(ty)))
    qqnorm(ty, xlab=paste("Transformed",xlab))
    qqline(ty)
  }

  tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(tt, "Box Cox Demo")

  img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)
  tcltk::tkpack(img, side='top')

  tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='lambda: '), side='left', anchor='s')
  tcltk::tkpack(tcltk::tkscale(fr, variable=lam, orient='horizontal',
                 command=function(...) tkrplot::tkrreplot(img,
                   hscale=as.numeric(tcltk::tclvalue(hsc)),
                   vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                 from=-2, to=4, resolution=.05), side='right')

  tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
  tcltk::tkpack(tcltk::tkbutton(tfr, text="Refresh", command=function() tkrplot::tkrreplot(img,
                                         hscale=as.numeric(tcltk::tclvalue(hsc)),
                                         vscale=as.numeric(tcltk::tclvalue(vsc)) ) ),
                  side='left',anchor='s')

  tcltk::tkpack(tcltk::tkbutton(tfr, text="Exit", command=function()tcltk::tkdestroy(tt)),
             side='right',anchor='s')

  tcltk::tkpack(tfr <- tcltk::tkframe(tt), side='bottom', fill='x')
  tcltk::tkpack(tcltk::tklabel(tfr,text="Hscale: "), side='left')
  tcltk::tkpack(tcltk::tkentry(tfr,textvariable=hsc,width=6), side='left')
  tcltk::tkpack(tcltk::tklabel(tfr,text="      Vscale: "), side='left')
  tcltk::tkpack(tcltk::tkentry(tfr,textvariable=vsc,width=6), side='left')

  if(wait) {
    tcltk::tkwait.window(tt)
    return( list(lambda = as.numeric(tcltk::tclvalue(lam)),
                 y = y,
                 ty = bct(y,as.numeric(tcltk::tclvalue(lam))) ) )
  } else {
    return(invisible(NULL))
  }
}

