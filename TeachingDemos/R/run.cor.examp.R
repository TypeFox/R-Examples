"run.old.cor.examp" <-
function(n=100,seed) {
  if (!missing(seed)){ set.seed(seed) }
  if(!requireNamespace('tcltk', quietly = TRUE)){stop('The tcltk package is needed')}

  x <- scale(matrix(rnorm(2*n,0,1), ncol=2))
  x <- x %*% solve( chol( cor(x) ) )
  xr <- range(x)

  cor.refresh <- function(...) {
    r <- slider(no=1)
    if ( r == 1 ) {
      cmat <- matrix( c(1,0,1,0),2 )
    } else if (r == -1) {
      cmat <- matrix( c(1,0,-1,0),2 )
    } else {
      cmat <- chol( matrix( c(1,r,r,1),2) )
    }
    new.x <- x %*% cmat

    plot(new.x, xlab='x',ylab='y', xlim=xr, ylim=xr)
    title(paste("r = ",round(cor(new.x[,1],new.x[,2]),3)))
  }

  slider( cor.refresh, 'Correlation', -1, 1, 0.01, 0,
         title="Correlation Demo")
  cor.refresh()
}

run.cor.examp <- function(n=100,seed,vscale=1.5,hscale=1.5,wait=FALSE) {

    if( !requireNamespace('tkrplot', quietly = TRUE) ) stop('This function depends on the tkrplot package being available')

    if(!missing(seed) ) set.seed(seed)

    x <- scale(matrix(rnorm(2*n,0,1), ncol=2))
    x <- x %*% solve( chol( cor(x) ) )
    xr <- range(x,-x)

    hsc <- tcltk::tclVar()
    tcltk::tclvalue(hsc) <- hscale
    vsc <- tcltk::tclVar()
    tcltk::tclvalue(vsc) <- vscale

    r <- tcltk::tclVar()
    tcltk::tclvalue(r) <- 0

    replot <- function(...) {
        tmp.r <- as.numeric(tcltk::tclvalue(r))
        if( tmp.r == 1 ) {
            cmat <- matrix( c(1,0,1,0),2 )
        } else if(  tmp.r == -1 ) {
            cmat <- matrix( c(1,0,-1,0),2 )
        } else {
            cmat <- chol( matrix( c(1,tmp.r,tmp.r,1),2) )
        }
        new.x <- x %*% cmat

        plot(new.x, xlab='x', ylab='y', xlim=xr, ylim=xr)
        title(paste("r =", round( tmp.r, 3)))
    }

    tt <- tcltk::tktoplevel()
    tcltk::tkwm.title(tt, "Cor2 Example")

    img <- tkrplot::tkrplot(tt, replot, vscale=vscale, hscale=hscale)
    tcltk::tkpack(img, side='top')


    tcltk::tkpack(fr <- tcltk::tkframe(tt), side='top')
    tcltk::tkpack(tcltk::tklabel(fr,text='r: '), side='left',anchor='s')
    tcltk::tkpack(tcltk::tkscale(fr, variable=r, orient='horizontal',
                   command=function(...) tkrplot::tkrreplot(img,
                    hscale=as.numeric(tcltk::tclvalue(hsc)),
                    vscale=as.numeric(tcltk::tclvalue(vsc)) ),
                   from=-1, to=1,
                   resolution=0.01), side='right')


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

    if(wait){
        tcltk::tkwait.window(tt)
        tmp.r <- as.numeric(tcltk::tclvalue(r))
        if( tmp.r == 1 ) {
            cmat <- matrix( c(1,0,1,0),2 )
        } else if(  tmp.r == -1 ) {
            cmat <- matrix( c(1,0,-1,0),2 )
        } else {
            cmat <- chol( matrix( c(1,tmp.r,tmp.r,1),2) )
        }
        new.x <- x %*% cmat
        return( list(x=new.x[,1], y=new.x[,2]) )
    } else {
        return(invisible(NULL))
    }
}
