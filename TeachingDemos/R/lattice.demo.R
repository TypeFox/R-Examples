"lattice.demo" <-
function(x,y,z, show3d=TRUE){

  if(!requireNamespace('tcltk', quietly = TRUE)){stop('The tcltk package is needed')}
  if(!requireNamespace('lattice', quietly=TRUE)){stop('the lattice package is needed')}
  if(!exists('slider.env')) slider.env <<- new.env()


  center <- mean(z); assign('center',tcltk::tclVar(center), envir=slider.env)
  width <- diff(range(z))/20*3; assign('width',tcltk::tclVar(width), envir=slider.env)

  s3d <- 1; assign('s3d', tcltk::tclVar(s3d), envir=slider.env)

  lattice.refresh <- function(...){
    center <- as.numeric(evalq(tcltk::tclvalue(center), envir=slider.env))
    width <- as.numeric(evalq(tcltk::tclvalue(width), envir=slider.env))

    s3d <- as.numeric(evalq(tcltk::tclvalue(s3d), envir=slider.env))

    shingle.min <- max(min(z), center-width/2)
    shingle.max <- min(max(z), center+width/2)

    shingle.scaled.range <- c( (shingle.min-min(z))/diff(range(z)),
                               (shingle.max-min(z))/diff(range(z))) - 0.5

    if(s3d){
      print(lattice::xyplot(y~x|shingle(z,rbind(range(z),c(shingle.min,shingle.max))),
                   index.cond=list(2),
                   strip=lattice::strip.custom(strip.names=TRUE,strip.levels=TRUE),
                   par.strip.text=list(cex=0.75)),
            split=c(1,1,1,2), more=T)

      print(lattice::cloud(y~z+x, panel=function(x,y,z,...){
        lattice::panel.cloud(x,y,z,panel.3d.cloud=function(x,y,z,groups,...){
          lattice::panel.3dscatter(x,y,z,
                          groups= factor(x>shingle.scaled.range[1] & x <shingle.scaled.range[2]),
                          ...)
          lattice::panel.3dwire(x=shingle.scaled.range,
                       y=c(-.5,.5), z=rep(-.5,4), at=c(-.57,.57), ...)
        },...) }), split=c(1,2,1,2),more=F)

    } else {
      print(lattice::xyplot(y~x|shingle(z,rbind(range(z),c(shingle.min,shingle.max))),
                   index.cond=list(2),
                   strip=lattice::strip.custom(strip.names=TRUE,strip.levels=TRUE),
                   par.strip.text=list(cex=0.75)
                   ),
            split=c(1,1,1,1), more=F)

    }

  }

  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m,'Trellis/Lattice Demo')
  tcltk::tkwm.geometry(m,'+0+0')

  # center
  tcltk::tkpack(fr <- tcltk::tkframe(m),side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='center', width='10'), side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=lattice.refresh, from=min(z), to=max(z),
                       orient='horiz',
                       resolution=diff(range(z))/25, showvalue=T),
         side='left')

  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=center),envir=slider.env)

  # width
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='width', width='10'), side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=lattice.refresh, from=diff(range(z))/20,
                       to=diff(range(z)),orient='horiz',
                       resolution=diff(range(z))/20, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=width), envir=slider.env)

  # show 3d
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(sc <- tcltk::tkcheckbutton(fr, command=lattice.refresh),
         side='left')
  tcltk::tkpack(tcltk::tklabel(fr, text='Show 3-D plot', width='25'),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=s3d), envir=slider.env)

  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=lattice.refresh),side='left')

  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')

}

