"vis.binom" <-
function(){

    if( !requireNamespace('tcltk', quietly = TRUE) ) stop('This function depends on the tcltk package')

  if(!exists('slider.env')) slider.env<<-new.env()


  n <- 10  ; assign('n',tcltk::tclVar(n),envir=slider.env)
  p <-  0.5; assign('p',tcltk::tclVar(p),envir=slider.env)

  sn <- 0  ; assign('sn',tcltk::tclVar(sn), envir=slider.env)
  sp <- 0  ; assign('sp',tcltk::tclVar(sp), envir=slider.env)


  binom.refresh <- function(...){
    n <- as.numeric(evalq(tcltk::tclvalue(n), envir=slider.env))
    p <- as.numeric(evalq(tcltk::tclvalue(p), envir=slider.env))

    sn <- as.numeric(evalq(tcltk::tclvalue(sn), envir=slider.env))
    sp <- as.numeric(evalq(tcltk::tclvalue(sp), envir=slider.env))

    mu <- p*n
    sd <- sqrt(n*p*(1-p))

    if(sn){
      xx <- seq(-1,n+1, length=250)
      plot(xx,dnorm(xx,mu,sd), type='l', col='green',
           ylim=range(0,dnorm(mu,mu,sd),dbinom( seq(0,n), n, p)),
           xlab='x', ylab='Probability')
      if(sp){
        points( seq(0,n), dpois( seq(0,n), mu ), type='h', col='blue')
        points( seq(0,n), dpois( seq(0,n), mu ), pch='-', col='blue',cex=2)
      }
      abline(h=0)
      lines(xx, dnorm(xx,mu,sd), col='green')
      points( seq(0,n), dbinom( seq(0,n), n, p), type='h' )
      points( seq(0,n), dbinom( seq(0,n), n, p), type='p' )

    } else {
      if(sp){
        plot( seq(0,n), dpois( seq(0,n), mu ), type='h', col='blue',
             xlim=c(-1,n+1), xlab='x', ylab='Probability',
             ylim=range(0,dpois( seq(0,n), mu), dbinom(seq(0,n),n,p)))
        points( seq(0,n), dpois( seq(0,n), mu ), pch='-', col='blue',cex=2)
        points( seq(0,n), dbinom( seq(0,n), n, p), type='h')
      } else {
        plot( seq(0,n), dbinom( seq(0,n), n, p), type='h',
             xlim=c(-1,n+1), xlab='x', ylab='Probability')
      }
      abline(h=0)
      points( seq(0,n), dbinom( seq(0,n), n, p) )
    }

    title( paste("Mean =",round(mu,3),"Std. Dev. =",round(sd,3)) )

  }

  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m,'Visualizing the Binomial Distribution')
  tcltk::tkwm.geometry(m,'+0+0')

  # n
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='n', width='10'), side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=binom.refresh, from=1, to=100,
                       orient='horiz',
                       resolution=1, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=n), envir=slider.env)

  # p
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='p', width='10'), side='right')
  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=binom.refresh, from=0, to=1,
                       orient='horiz',
                       resolution=0.01, showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=p), envir=slider.env)

  # show normal
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(sc <- tcltk::tkcheckbutton(fr, command=binom.refresh),
         side='left')
  tcltk::tkpack(tcltk::tklabel(fr, text='Show Normal Approximation',width='25'),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sn), envir=slider.env)

  # show poisson
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(sc <- tcltk::tkcheckbutton(fr, command=binom.refresh),
         side='left')
  tcltk::tkpack(tcltk::tklabel(fr, text='Show Poisson Approximation',width='25'),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sp), envir=slider.env)

  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=binom.refresh), side='left')
  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')

}

