"mle.demo" <-
function(x=rnorm(10, 10, 2), start.mean = mean(x)-start.sd,
                     start.sd = 1.2* sqrt(var(x)) ){

    if( !requireNamespace('tcltk', quietly = TRUE) ) stop('This function depends on the tcltk package.')

  if(!exists('slider.env')) slider.env <<- new.env()

  mu <- start.mean; assign('mu',tcltk::tclVar(mu),envir=slider.env)
  sig <- start.sd;  assign('sig',tcltk::tclVar(sig),envir=slider.env)

  .mu <- .sig <- .ll <- numeric(0)

  mle.refresh <- function(...){
    mu <- as.numeric(evalq(tcltk::tclvalue(mu), envir=slider.env))
    sig <- as.numeric(evalq(tcltk::tclvalue(sig), envir=slider.env))

    old.par <- par(no.readonly=T)
    on.exit(par(old.par))

    par(mar=c(5,4,0,2)+.1)

    .mu <<- c(.mu, mu)
    .sig <<- c(.sig, sig)

    ll <- sum( dnorm(x, mu, sig, TRUE) )
    .ll <<- c(.ll,ll)

    layout( matrix( 1:3, ncol=1 ), heights=c(2,1,1))

    xx <- seq( min(x) - 1.2 * (max(x)-min(x)),
               max(x) + 1.2 * (max(x)-min(x)), length=250)
    plot(xx, dnorm(xx, mu, sig), type='l',
         ylim=c(0,dnorm(0,0,0.5*sqrt(var(x)))),xlab='x', ylab='Likelihood')
    segments(x, 0, x, dnorm(x, mu, sig))
    points(x,dnorm(x, mu, sig))
    points(x,rep(0,length(x)))

    text(xx[1], dnorm(0,0,0.5*sqrt(var(x)))*.9, paste("Log Likelihood =",
                                                      format(ll, digit=5)),
         adj=0,cex=3)

    plot(.mu, .ll, xlab=expression(mu), ylab='Log Likelihood')
    points(mu,ll, pch=16, col='red')
    plot(.sig, .ll, xlab=expression(sigma), ylab='Log Likelihood')
    points(sig, ll, pch=16, col='red')

  }

  m <- tcltk::tktoplevel()
  tcltk::tkwm.title(m,'Maximum Likelihood Estimation')
  tcltk::tkwm.geometry(m, '+0+0')

  # mu
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='mu', width='10'), side='right')

  tmp <- pretty( c( start.mean - 2*start.sd, start.mean + 3*start.sd),
                100)

  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=mle.refresh, from=min(tmp),
                       to=max(tmp),
                       orient='horiz',
                       resolution=tmp[2] - tmp[1],showvalue=T),
         side='left')
  assign('sc',sc, envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=mu), envir=slider.env)

  # sigma
  tcltk::tkpack(fr <- tcltk::tkframe(m), side='top')
  tcltk::tkpack(tcltk::tklabel(fr, text='sigma', width='10'), side='right')

  tmp <- pretty( c( 0.5*start.sd, 2*start.sd), 100)

  tcltk::tkpack(sc <- tcltk::tkscale(fr, command=mle.refresh, from=min(tmp),
                       to=max(tmp),
                       orient='horiz',
                       resolution=tmp[2]-tmp[1], showvalue=T),
         side='left')
  assign('sc',sc,envir=slider.env)
  evalq(tcltk::tkconfigure(sc, variable=sig), envir=slider.env)



  tcltk::tkpack(tcltk::tkbutton(m, text="Refresh", command=mle.refresh), side='left')
  tcltk::tkpack(tcltk::tkbutton(m, text="Exit", command=function()tcltk::tkdestroy(m)),
         side='right')

  return(invisible(x))
}

