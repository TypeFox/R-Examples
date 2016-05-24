"run.ci.examp" <-
function(reps=100,seed, method='z',n=25) {

  if(!requireNamespace('tcltk', quietly = TRUE)){stop('The tcltk package is needed')}
    if (!missing(seed)){ set.seed(seed) }
  data <- matrix( rnorm( n*reps, 100, 10), ncol=n)
  rmeans <- rowMeans(data)

  ci.refresh <- function(...) {
    conf.level=slider(no=1)
    switch(method, Z=,z={
      lower <- qnorm( (1-conf.level)/2, rmeans, 10/sqrt(n))
      upper <- qnorm( 1-(1-conf.level)/2, rmeans, 10/sqrt(n))
    },
           T=,t= {
             cv.l <- qt((1-conf.level)/2, n-1)
             cv.u <- qt(1-(1-conf.level)/2, n-1)
             rsds <- sqrt(apply(data,1,var))/sqrt(n)

             lower <- rmeans+cv.l*rsds
             upper <- rmeans+cv.u*rsds
           },
           BOTH=, Both=, both={
             lz <- qnorm( (1-conf.level)/2, rmeans, 10/sqrt(n))
             uz <- qnorm( 1-(1-conf.level)/2, rmeans, 10/sqrt(n))

             cv.l <- qt((1-conf.level)/2, n-1)
             cv.u <- qt(1-(1-conf.level)/2, n-1)
             rsds <- sqrt(apply(data,1,var))/sqrt(n)

             lt <- rmeans+cv.l*rsds
             ut <- rmeans+cv.u*rsds

             lower <- c(rbind(lt,lz,100))
             upper <- c(rbind(ut,uz,100))

             reps <- reps*3
             rmeans <- rep(rmeans, each=3)
             rmeans[c(F,F,T)] <- NA

           },
           stop("method must be z, t, or both") )


    xr <- 100 + 4.5*c(-1,1)*10/sqrt(n)

    plot(lower,seq(1,reps), type="n", xlim=xr, xlab="Confidence Interval",
         ylab="Index")

    abline( v= qnorm(c((1-conf.level)/2,1-(1-conf.level)/2), 100,
              10/sqrt(n)), col='lightgreen')

    if( method=="both" || method=="Both" || method=="BOTH"){
      title( main="Confidence intervals based on both distributions",
            sub="Upper interval is Z in each pair")
    } else {
      title( main=paste("Confidence intervals based on",method,"distribution"))
    }

    colr <- ifelse( lower > 100, 'blue', ifelse( upper < 100, 'red', 'black') )

    abline(v=100)

    segments(lower,1:reps,upper,1:reps, col=colr)
    points( rmeans, seq(along=rmeans), pch='|', col=colr )
    invisible(NULL)
  }

  slider( ci.refresh, 'Confidence Level', 0.5, 0.995, 0.005, 0.95,
         title="Confidence Interval Demo")

}

