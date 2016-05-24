`Waldci` <-
function(cmat, estp, varp, varcor, alternative="two.sided", conf.level=0.95, dist="MVN")

  {

k=ncol(cmat)
m=nrow(cmat)

if(any( c(length(estp), length(varp), length(varcor) )!=k ))
 { stop("lengths of estp, nadj, varp, varcor, and ncol(cmat) must be the same")}

    estC <- cmat %*% estp

    CorrMat<-corrMatgen(CM=cmat, varp=varcor)

varC <- (cmat^2) %*% (varp) 
  
dist<-match.arg(dist, choices=c("MVN", "N"))

switch(dist,

"MVN"={
 
    # Berechnung des MVNORM quantils

if(alternative=="two.sided")
 { quanti <- qmvnorm(p=conf.level, sigma=CorrMat, tail="both.tails")$quantile
   stderr <- sqrt(varC)
   lCI <- estC-quanti*stderr
   uCI <- estC+quanti*stderr
 }

else
 {
  if(alternative=="less")
   { quanti <- qmvnorm(p=conf.level, sigma=CorrMat, tail="lower.tail")$quantile
     stderr <- sqrt(varC)
     lCI <- rep(-Inf, m)
     uCI <- estC+quanti*stderr
   }

  else
   {
    if(alternative=="greater")
     {quanti <- qmvnorm(p=conf.level, sigma=CorrMat, tail="upper.tail")$quantile
      stderr <- sqrt(varC)
      lCI <- estC+quanti*stderr
      uCI <- rep(Inf, m)
     }
   }
 }
},

"N"={

if(alternative=="two.sided")
 { quanti <- qnorm(p=1-(1-conf.level)/2)
   stderr <- sqrt(varC)
   lCI <- estC-quanti*stderr
   uCI <- estC+quanti*stderr
 }

else
 {
  if(alternative=="less")
   { quanti <- qnorm(p=conf.level)
     stderr <- sqrt(varC)
     lCI <- rep(-Inf, m)
     uCI <- estC+quanti*stderr
   }

  else
   {
    if(alternative=="greater")
     {quanti <- qnorm(p=1-conf.level)
      stderr <- sqrt(varC)
      lCI <- estC+quanti*stderr
      uCI <- rep(Inf, m)
     }
   }
 }

})

conf.int<-cbind(lCI,uCI)
colnames(conf.int)<-c("lower","upper")

quantile<-quanti
attr(x=quantile, which="dist")<-dist

out<-list(conf.int=conf.int,
 alternative=alternative,
 conf.level=conf.level,
 quantile=quantile,
 corrmat=CorrMat,
 dist=dist)

class(out)<-"sci"

return(out)

  }

