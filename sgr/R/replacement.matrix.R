#rm(list=ls())
#source("~/lavori/Rdevel/SGR1.0/R/pfakegood.R")
#source("~/lavori/Rdevel/SGR1.0/R/pfakebad.R")
#source("~/lavori/Rdevel/SGR1.0/R/dgBetaD.R")
#source("~/lavori/Rdevel/sgr1.0/R/pfake.R")
#Q <- 5
#p <- c(0,.5)
#fake <- "bad"

replacement.matrix <- function(Q=5,p=c(0,0),gam=c(1,1),del=c(1,1),
    fake.model=c("uninformative","average","slight","extreme")) {
  
  if (length(p)<2) p <- rep(p,2)
  if (length(gam)<2) gam <- rep(gam,2)
  if (length(del)<2) del <- rep(del,2)
  
  fake.model <- match.arg(fake.model)
  if (fake.model!="uninformative") {
    MF <- model.fake.par(fake.model)
    gam <- MF$gam; del <- MF$del
  }
  
  if (sum(p)>1) {
    warning("sum(p) must be not greater than 1")
    p <- p/sum(p) 
  }
  
  fake <- "both"
  if ((p[1]==0)&(p[2]>0)) fake <- "bad"
  if ((p[2]==0)&(p[1]>0)) fake <- "good"  
    
  if (fake=="good") {
    FUN <- get("pfakegood")  
    p <- p[1]
    gam <- gam[1]
    del <- del[1]
  }
  if (fake=="bad")  {
    FUN <- get("pfakebad")  
    p <- p[2]
    gam <- gam[2]
    del <- del[2]
  }
  if (fake=="both") FUN <- get("pfake")
  
  R <- matrix(,Q,Q)
  for (i in 1:Q) {
    for (j in 1:Q){
      R[i,j] <- FUN(j,i,p,Q,gam,del)
    }
  }
  
  #return(list(gam,del))
  return(R)
}

## example
#replacement.matrix(Q=7) 
#replacement.matrix(Q=7,p=c(.1,.1))
#replacement.matrix(Q=7,p=c(.5,0),gam=8,del=2.5)

#replacement.matrix("bad",Q=7,p=.5) 
#replacement.matrix(Q=7,p=.5)
#replacement.matrix(Q=7,p=.5,gam=8,del=2.5)

#replacement.matrix(Q=7,p=c(0,.4),fake.model="extreme")
#replacement.matrix(Q=7,p=c(.4,0),fake.model="extreme")
#replacement.matrix(Q=7,p=c(.4,.4),fake.model="slight")
