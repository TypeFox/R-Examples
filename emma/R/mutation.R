mutation  <- function(mom,x,pr.mut,b,C,t)
{
  dx <- ncol(x$xspace)

  limits <- c()
  for(i in 1:dx){
    limits <- rbind(limits,range(x$xspace[,i]))
  }
  probab <- pr.mut[t]

  for(i in 1:nrow(mom)){
    x.mut <- sample(c(0,1),1,prob=c(1-probab,probab))
    if(x.mut==1){
      var <- sample(1:dx,1)
      flip <- sample(c(0,1),1)
      if(flip==0){
        UB <- limits[var,2]
        r <- runif(1)
        arg <- (UB-mom[i,var])*(1-r^((1-t/C)^b))
        mom[i,var] <- mom[i,var]+arg
      }
      if(flip==1){
        LB <- limits[var,1]
        r <- runif(1)
        arg <- (mom[i,var]-LB)*(1-r^((1-t/C)^b))
        mom[i,var] <- mom[i,var]-arg
      }
    }
  }
  return(mom)
}

