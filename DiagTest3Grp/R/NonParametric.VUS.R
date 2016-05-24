NonParametric.VUS <- function(x,y,z,alpha=0.05,NBOOT=50,FisherZ=FALSE)
  {
    ######################################################################################################################################################
    ########This function provides Nonparametric estimate on VUS =Pr(X1<=X2<=X3)
    ########
    ####1. Input:
    ###x,y,z:vectors are the test marker measurements for the 3 ordinal group (D-,D0,D+)
    ####2. Output: nonparametric VUS estimate 
    ######################################################################################################################################################
    x <- na.exclude(x)
    y <- na.exclude(y)
    z <- na.exclude(z)

    ###sample size
    n.minus <- length(x)
    n0 <- length(y)
    n.plus <- length(z)

    ##normal sample means
    mu.minus <- mean(x,na.rm=TRUE)
    mu0 <- mean(y,na.rm=TRUE)
    mu.plus <- mean(z,na.rm=TRUE)

    ##normal sample SDs
    s.minus <- sd(x,na.rm=TRUE)
    s0 <- sd(y,na.rm=TRUE)
    s.plus <- sd(z,na.rm=TRUE)
    
    ###data summary
    dat.summary <- data.frame(n=c(n.minus,n0,n.plus),mu=c(mu.minus,mu0,mu.plus),sd=c(s.minus,s0,s.plus),row.names=c("D-","D0","D+"))

    ###non-parametric VUS estimate by ranking the triplets, each from a diagnosis group
    temp.merge1 <- merge(data.frame(ID=rep(1,n.minus),x=x),data.frame(ID=rep(1,n0),y=y),by.x=1,by.y=1,all.x=T,all.y=T)
    temp.merge2 <- merge(temp.merge1,data.frame(ID=rep(1,n.plus),z=z),by.x=1,by.y=1,all.x=T,all.y=T)
    flag <- 1*(temp.merge2$x<=temp.merge2$y & temp.merge2$y<=temp.merge2$z)###if no ties,weight 1
    flag[temp.merge2$x==temp.merge2$y | temp.merge2$y==temp.merge2$z] <- 0.5##ties have weight 0.5
    VUS <- mean(flag)
    if(FisherZ) VUS <- FisherZ(VUS)
    ###bootstrap variance and CI has to be obtained by bootstrapping outside the function
       
    return(list(dat=list(x=x,y=y,z=z),dat.summary=dat.summary,estimate=VUS))
  }

