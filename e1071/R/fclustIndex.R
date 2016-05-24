fclustIndex <- function ( y, x, index= "all" )
{

  clres <- y
###########################################################################
################SESSION 1: MEASURES#########################################
###########################################################################

  gath.geva <- function (clres,x)#for m=2
    {
      xrows <- dim(clres$me)[1]
      xcols <- dim(clres$ce)[2]
      ncenters <- dim(clres$centers)[1]
      scatter <- array(0.0, c(xcols, xcols, ncenters))
      scatternew <- array(0.0, c(xcols, xcols, ncenters))
      fhv <-as.double(0)
      apd <-as.double(0)
      pd <-  as.double(0)
      control <- as.double(0)

      for (i in 1:ncenters){
        paronomastis <- as.double(0)
        paronomastis2 <- as.double(0)

        for (j in 1:xrows){
          paronomastis <- paronomastis+clres$me[j,i]

          diff <- x[j,]-clres$ce[i,]
          scatternew[,,i] <- clres$me[j,i]*(t(t(diff))%*%t(diff))
          scatter[,,i] <- scatter[,,i]+scatternew[,,i]
        }#xrows
        scatter[,,i] <- scatter[,,i]/paronomastis


        for (j in 1:xrows){
          diff <- x[j,]-clres$ce[i,]
          control <- (t(diff)%*%solve(scatter[,,i]))%*%t(t(diff))
          if (control<1.0)
            paronomastis2 <- paronomastis2+clres$me[j,i]
          ##   else
          ##     cat("...")
        }#xrows
        fhv <- fhv+sqrt(det(scatter[,,i]))

        apd <- apd+paronomastis2/sqrt(det(scatter[,,i]))

        pd <- pd+paronomastis2

      }#ncenters
      pd <- pd/fhv
      apd <- apd/ncenters

      retval <- list(fuzzy.hypervolume=fhv,average.partition.density=apd, partition.density=pd)
      return(retval)
    }

  xie.beni <- function(clres){#for all m
    xrows <- dim(clres$me)[1]
    minimum<--1
    error <- clres$within#sd
    ncenters <- dim(clres$centers)[1]
    for (i in 1:(ncenters-1)){
      for (j in (i+1):ncenters){
        diff<- clres$ce[i,]-clres$ce[j,]
        diffdist <- t(diff)%*%t(t(diff))
        if (minimum==-1)
          minimum <- diffdist
        if (diffdist<minimum)
          minimum <- diffdist
      }}
    xiebeni <- error/(xrows*minimum)
    return(xiebeni)
  }

  fukuyama.sugeno <- function(clres){#for all m
    xrows <- dim(clres$me)[1]
    ncenters <-dim(clres$centers)[1]
    error <- clres$within#sd
    k2<-as.double(0)

    meancenters <- apply(clres$ce,2,mean)
    for (i in 1:ncenters){
      paronomastis3 <- as.double(0)
      for (j in 1:xrows){
        paronomastis3 <- paronomastis3+(clres$me[j,i]^2)}

      diff <- clres$ce[i,]-meancenters
      diffdist <- t(diff)%*%t(t(diff))
      k2 <- k2 +paronomastis3*diffdist
    }#ncenters

    fukuyamasugeno<-error-k2
    return(fukuyamasugeno)
  }

  partition.coefficient <- function(clres){
    xrows <- dim(clres$me)[1]

    partitioncoefficient <- sum(apply(clres$me^2,1,sum))/xrows
    return(partitioncoefficient)
  }

  partition.entropy <- function(clres){
    xrows <- dim(clres$me)[1]
    ncenters <- dim(clres$centers)[1]
    partitionentropy <- 0.0
    for (i in 1:xrows){
      for (k in 1:ncenters){
        if (clres$me[i,k]!=0.0)
          partitionentropy<- partitionentropy+(clres$me[i,k]*log(clres$me[i,k]))
      }}
    partitionentropy<-partitionentropy/((-1)*xrows)
    ##partitionentropy <- sum(apply(clres$me*log(clres$me),1,sum))/((-1)*xrows)
    return(partitionentropy)
  }

  separation.index <- function(clres, x)
    {
      xrows <- dim(clres$me)[1]
      xcols <- dim(x)[2]
      ncenters <- dim(clres$centers)[1]
      maxcluster <- double(ncenters)
      minimum <- -1.0
      ##hardpartition <- matrix(0,xrows,ncenters)

      ##for (i in 1:xrows)
      ## hardpartition[i,clres$cl[i]] <- 1
      for (i in 1:ncenters){
        maxcluster[i] <- max(dist(matrix(x[clres$cl==i],ncol=xcols)))
      }
      maxdia <- maxcluster[rev(order(maxcluster))[1]]
      for (i in 1:(ncenters-1)){
        for (j in (i+1):(ncenters)){
          for (m in 1:xrows){
            if (clres$cl[m]==i){
              for (l in 1:xrows){
                if (clres$cl[l]==j){
                  diff <- x[m,]-x[l,]
                  diffdist <- sqrt(t(diff)%*%t(t(diff)))
                  fraction <- diffdist/maxdia
                  if (minimum==-1)
                    minimum <- fraction
                  if (fraction<minimum)
                    minimum <- fraction
                }}}}}}
      return(minimum)
    }


  proportion.exponent <- function(clres)
    {
      k <- dim(clres$centers)[2]
      xrows <- dim(clres$me)[1]

      bexp <- as.integer(1)
      for (j in 1:xrows){
        greatint <- as.integer(1/max(clres$me[j,]))
        aexp <- as.integer(0)
        for (l in 1:greatint){
          aexp <- aexp + (-1)^(l+1)*(gamma(k+1)/(gamma(l+1)*gamma(k-l+1)))*(1-l*max(clres$me[j,]))^(k-1)
          ##if (aexp==0.0){
          ##print("aexp=0.0")
          ##print(j)
          ##}
        }
        bexp <- bexp * aexp

      }
      proportionexponent <- -log(bexp)
      return(proportionexponent)
    }

###########################################################################
################SESSION 2: MAIN PROGRAM####################################
###########################################################################

  index <- pmatch(index, c("gath.geva", "xie.beni", "fukuyama.sugeno",
                           "partition.coefficient", "partition.entropy",
                           "proportion.exponent", "separation.index", "all"))


  if (is.na(index))
    stop("invalid clustering index")
  if (index == -1)
    stop("ambiguous index")

   vecallindex <- numeric(9)

  if (any(index==1) || (index==8)){
    gd <- gath.geva(clres, x)
    vecallindex[1] <- gd$fuzzy
    vecallindex[2] <- gd$average
    vecallindex[3] <- gd$partition}
  if (any(index==2) || (index==8))
    vecallindex[4] <- xie.beni(clres)
  if (any(index==3) || (index==8))
    vecallindex[5] <- fukuyama.sugeno(clres)
  if (any(index==4) || (index==8))
    vecallindex[6] <- partition.coefficient(clres)
  if (any(index==5) || (index==8))
    vecallindex[7] <- partition.entropy(clres)
  if (any(index==6) || (index==8))
    vecallindex[8] <- proportion.exponent(clres)
  if (any(index==7) || (index==8))
      vecallindex[9] <- separation.index(clres,x)

  names(vecallindex) <- c("fhv", "apd", "pd", "xb", "fs", "pc", "pe",
                          "pre", "si")


  if (index < 8) {
      if (index == 1)
          vecallindex <- vecallindex[1:3]
      else
          vecallindex <- vecallindex[index + 2]
  }

  return(vecallindex)
}


