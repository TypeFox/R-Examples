
clustIndex <- function( y, x, index="all" )
  {
    clres <- y
################################################################### 
#########SESSION 1: Necessary non-commented functions##############
###################################################################    
    
    withinss <- function(clobj, x){
      
      retval <- rep(0, nrow(clobj$centers))
      x <- (x - clobj$centers[clobj$cluster, ])^2
      for(k in 1:nrow(clobj$centers)){
        retval[k] <- sum(x[clobj$cluster==k,])
      }
      retval
    }

    varwithinss <- function(x, centers, cluster)
      {
        nrow<-dim(centers)[1]
        nvar<-dim(x)[2]
        varwithins <- matrix(0, nrow,nvar)
        x <- (x - centers[cluster, ])^2
        for(l in 1:nvar){
          for(k in 1:nrow){
            varwithins[k,l] <- sum(x[cluster==k,l])
          }
        }
        return(varwithins)
      }

    maxmindist <-  function (clsize, distscen)
      {
        ##only for binary data
        
        ncl <- length(clsize)
        
        npairs <- 0
        for (i in 1:ncl)
          npairs <- npairs + clsize[i]*(clsize[i]-1)/2
        
        ##minimum distance
        mindw <- 0
        nfound <- distscen[1]
        i <- 1
        while (nfound < npairs)
          {
            if ((nfound+distscen[i+1]) < npairs)
              {
                mindw <- mindw + i*distscen[i+1]
                nfound <- nfound+distscen[i+1]
              }
            else
              {
                mindw <- mindw + i*(npairs-nfound)
                nfound <- nfound+distscen[i+1]
              }
            i <- i+1
          }
        
        ##maximum distance
        maxdw <- 0
        nfound <- 0
        i <- length(distscen) - 1
        while (nfound < npairs)
          {
            if ((nfound+distscen[i+1]) < npairs)
              {
                maxdw <- maxdw + i*distscen[i+1]
                nfound <- nfound+distscen[i+1]
              }
            else
              {
                maxdw <- maxdw + i*(npairs-nfound)
                nfound <- nfound+distscen[i+1]
              }
            i <- i-1
          }
        
        minmaxd <-  list(mindw=mindw, maxdw=maxdw)
        
        return(minmaxd)
      }

    gss <- function(x, clsize, withins)
      {
        n <- sum(clsize)
        k <- length(clsize)
        allmean <- apply(x,2,mean)
        dmean <- sweep(x,2,allmean,"-")
        allmeandist <- sum(dmean^2)
        wgss <- sum(withins)
        bgss <- allmeandist - wgss
        
        zgss <- list(wgss=wgss, bgss=bgss)
        return(zgss)
      }

    vargss <- function(x, clsize, varwithins)
      {
        nvar<-dim(x)[2]
        n <- sum(clsize)
        k <- length(clsize)
        varallmean<-rep(0,nvar)
        varallmeandist<-rep(0,nvar)
        varwgss<-rep(0,nvar)
        for (l in 1:nvar)
          varallmean[l] <- mean(x[,l])
        vardmean <- sweep(x,2,varallmean,"-")
        for (l in 1:nvar)
          {
            varallmeandist[l] <- sum((vardmean[,l])^2)
            varwgss[l] <- sum(varwithins[,l])
          }
        varbgss <- varallmeandist - varwgss
        vartss<-varbgss+varwgss
        zvargss <- list(vartss=vartss, varbgss=varbgss)
        return(zvargss)
      }

    count <- function(x)
      {
        nr <- nrow(x)
        nc <- ncol(x)
        d <- integer(nc+1)
        
        retval <- .C("count", xrows=nr, xcols=nc, x=as.integer(x), d=d,
                     PACKAGE="cclust")
        
        d <- retval$d
        
        names(d) <- 0:nc
        return(d)
      }

    ttww <- function(x, clsize, cluster)
      {
        n <- sum(clsize)
        k <- length(clsize)
        w<-0
        tt <- cov(x)*n
        for (l in 1:k)
          w<- w+cov(x[cluster==l,])*clsize[l]
        zttw <- list(tt=tt, w=w)
        return(zttw)
      }
################# END SESSION 1 ##########################################
###########################################################################
################SESSION 2: INDEXES#########################################
    calinski <- function(zgss, clsize)
      {
        n <- sum(clsize)
        k <- length(clsize)
        vrc <- (zgss$bgss/(k-1))/(zgss$wgss/(n-k))
        return(vrc=vrc)
      }
    
    cindex <- function (withins, minmaxd, clsize)
      {
        dw <- sum(withins*clsize)    
        ##c-index
        cindex <- (dw -minmaxd$mindw)/(minmaxd$maxdw - minmaxd$mindw)
        return(cindex)
      }
    
    
    
    db <- function(withins, centers, cluster)
      {
        mse <- withins/table(cluster)
        r <- outer(mse, mse, "+") / as.matrix(dist(centers, diag=TRUE))
        diag(r) <- 0
        db <- mean(apply(r,1,max))
        return(db)
      }
    
    
    hartigan <- function(zgss)
      {
        hart <- log(zgss$bgss/zgss$wgss)
        return(hart)
      }
    
    
    ratkowsky <- function(zvargss, clsize)
      {
        k <- length(clsize)
        ratio<-mean(sqrt(zvargss$varbgss/zvargss$vartss))
        rat <- ratio/sqrt(k)
        return(rat)
      }
    
    
    scott <- function(zttw, clsize)
      {
        n <- sum(clsize)
        dettt<-prod(eigen(zttw$tt)$values)
        detw<-prod(eigen(zttw$w)$values)
        scott <- n * log(dettt/detw)
        return(scott)
      }
    
    marriot <- function(zttw, clsize)
      {
        k <- length(clsize)
        detw<-prod(eigen(zttw$w)$values)
        mar <- (k**2) * detw
        return(mar)
      }
    
    
    ball <- function(withins, clsize)
      {
        ball <- sum(withins)/length(clsize)
      }
    
    
    tracecovw <- function(zttw)
      {
        trcovw <- sum(diag(cov(zttw$w)))
        return(trcovw)
      }
    
    
    tracew <- function(zttw)
      {
        tracew <- sum(diag(zttw$w))
        return(tracew)
      }
    
    
    friedman <- function(zttw)
      {
        b <- zttw$tt-zttw$w
        fried <- sum(diag(solve(zttw$w)%*%b))
        return(fried)
      }
    
    
    rubin <- function(zttw)
      {
        dettt<-prod(eigen(zttw$tt)$values)
        detw<-prod(eigen(zttw$w)$values)
        friedm <- dettt/detw
        return(friedm)
      }
    
    
    ssi <- function (centers, clsize)
      {
        ncl <- dim(centers)[1]
        nvar <- dim(centers)[2]
        n <- sum(clsize)
        
        cmax <- apply(centers, 2, max)
        cmin <- apply(centers, 2, min)
        cord <- apply(centers, 2, order)
        cmaxi <- cord[ncl,]
        cmini <- cord[1,]
        
        meanmean <- mean(centers)
        absmdif <- abs(apply(centers, 2, mean) - meanmean)
        span <- cmax - cmin
        csizemax <- clsize[cmaxi]
        csizemin <- clsize[cmini]

        hiest <- nvar
        hiestw <- hiest * max(span) * max(max(csizemax), max(csizemin)) * exp(-min(absmdif))
        
        sist <- sum(span)/hiest
        
        sistw <- (span * exp(-absmdif)) %*% sqrt(csizemax*csizemin) / hiestw
        
        return(list(ssi=sist, ssiw=sistw))
      }
    
    likelihood <- function (x, centers, cluster)
      {
        n <- nrow(x)
        l <- 0
        
        for (i in 1:n)
          l <- l - log(prod(x[i,]*centers[cluster[i],] +
                            (1-x[i,])*(1-centers[cluster[i],])))
        
        return(l)
      }
    
    xu <- function(x, clsize,  zgss)
      {
        n <- sum(clsize)
        k <- length(clsize)
        d <- dim(x)[2]
        
        xuindex <- d * log(sqrt(zgss$wgss/(d*(n^2)))) + log(k) 
        return(xuindex)
      }
##################END SESSION 2###################################
##################################################################
##################SESSION 3: MAIN PROGRAM#########################    


###needed measures 
    ##withins <- withinss1(x, centers, cluster)
    ## varwithins <- varwithinss(x, clres$centers, clres$cluster)
    zgss <- gss(x,clres$size, clres$withins)
    ## zvargss <- vargss(x, clres$size,varwithins)
    zttw <- ttww(x, clres$size, clres$cluster)
    ## minmaxd <-maxmindist(clres$size,distdata)
    ## distdata <- countdist(x)

    
    
###indexes calculations
    index <- pmatch(index, c("calinski", "cindex", "db", "hartigan",
                          "ratkowsky", "scott", "marriot", "ball",
                          "trcovw", "tracew", "friedman",
                          "rubin","ssi","likelihood","xuindex","all"))

    if (is.na(index)) 
      stop("invalid clustering index")
    if (index == -1) 
      stop("ambiguous index")

    vecallindex <- numeric(15)

    if (any(index==1) || (index==16))
      vecallindex[1] <- calinski(zgss,clres$size)
    if (any(index==2) || (index==16))
      {
        if (length(unique(x))==2)
          {
            distdata <- count(x)
            minmaxd <- maxmindist(clres$size,distdata)
            vecallindex[2] <- cindex(clres$withins, minmaxd, clres$size)
          }
        else  vecallindex[2] <- NA
      }
    if (any(index==3) || (index==16))
      vecallindex[3] <- db(clres$withins, clres$centers, clres$cluster)
    if (any(index==4) || (index==16))
      vecallindex[4] <- hartigan(zgss)
    if (any(index==5) || (index==16))
      {
        varwithins <- varwithinss(x, clres$centers, clres$cluster)
        zvargss <- vargss(x, clres$size,varwithins)
        vecallindex[5] <-ratkowsky(zvargss, clres$size)
      }
    if (any(index==6) || (index==16))
       vecallindex[6] <- scott(zttw, clres$size)
    if (any(index==7) || (index==16))
       vecallindex[7] <- marriot(zttw,clres$size)
    if (any(index==8) || (index==16))
       vecallindex[8] <- ball(clres$withins, clres$size)
    if (any(index==9) || (index==16))
       vecallindex[9] <- tracecovw(zttw)
    if (any(index==10) || (index==16))
       vecallindex[10] <- tracew(zttw)
    if (any(index==11) || (index==16))
       vecallindex[11] <- friedman(zttw)
    if (any(index==12) || (index==16))
       vecallindex[12] <- rubin(zttw)
    if (any(index==13) || (index==16))
       vecallindex[13] <- ssi(clres$centers, clres$size)$ssiw 
    if (any(index==14) || (index==16))
       vecallindex[14] <- likelihood(x,clres$centers, clres$cluster)
    if (any(index==15) || (index==16))
       vecallindex[15] <- xu(x, clres$size, zgss)
    
    names(vecallindex) <- c("calinski", "cindex", "db", "hartigan",
                              "ratkowsky", "scott", "marriot", "ball",
                              "trcovw", "tracew", "friedman",
                              "rubin","ssi","likelihood","xuindex")

    if (index < 16)
      vecallindex <- vecallindex[index]

    return(vecallindex)
  }









