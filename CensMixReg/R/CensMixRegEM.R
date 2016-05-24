########################################################################################################
#######   Begin the EM algorithm for mixture the normal and t-student with censored      ###############
#######                                    update in 2015-09-10                          ###############
########################################################################################################

CensMixRegEM <- function(cc, y, x1, Abetas = NULL, medj= NULL, sigma2 = NULL, nu=NULL, pii = NULL, g = NULL, get.init = TRUE, criteria = TRUE, group = FALSE, family = "T", error = 0.00001, iter.max = 100, obs.prob= FALSE, kmeans.param = NULL)
{

  if(ncol(as.matrix(y)) > 1) stop("This function is only for univariate response y!")
  if((family != "T") && (family != "Normal")) stop(paste("Family",family,"not recognized.\n",sep=" "))
  if((length(g) == 0) && ((length(medj)==0) || (length(sigma2)==0) ||  (length(pii)==0)))  stop("The model is not specified correctly.\n")
  if(get.init == FALSE)
  {
    g <- length(medj)
    if((length(medj) != length(sigma2)) || (length(medj) != length(pii))) stop("The size of the initial values are not compatibles.\n")
    if((family == "Normal" || family == "T")) stop("Family not recognized.\n")
    if(sum(pii) != 1) stop("probability of pii does not sum to 1.\n")
  }

  if((length(g)!= 0) && (g < 1)) stop("g must be greater than 0.\n")

  if (get.init == TRUE)
  {
    if(length(g) == 0) stop("g is not specified correctly.\n")

    k.iter.max <- 50
    n.start    <- 1
    algorithm  <- "Hartigan-Wong"

    if(length(kmeans.param) > 0)
    {
      if(length(kmeans.param$iter.max)  > 0) {k.iter.max <- kmeans.param$iter.max}
      if(length(kmeans.param$n.start)   > 0) {n.start    <- kmeans.param$n.start}
      if(length(kmeans.param$algorithm) > 0) {algorithm  <- kmeans.param$algorithm}
    }

    if(g > 1)
    {
      Abetas <- solve(t(x1)%*%x1)%*%t(x1)%*%y
      yr     <- y-x1%*%Abetas
      init   <- kmeans(yr,g,k.iter.max,n.start,algorithm)
      pii    <- init$size/length(yr)
      medj   <- as.vector(init$centers)
      sigma2 <- init$withinss/init$size
      nu     <- 4
    }

    else{
      Abetas <- solve(t(x1)%*%x1)%*%t(x1)%*%y
      yr     <- y-x1%*%Abetas
      medj   <- mean(yr)
      sigma2 <- var(yr)
      pii    <- 1
      nu     <- 4
    }
  }

  ################################################################################
  ###                       Mixture Student-t                                  ###
  ################################################################################
  if (family == "T")
  {
    start.time  <- Sys.time()
    ERRO        <- 1e-6
    TOLERANCIA  <- 1e-6
    MAX_NU      <- 20
    MIN_NU      <- 2.01
    n           <- length(y)
    p           <- ncol(x1)-1
    x           <- as.matrix(x1[,2:(p+1)])

    beta0       <- Abetas[1] ## intercepto
    betas       <- as.matrix(Abetas[2:(p+1)])   ### parameters of regression with dimension "p"
    varphi      <- rep(0,g)  ### alphas
    mu          <- mu1 <- matrix(0,n,g)

    for (k in 1:g)
    {
      varphi[k] <- beta0+medj[k]
      mu1[,k]   <- x%*%betas
      mu[,k]    <- mu1[,k]+varphi[k]
    }

    teta        <- c(Abetas, medj, sigma2, pii,nu)

    criterio    <- 1
    count       <- 0

    lk          <- sum(log(d.mixedT(cc, y, pii, mu, sigma2, nu)))

    while((criterio > error) && (count <= iter.max))
    { #Start the While
      count   <- count + 1
      tal     <- matrix(0, n, g)
      soma1   <- matrix(0, p,1)
      soma2   <- matrix(0, p, p)

      #u1<-y
      #u2<-y^2

      for (j in 1:g)
      {
        ### E-step: Computing the moments
        aux1t  <- CalMom(mu[,j],sigma2[j],nu,y,family)
        d      <- ((y-mu[,j])^2)/sigma2[j]
        u0     <- (nu+1)/(nu+d)
        u1     <- y*(nu+1)/(nu+d)
        u2     <- y^2*(nu+1)/(nu+d)

        u1[cc==1] <- aux1t$Ey1[cc==1]
        u2[cc==1] <- aux1t$Ey2[cc==1]
        u0[cc==1] <- aux1t$Ey0[cc==1]

        d1      <- dT(cc, y, mu[,j], sigma2[j],nu)

        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin

        d2      <- d.mixedT(cc, y, pii, mu, sigma2,nu)

        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j] <- d1*pii[j]/d2

        ### M-step: update the parameters

        pii[j]    <- (1/n)*sum(tal[,j])
        sigma2[j] <- sum(tal[,j]*(u2+u0*varphi[j]^2+u0*mu1[,j]^2-2*u1*varphi[j]-2*u1*mu1[,j]+2*u0*varphi[j]*mu1[,j])) / sum(tal[,j])
        varphi[j] <- sum(tal[,j]*(u1-u0*mu1[,j]))/ sum(u0*tal[,j])

        soma1     <- soma1 + t(x)%*%diag(tal[,j]/sigma2[j])%*%(u1-u0*varphi[j])
        soma2     <- soma2 + t(x)%*%diag(tal[,j]*c(u0)/sigma2[j])%*%x
      }

      betas      <- solve(soma2)%*%soma1

      pii[g]     <- 1 - (sum(pii) - pii[g])

      zero.pos   <- NULL
      zero.pos   <- which(pii == 0)

      if(length(zero.pos) != 0)
      {
        pii[zero.pos] <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }

      if (pii[1]< 0.5)
      {
        # mu     <- cbind(mu[,2],mu[,1])
        # medj   <- as.vector(c(medj[2],medj[1]))
        varphi   <- as.vector(c(varphi[2], varphi[1]))
        pii      <- as.vector(c(pii[2], pii[1]))
        sigma2   <- as.vector(c(sigma2[2], sigma2[1]))
      }

      beta0      <- sum(pii*varphi)

      for (k in 1:g)
      {
        mu1[,k]  <- x%*%betas
        mu[,k]   <- mu1[,k]+varphi[k]
        medj[k]  <- varphi[k]-beta0
      }

      Abetas     <- c(beta0,betas)
      nuold      <- nu

      ft <- function(nu) sum(log(d.mixedT(cc, y, pii, mu, sigma2,nu) ))
      nu <- optimize(f=ft, interval=c(MIN_NU,MAX_NU),maximum=TRUE,tol=TOLERANCIA)$maximum

      teta   <- c(Abetas, medj, sigma2, pii,nu)
      auxlog <- d.mixedT(cc, y, pii, mu, sigma2,nu)

      if(length(which(auxlog == 0)) > 0) auxlog[which(auxlog == 0)] <- .Machine$double.xmin
      lk1      <- sum(log(auxlog))
      criterio <- abs(lk1/lk-1)
      lk       <- lk1
    } #End the While

    if(criteria == TRUE)
    {
      cl  <- apply(tal, 1, which.max)
      icl <- 0
      for (j in 1:g) { icl <- icl+sum(log(pii[j]*dT(cc, y, mu[,j], sigma2[j],nu))) }
      #icl <- -2*icl+(4*g-1)*log(n)
    }

    #### grafico ajuste
    ##hist(y, breaks = 40,probability=T,col="grey")
    ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
    ##lines(xx,d.mixedSN(xx, pii, mu, sigma2, shape),col="red")
    ######################


    teta            <- c(Abetas, medj, sigma2, pii,nu)
    end.time        <- Sys.time()
    time.taken      <- end.time - start.time

  }

  ################################################################################
  ###                           Mixture of Normal
  ################################################################################
  if (family == "Normal")
  {
    start.time  <- Sys.time()
    n           <- length(y)
    p           <- ncol(x1)-1
    x           <- as.matrix(x1[,2:(p+1)])

    beta0       <- Abetas[1]      ## slope
    betas       <- as.matrix(Abetas[2:(p+1)])   ### parameters of regression the dimension "p"
    varphi      <- rep(0,g) ### alphas
    mu          <- mu1 <- matrix(0,n,g)

    for (k in 1:g)
    {
      varphi[k] <- beta0+medj[k]
      mu1[,k]   <- x%*%betas
      mu[,k]    <- mu1[,k] +varphi[k]
    }

    teta        <- c(Abetas, medj, sigma2, pii)

    criterio    <- 1
    count       <- 0

    lk          <- sum(log(d.mixedN(cc, y, pii, mu, sigma2)))

    while((criterio > error) && (count <= iter.max))
    {
      count <- count + 1
      tal   <- matrix(0, n, g)
      soma1 <- matrix(0, p,1)
      soma2 <- matrix(0, p, p)

      u1    <- y
      u2    <- y^2

      for (j in 1:g)
      {
        ### E-step: computing the moments

        #aux1<-MomN(mu[,j],sigma2[j],y)
        aux1<-CalMom(mu[,j],sigma2[j],nu,y,family)

        u1[cc==1] <- aux1$Ey1[cc==1]
        u2[cc==1] <- aux1$Ey2[cc==1]

        d1 <- dNormal(cc, y, mu[,j], sigma2[j])

        if(length(which(d1 == 0)) > 0) d1[which(d1 == 0)] <- .Machine$double.xmin

        d2 <- d.mixedN(cc, y, pii, mu, sigma2)

        if(length(which(d2 == 0)) > 0) d2[which(d2 == 0)] <- .Machine$double.xmin

        tal[,j] <- d1*pii[j]/d2

        # S1[,j] <- u1
        # S2[,j] <- u2


        ### M-step: update the parameters
        pii[j]    <- (1/n)*sum(tal[,j])
        sigma2[j] <- sum(tal[,j]*(u2+varphi[j]^2+mu1[,j]^2-2*u1*varphi[j]-2*u1*mu1[,j]+2*varphi[j]*mu1[,j])) / sum(tal[,j])
        varphi[j] <- sum(tal[,j]*(u1-mu1[,j]))/ sum(tal[,j])

        soma1     <- soma1+ t(x)%*%diag(tal[,j]/sigma2[j])%*%(u1-varphi[j])
        soma2     <- soma2+ t(x)%*%diag(tal[,j]/sigma2[j])%*%x
      }

      betas      <- solve(soma2)%*%soma1
      pii[g]     <- 1 - (sum(pii) - pii[g])

      zero.pos   <- NULL
      zero.pos   <- which(pii == 0)

      if(length(zero.pos) != 0)
      {
        pii[zero.pos]               <- 1e-10
        pii[which(pii == max(pii))] <- max(pii) - sum(pii[zero.pos])
      }


      if (pii[1]< 0.5)
      {
        # mu   <-cbind(mu[,2],mu[,1])
        # medj <- as.vector(c(medj[2],medj[1]))
        varphi <- as.vector(c(varphi[2], varphi[1]))
        pii    <- as.vector(c(pii[2], pii[1]))
        sigma2 <- as.vector(c(sigma2[2], sigma2[1]))
      }

      beta0<-sum(pii*varphi)

      for (k in 1:g)
      {
        mu1[,k]  <- x%*%betas
        mu[,k]   <- mu1[,k]+varphi[k]
        medj[k]  <- varphi[k]-beta0
      }

      Abetas     <- c(beta0,betas)
      teta       <- c(Abetas, medj, sigma2, pii)

      lk1        <- sum(log( d.mixedN(cc, y, pii, mu, sigma2) ))
      criterio   <- abs(lk1/lk-1)
      lk         <- lk1

      end.time        <- Sys.time()
      time.taken      <- end.time - start.time
    }

    if (criteria == TRUE)
    {
      cl  <- apply(tal, 1, which.max)
      icl <- 0

      for (j in 1:g) icl<-icl+sum(log(pii[j]*dNormal(cc, y, mu[,j], sigma2[j])))
      #icl <- -2*icl+(4*g-1)*log(n)
    }

    #### grafico ajuste
    ##hist(y, breaks = 40,probability=T,col="grey")
    ##xx=seq(min(y),max(y),(max(y)-min(y))/1000)
    ##lines(xx,d.mixedSN(xx, pii, mu, sigma2, shape),col="red")
    ######################
  }


  if (criteria == TRUE)
  {
    if((family == "Normal")) d <- p+1+ g*2 + (g-1) #Abeta + varphi, Sigma + pi
    else d  <- p+1+ g*2 + (g-1) +1 #Abeta + varphi+ Sigma +pi+nu
    aic     <- -2*lk + 2*d
    bic     <- -2*lk + log(n)*d
    edc     <- -2*lk + 0.2*sqrt(n)*d
    icl     <- -2*icl + log(n)*d

    obj.out <- list(Abetas = Abetas, medj=medj, sigma2 = sigma2, nu=nu,  pii = pii,  aic = aic, bic = bic, edc = edc, icl= icl, loglik=lk, iter = count, n = length(y), group = cl, time = time.taken)#, im.sdev=NULL)
  } else {
    obj.out <- list(Abetas = Abetas, medj=medj , sigma2 = sigma2, nu=nu, pii = pii, iter = count, n = length(y), group = apply(tal, 1, which.max), time = time.taken)#, im.sdev=NULL)
       }
  #if (group == FALSE) obj.out <- obj.out[-(length(obj.out)-1)]
  if (group == FALSE) obj.out <- obj.out[-(length(obj.out))]
  if (obs.prob == TRUE)
  {
    nam <- c()
    for (i in 1:ncol(tal))  {nam           <- c(nam,paste("Group ",i,sep=""))}
    if(ncol(tal) == 1) {dimnames(tal)[[2]] <- list(nam)}
    if(ncol(tal) > 1)  {dimnames(tal)[[2]] <- nam}
    obj.out$obs.prob <- tal
    if((ncol(tal) - 1) > 1) obj.out$obs.prob[,ncol(tal)]                <- 1 - rowSums(obj.out$obs.prob[,1:(ncol(tal)-1)])
    else obj.out$obs.prob[,ncol(tal)]                                   <- 1 - obj.out$obs.prob[,1]
    obj.out$obs.prob[which(obj.out$obs.prob[,ncol(tal)] < 0),ncol(tal)] <- 0.0000000000
    obj.out$obs.prob <- round(obj.out$obs.prob,10)
  }

  ################################################################################
  #Standard erro
  #Following the method of Basford et al. (1997)
  #Basford, K., Greenway, D., McLachlan, G. & Peel, D. (1997). Standard errors of fitted component
  #means of normal mixtures. Computational Statistics, 12(1)



  ################################################################################

  parameters       <- rbind(teta[1:(length(teta)-2)])
  parameters       <- c(parameters,teta[length(teta)])
  ttable           <- data.frame(as.numeric(parameters))

  #Row names of ttable
  namesrowAbetas <- c(); for(i in 1:length(Abetas)){namesrowAbetas[i] <- paste("beta",i-1,sep="")}
  namesrowMu     <- c(); for(i in 1:g){namesrowMu[i]<- paste("mu",i,sep="")}
  namesrowSigmas <- c(); for(i in 1:g){namesrowSigmas[i]<- paste("sigma",i,sep="")}
  rownames(ttable) <- c(namesrowAbetas,namesrowMu,namesrowSigmas,"pii1","nu")
  colnames(ttable) <- c("Estimate")


  if (criteria == TRUE){
   obj.out <- list(Abetas = Abetas, medj=medj, sigma2 = sigma2, nu=nu,  pii = pii, ttable=ttable,aic = aic, bic = bic, edc = edc, icl= icl, loglik=lk, iter = count, n = length(y), group = cl, time = time.taken)#, im.sdev=NULL)
   class(obj.out)   <- family
  } else {
   obj.out <- list(Abetas = Abetas, medj=medj , sigma2 = sigma2, nu=nu, pii = pii, ttable=ttable, iter = count, n = length(y), group = apply(tal, 1, which.max), time = time.taken)#, im.sdev=NULL)
   class(obj.out)   <- family
    }
  class(obj.out)   <- family
  return(obj.out)
}

#######   End the EM algorithm for mixture the normal and t-student with censored      ###############
######################################################################################################


