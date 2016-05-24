ridgeline <- function(alpha, mu1, mu2, Sigma1, Sigma2){
  Sigma1i <- solve(Sigma1)
  Sigma2i <- solve(Sigma2)
  alpha2 <- 1-alpha
  out <- solve(alpha*Sigma1i+alpha2*Sigma2i)%*%
      (alpha*(Sigma1i %*% mu1)+alpha2*(Sigma2i %*% mu2))
  out
}

dridgeline <- function(alpha=seq(0,1,0.001), prop,
                          mu1, mu2, Sigma1, Sigma2, showplot=FALSE, ...){
#  require(mvtnorm)
  out <- numeric(0)
  for (alpha1 in alpha){
    ralpha <- ridgeline(alpha1, mu1, mu2, Sigma1, Sigma2)
    if (dim(Sigma1)[1]==1)
      out <- c(out,prop*dnorm(ralpha,mu1,sqrt(Sigma1))+
             (1-prop)*dnorm(ralpha,mu2,sqrt(Sigma2)))
    else  
      out <- c(out,prop*dmvnorm(t(ralpha),mu1,Sigma1)+
             (1-prop)*dmvnorm(t(ralpha),mu2,Sigma2))
  }
  if (showplot)
    plot(alpha,out,type="l",ylab="density", ylim=c(0,max(out)), ...)
  invisible(out)
}

piridge <- function(alpha, mu1, mu2, Sigma1, Sigma2, showplot=FALSE){
#  require(mvtnorm)
  out <- numeric(0)
  for (alpha1 in alpha){
    ralpha <- ridgeline(alpha1, mu1, mu2, Sigma1, Sigma2)
    alpha2 <- 1-alpha1
    if (dim(Sigma1)[1]==1)
      out <- c(out,1/
            (1+alpha2*dnorm(ralpha,mu1,sqrt(Sigma1))/
             (alpha1*dnorm(ralpha,mu2,sqrt(Sigma2)))))
    else
      out <- c(out,1/
            (1+alpha2*dmvnorm(t(ralpha),mu1,Sigma1)/
             (alpha1*dmvnorm(t(ralpha),mu2,Sigma2))))
  }
  out[is.na(out) & alpha<1e-8] <- 0
  out[is.na(out) & alpha>1-1e-8] <- 1  
  if (showplot)
    plot(alpha,out,type="l",ylab="piridge")
  invisible(out)
}
    
piridge.zeroes <- function(prop, mu1, mu2, Sigma1, Sigma2, alphamin=0,
                          alphamax=1,by=0.001){
  alpha <- seq(alphamin,alphamax,by=by)
  la <- length(alpha)
#  print(Sigma1)
#  print(Sigma2)
  piridgezero <- piridge(alpha,mu1, mu2, Sigma1, Sigma2)-prop
#  print(piridgezero)
  piridgezero[is.na(piridgezero)] <- 0
  piridgesign <- sign(piridgezero)
  zeroplaces1 <- piridgesign==0
  zeroplaces2 <- abs(piridgesign[2:la]-piridgesign[1:(la-1)])==2
  number.zeros <- sum(zeroplaces1)+sum(zeroplaces2)
  estimated.roots <- sort(c(alpha[zeroplaces1],alpha[2:la][zeroplaces2]-by/2))
  out <- list(number.zeroes=number.zeros,estimated.roots=estimated.roots)
  out
}



# ridgelineplot can be "none", "matrix", "pairwise"
ridgeline.diagnosis <- function(propvector,muarray,Sigmaarray,
                                k=length(propvector),
                                ipairs="all", compute.ratio=TRUE,by=0.001,
                                ratiocutoff=NULL,ridgelineplot="matrix"){
#  require(prabclus)
  comat <- diag(k)
  if (compute.ratio)
    ratiomatrix <- diag(k)
  else
    ratiomatrix <- NULL
  pairlist <- ipairs
  if (identical(ipairs,"all")){
    pairlist <- list()
    m <- 1
    for (i in 1:(k-1))
      for (j in (i+1):k){
        pairlist[[m]] <- c(i,j)
        m <- m+1
      }
  }
#  print(pairlist)
  m <- length(pairlist)
  if (ridgelineplot=="matrix")
    par(mfrow=c(k-1,k))
  for (q in 1:m){
    if (q==1) ia <- TRUE
    else
      ia <- !(pairlist[[q]][1]==pairlist[[q-1]][1])
    i <- pairlist[[q]][1]
    j <- pairlist[[q]][2]
    if (ridgelineplot=="matrix" & ia)
      for (l in (1:i))
        plot(1,1,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
    propij <- propvector[i]/(propvector[i]+propvector[j])
    if (ridgelineplot!="none")
      dridgeline(seq(0,1,by=by),propij,muarray[,i],muarray[,j],
                   as.matrix(Sigmaarray[,,i]),as.matrix(Sigmaarray[,,j]),
                 showplot=TRUE,
                   main=paste("Components",i,"and",j))
#    cat(i," ",j," ",propij,muarray[,i],muarray[,j],Sigmaarray[,,i],
#        Sigmaarray[,,j],"\n")
    zerosij <- piridge.zeroes(propij,muarray[,i],muarray[,j],
                   as.matrix(Sigmaarray[,,i]),as.matrix(Sigmaarray[,,j])
                              ,by=by)
    if (ridgelineplot!="none")
      title(sub=paste("number of local optima=",zerosij$number.zeroes))
#    print(zerosij$number.zeroes)
    if (is.null(ratiocutoff) | !compute.ratio)
      comat[i,j] <- comat[j,i] <- zerosij$number.zeroes<2
    if (compute.ratio){
#      print(zerosij$estimated.roots)
       densij <- dridgeline(zerosij$estimated.roots,propij,
                             muarray[,i],muarray[,j],
                   as.matrix(Sigmaarray[,,i]),as.matrix(Sigmaarray[,,j]))
#      print(densij)
      if (min(c(densij[1],densij[zerosij$number.zeroes]))>0)
        ratiomatrix[i,j] <- ratiomatrix[j,i] <- min(densij)/
          min(c(densij[1],densij[zerosij$number.zeroes]))
      else
        ratiomatrix[i,j] <- 1
#      print(ratiomatrix)
#      print(ratiocutoff)
       if (!is.null(ratiocutoff))
          comat[i,j] <- comat[j,i] <- ratiomatrix[i,j]>=ratiocutoff
     } # if (compute.ratio)
  } # for q
# print(sum(comat))
  merged.clusters <- con.comp(comat)
  out <- list(merged.clusters=merged.clusters,connection.matrix=comat,
              ratiomatrix=ratiomatrix)
  out
}
    
# Misclassification prob. is majorized by exp(-bhat.dist.)
bhattacharyya.dist <- function(mu1, mu2, Sigma1, Sigma2){
  aggregatesigma <- (Sigma1+Sigma2)/2
  d1 <- mahalanobis(mu1,mu2,aggregatesigma)/8
  d2 <- log(det(as.matrix(aggregatesigma))/sqrt(det(as.matrix(Sigma1))*
                                                det(as.matrix(Sigma2))))/2
  out <- d1+d2
  out
}

bhattacharyya.matrix <- function(muarray,Sigmaarray,ipairs="all", 
                                 misclassification.bound=TRUE){
  k <- ncol(muarray)
#  print(muarray)
#  print(k)
  outmatrix <- matrix(NA,ncol=k,nrow=k)
  pairlist <- ipairs
  if (identical(ipairs,"all")){
    pairlist <- list()
    m <- 1
    for (i in 1:(k-1))
      for (j in (i+1):k){
        pairlist[[m]] <- c(i,j)
        m <- m+1
      }
  }
  m <- length(pairlist)
#  print(pairlist)
  for (q in 1:m){
      i <- pairlist[[q]][1]
      j <- pairlist[[q]][2]
#      print("bhatmatrix")
#      print(i)
#      print(j)
#      print(Sigmaarray)
      outmatrix[i,j] <- outmatrix[j,i] <-
        bhattacharyya.dist(muarray[,i],muarray[,j],as.matrix(Sigmaarray[,,i]),
                           as.matrix(Sigmaarray[,,j]))
  }
  if (misclassification.bound) outmatrix <- exp(-outmatrix)
  outmatrix
}

confusion <- function(z,pro,i,j,adjustprobs=FALSE){
  n <- nrow(z)
  probs <- pro[c(i,j)]
  zz <- z[,c(i,j)]
  if (adjustprobs){
    probs <- probs/sum(probs)
    zz <- zz/apply(zz,1,sum)
  }
  nclustering <- apply(zz,1,which.max)
  sum(zz[nclustering==1,2])/(n*probs[2])
}
  
# summary can be "max" or "mean"
# if symmetric=FALSE, [i,j] is estimated misclassification prob. of
# true j into i (i given true j)
zmisclassification.matrix <- function(z,pro=NULL,clustering=NULL,
                                      ipairs="all",symmetric=TRUE,
                                      stat="max"){
  k <- ncol(z)
  if (is.null(pro))
    pro <- apply(z,2,mean)
  if (is.null(clustering))
    clustering <- apply(z,1,which.max)
# print(str(z))
# print(k)
  outmatrix <- matrix(0,ncol=k,nrow=k)
  if (symmetric) outmatrix2 <- outmatrix
  pairlist <- ipairs
  if (identical(ipairs,"all")){
    pairlist <- list()
    m <- 1
    for (i in 1:(k-1))
      for (j in (i+1):k){
        pairlist[[m]] <- c(i,j)
        m <- m+1
      }
  }
  m <- length(pairlist)
# print(pairlist)
  for (q in 1:m){
      i <- pairlist[[q]][1]
      j <- pairlist[[q]][2]
#      print("bhatmatrix")
#      print(i)
#      print(j)
#      print(Sigmaarray)
      outmatrix[i,j] <- confusion(z,pro,i,j)
      outmatrix[j,i] <- confusion(z,pro,j,i)
      if (symmetric){
        val2 <- c(outmatrix[i,j],outmatrix[j,i])
        outmatrix2[i,j] <- outmatrix2[j,i] <- switch(stat,
                     max=max(val2),mean=mean(val2))
#        print(val2)
      }
  }
  if (symmetric) outmatrix <- outmatrix2  
  outmatrix
}


unimodal.ind <- function(y){
  ly <- length(y)
  t1 <- (y[2]>=y[1])
  out <- TRUE
  if (t1){
    q <-(2:(ly-1))[y[3:ly]<y[2:(ly-1)]]
    if (length(q)>0){
      k <- min((2:(ly-1))[y[3:ly]<y[2:(ly-1)]])
      out <- all(y[(k+1):ly]<=y[k:(ly-1)])
    }
    else
      out <- TRUE
  }
  else{
    k <- 2
#    k <- min((2:(ly-1))[y[3:ly]>y[2:(ly-1)]])    
    out <- all(y[(k+1):ly]<=y[k:(ly-1)])
  }
  out
}
    
dipp.tantrum <- function(xdata,d,M=100){
#  require(diptest)
  xfit <- density(xdata)
  n <- length(xdata)
  bw <- bwnew <- xfit$bw
  inc <- bw/10
  repeat{
    if (unimodal.ind(xfit$y))
      break
    else{
      bwnew <- bwnew+inc
      xfit <- density(xdata,bw=bwnew)
    }
  }
  dv <- numeric(0)
  for (i in 1:M){
    x.new <- rnorm(n, sample(xdata, size = n, replace = TRUE), bwnew)
    dv[i] <- dip(x.new)
  }
  p.value <- (1 + sum(dv >= d))/(1 + M)
  out <- list(p.value=p.value,bw=bwnew,dv=dv)
  out
}

diptest.multi <-  function(xdata,class,pvalue="uniform",M=100){
#  require(diptest)
  if (ncol(as.matrix(xdata))==1)
    xuni <- xdata
  else{
    dc <- discrcoord(xdata,class)
    xuni <- dc$proj[,1]
  }
  xd <- dip(xuni)
  if (pvalue=="uniform"){
    out <- dip.test(xd,simulate.p.value = FALSE)$p.value
#    out <- dip.pvalue(xd,nrow(as.matrix(xdata)))
#    print("diptest pvalue")
#    print(out)
  }
  if (pvalue=="tantrum")
   xd <- dip(xuni)
   out <- dipp.tantrum(xuni,xd,M=M)$p.value
  out
}

# Original taken from SLmisc, M. Kohl
# Methods: mclust, kmeans
# teststats: gap, bic (latter only for mclust)
# gaptestflex <-  function (xdata, class = rep(1, nrow(xdata)), M = 500,
#                              G=length(levels(as.factor(class))),
#                           gmethod="mclust", rotate=TRUE,
#                           teststat="gap", databic=NULL,
#                              modelNames=NULL, prior=NULL,
#                              control=emControl(), 
#                initialization=list(hcPairs=NULL, subset=NULL, noise=NULL), 
#                Vinv=NULL, warn=FALSE, plotnew=FALSE) 
# {
#     require(mclust)
#     if (!(length(class) == nrow(xdata))) 
#         stop("Length of class vector differs from nrow of data")
#     if (M <= 0) 
#         stop("'M' has to be a positive integer")
#     xdata <- as.matrix(xdata)
#     xdata <- scale(xdata, center = TRUE, scale = FALSE)
#     M <- trunc(M)
#     pw.dist <- function(x) {
#         sum(dist(x)/ncol(x))/2
#     }
#     if (teststat=="gap"){
#       temp1 <- log(sum(by(xdata, factor(class), pw.dist)))
#       temp11 <- log(sum(by(xdata, factor(rep(1,nrow(xdata))), pw.dist)))
#       datadiff <- temp11-temp1
#     }
#     else{
#       temp1 <- NULL
#       temp11 <- NULL
#     }
#     veigen <- svd(xdata)$v
#     if (rotate)
#       x1 <- crossprod(t(xdata), veigen)
#     else
#       x1 <- xdata
#     z1 <- matrix(NA, nrow = nrow(x1), ncol = ncol(x1))
#     tots <- tots1 <- vector(length = M)
#     for (k in 1:M) {
#         for (j in 1:ncol(x1)) {
#             min.x <- min(x1[, j])
#             max.x <- max(x1[, j])
#             z1[, j] <- runif(nrow(x1), min = min.x, max = max.x)
#         }
#         if (rotate)
#           z <- crossprod(t(z1), t(veigen))
#         else
#           z <- z1
#         if (gmethod=="mclust"){
#           zclus <- mclustBIC(z,G, modelNames, prior, control, 
#                initialization, 
#                Vinv, warn)
#           zs <- summary(zclus,z)
#           zclass <- zs$classification
#           if (plotnew){
#             pairs(z1,col=zclass)
#             print(summary(zclus,z))
#           }
#           if (teststat=="bic"){
#             tots[k] <- zs$bic[1]
#             zm1 <- mclustBIC(z,1, modelNames, prior, control, 
#                initialization, 
#                Vinv, warn)
# #            print(summary(zm1,z)$bic)
#             tots1[k] <- summary(zm1,z)$bic[1]
#           }
#         }
#         if (gmethod=="kmeans")
#           zclass <- kmeans(z,G)$cluster
#         if (teststat=="gap"){
#           tots[k] <- log(sum(by(z, factor(zclass), pw.dist)))
#           tots1[k] <- log(sum(by(z, factor(rep(1,nrow(xdata))), pw.dist)))
#         }
#     }
#     if (teststat=="bic"){
#       data1 <- mclustBIC(xdata,1,modelNames, prior, control, 
#                initialization, 
#                Vinv, warn)
#       if (is.null(databic)){
#         datam <- mclustBIC(xdata,G, modelNames, prior, control, 
#                initialization, Vinv, warn)
#         databic <- summary(datam,xdata)$bic[1]
#       }
#       datadiff <- databic-summary(data1,xdata)$bic[1]
#     }
#     tdiff <- tots1-tots
#     tdiffmean <- mean(tdiff)
#     tsd <- sqrt(1 + 1/M) * sd(tdiff)
#     p.value.sim <- (1 + sum(tdiff >= datadiff))/(1 + M)
#     p.value.sd <- 1-pnorm(datadiff,tdiffmean,tsd)
#     out <- list(gap=mean(tots) - temp1, sk=sqrt(1 + 1/M) * sd(tots),
#                 datadiff=datadiff, simulated.diff=tdiff, diffmean=tdiffmean,
#                 diffsd=tsd, gap1=mean(tots1)-temp11, simulatedt11=tots1,
#                 simulatedt1=tots,
#                 p.value.sim=p.value.sim,p.value.sd=p.value.sd)
#     return(out)
# }
# 
# # pvalue.mode is "simulation" or "sd"
# pairwise.gaptest <-  function(xdata, clustering, M = 200, ipairs="all",
#                               pvalue.mode="simulation", gmethod="mclust",
#                               teststat="gap",
#                               modelNames=NULL, prior=NULL,
#                              control=emControl(), 
#                initialization=list(hcPairs=NULL, subset=NULL, noise=NULL), 
#                Vinv=NULL, warn=FALSE, ...){
#   k <- max(clustering)
# #  print(k)
#   outmatrix <- matrix(0,ncol=k,nrow=k)
#   pairlist <- ipairs
#   if (identical(ipairs,"all")){
#     pairlist <- list()
#     m <- 1
#     for (i in 1:(k-1))
#       for (j in (i+1):k){
#         pairlist[[m]] <- c(i,j)
#         m <- m+1
#       }
#   }
#   m <- length(pairlist)
#   for (q in 1:m){
#       i <- pairlist[[q]][1]
#       j <- pairlist[[q]][2]
#       involvedx <- clustering==i | clustering==j
# #      pairs(xdata[involvedx,],col=clustering[involvedx])
#       gapobject <- 
#         gaptestflex(xdata[involvedx,], class = clustering[involvedx],
#                        M = M, gmethod=gmethod, teststat=teststat,
#                     G=2, modelNames=modelNames, prior=prior,
#                              control=control, 
#                initialization=initialization, 
#                Vinv=Vinv, warn=warn)
#       if(pvalue.mode=="simulation")
#         outmatrix[i,j] <- outmatrix[j,i] <- gapobject$p.value.sim
#       else
#         outmatrix[i,j] <- outmatrix[j,i] <- gapobject$p.value.sd
#   }
#   outmatrix
# }
# 
## already in zmisclassification.matrix, confusion
# aggregate can be "class1", "class2", "average", "max", "weighted", "both"
# assume 2 probs and a 2-d z-vector.
# misclassprob <- function(probs,z,aggregate="max",adjustprobs=TRUE){
#   n <- nrow(z)
#   if (adjustprobs){
#     probs <- probs/sum(probs)
#     z <- z/apply(z,1,sum)
#   }
#   classz <- apply(z,1,which.max)
#   mc1 <- mc2 <- 0
#   if (aggregate!="class1")
#     mc2 <- sum((classz==1)*z[,2])/(n*probs[2])
#   if (aggregate!="class2")
#     mc1 <- sum((classz==2)*z[,1])/(n*probs[1])
#   if (aggregate=="both")
#     mc <- c(mc1,mc2)
#   else
#     mc <- switch(aggregate,
#                class1=mc1,
#                class2=mc2,
#                average=mean(c(mc1,mc2)),
#                max=max(c(mc1,mc2)),
#                weighted=probs[1]*mc1+probs[2]*mc2)
#   mc
# }

mixdens <- function(modelName,data,parameters){
#  require(mvtnorm)
#  print("begin mixdens")
  G <- parameters$variance$G
  if (modelName %in% c("E","V")){
    out <- 0
    for (i in 1:G){
      if (modelName=="E")
        out <- out+parameters$pro[i]*dnorm(data,parameters$mean[i],
                                         sqrt(parameters$variance$sigmasq))
      if (modelName=="V")
        out <- out+parameters$pro[i]*dnorm(data,parameters$mean[i],
                                         sqrt(parameters$variance$sigmasq[i]))
    }
  }
  else{
    out <- 0
    if (G==1)
      out <- dmvnorm(data,parameters$mean,parameters$variance$sigma)
    else
      for (i in 1:G)
        out <- out+parameters$pro[i]*dmvnorm(data,
                                             as.matrix(parameters$mean)[,i],
                                         parameters$variance$sigma[,,i])
  }
#  print("end mixdens")
  out
}
    

                                         
                                           
    

extract.mixturepars <- function(mclustsum,compnumbers,noise=FALSE){
  d <- mclustsum$parameters$variance$d
  pcompnumbers <- compnumbers
  if (noise)
    pcompnumbers <- compnumbers+1
  parameters <- list()
  parameters$pro <- mclustsum$parameters$pro[pcompnumbers]/
        sum(mclustsum$parameters$pro[pcompnumbers])
  if (d>1)
    parameters$mean <- mclustsum$parameters$mean[,compnumbers]
  else
    parameters$mean <- mclustsum$parameters$mean[compnumbers]
  parameters$variance <- mclustsum$parameters$variance
  parameters$variance$G <- length(compnumbers)
  vm <- parameters$variance$modelName
  if (vm %in% c("V","VII"))
    parameters$variance$sigmasq <-
      mclustsum$parameters$variance$sigmasq[compnumbers]
  if (!(vm %in% c("E","V")))
    parameters$variance$sigma <-
      mclustsum$parameters$variance$sigma[,,compnumbers]
  if (vm=="VVV")
    parameters$variance$cholsigma <-
      mclustsum$parameters$variance$cholsigma[,,compnumbers]
  if (vm %in% c("EVI","VVI"))
    parameters$variance$shape <-
      mclustsum$parameters$variance$shape[,compnumbers]
  if (vm=="VEV")
    parameters$variance$orientation <-
      mclustsum$parameters$variance$orientation[,,compnumbers]
  parameters
}

## Documented until here. 



mixpredictive <- function(xdata, Gcomp, Gmix, M=50, ...){
  xdata <- as.matrix(xdata)
  n <- nrow(xdata)
#  print(str(xdata))
#  print(n)
  p <- ncol(xdata)
  nf <- c(floor(n/2),n-floor(n/2))
  indvec <- clusterings <- sclusterings <- jclusterings <-
    classifications <- list()
  prederr <- numeric(0)
  for (l in 1:M){
    nperm <- sample(n,n)
#    cat("Mixpredrun ",l,"\n")
    indvec[[l]] <- list()
    indvec[[l]][[1]] <- nperm[1:nf[1]]
    indvec[[l]][[2]] <- nperm[(nf[1]+1):n]
    for (i in 1:2){
#    print("pred mclust")
      clusterings[[i]] <- mclustBIC(xdata[indvec[[l]][[i]],],G=Gcomp,...)
      sclusterings[[i]] <- summary(clusterings[[i]],xdata[indvec[[l]][[i]],])
#    cvmodels[[i]] <- mclustModel(xdata[indvec[[l]][[i]],],clusterings[[i]])
      jclusterings[[i]] <- mergenormals(xdata[indvec[[l]][[i]],],
                                      sclusterings[[i]],
                                      method="demp", numberstop=Gmix)
#      pairs(xdata[indvec[[l]][[i]],],col=jclusterings[[i]]$clustering,main=Gmix)
#    print(Gmix)
#    print(jclusterings[[i]]$clusternumbers)
      j <- 3-i
      densvalues <- list()
      densvalues[[j]] <- matrix(0,nrow=nf[j],ncol=Gmix)
      for (k in 1:Gmix){
#      print(j)
#      print(k)
#      print(jclusterings[[i]]$clusternumbers[k])
#      print(jclusterings[[i]]$probs)
#      print(indvec[[l]][[j]])
#      print(jclusterings[[i]]$parameters[[k]])
        densvalues[[j]][,k] <-
          jclusterings[[i]]$probs[jclusterings[[i]]$clusternumbers[k]]*mixdens(
                    jclusterings[[i]]$parameters[[k]]$variance$modelName,
                    xdata[indvec[[l]][[j]],],
                    parameters=jclusterings[[i]]$parameters[[k]])
#      print(densvalues[[j]][,k])
      } # for k
      classifications[[j]] <- apply(densvalues[[j]],1,which.max)
#    print(classifications[[j]])
    } # for i
    ps <- matrix(0,nrow=2,ncol=Gmix)
    for (i in 1:2){
      ctable <- table(jclusterings[[i]]$clustering,classifications[[i]])
      for (k in 1:Gmix){
        ps[i,k] <- sum(ctable[k,]^2-ctable[k,]) 
        nik <- sum(jclusterings[[i]]$clustering==k)
        if (nik>1)
          ps[i,k] <- ps[i,k]/(nik*(nik-1))
        else
          ps[i,k] <- 1
#        if (nik>1){
#          for (j1 in (1:(nf[i]-1))[jclusterings[[i]]$clustering[1:(nf[i]-1)]==
#                 k]){
##        print(j1)
##        print(nf[i])
#            for (j2 in (j1+1):nf[i])
#              if (jclusterings[[i]]$clustering[j2]==k)
#                ps[i,k] <- ps[i,k]+(classifications[[i]][j1]==
#                              classifications[[i]][j2])
#          } # for j1
#          ps[i,k] <- 2*ps[i,k]/(nik*(nik-1))
#        } # if nik>0
      } # for k
#      qq <- which.min(ps[i,])
#      pairs(xdata[indvec[[l]][[i]],],
#            col=(jclusterings[[i]]$clustering==qq)+
#            (jclusterings[[i]]$clustering!=qq)*
#            (jclusterings[[i]]$clustering+1),main=Gmix)
#      title(sub=min(ps[i,]))

    } # for i
    prederr[l] <- mean(c(min(ps[1,]),min(ps[2,])))
#    title(sub=prederr[l])
#    print(prederr[l])    
  } # for l
  out <- list(predcorr=prederr,mean.pred=mean(prederr))
#  print(ps)
  out
}

prediction.strength <- function(xdata, Gmin=2, Gmax=10,M=50,
                                clustermethod=kmeansCBI,
                                classification="centroid",
                                cutoff=0.8,nnk=1,
                                distances=inherits(xdata,"dist"),
                                count=FALSE,...){
#  require(cluster)
#  require(class)
  xdata <- as.matrix(xdata)
  n <- nrow(xdata)
#  print(str(xdata))
#  print(n)
  nf <- c(floor(n/2),n-floor(n/2))
  indvec <- clcenters <- clusterings <- jclusterings <-
    classifications <- list()
  corrpred <- list()
  for (k in Gmin:Gmax){
    if (count)
      cat(k," clusters\n")
    corrpred[[k]] <- numeric(0)
    for (l in 1:M){
      nperm <- sample(n,n)
      if (count)
        cat(" Run ",l,"\n")
      indvec[[l]] <- list()
      indvec[[l]][[1]] <- nperm[1:nf[1]]
      indvec[[l]][[2]] <- nperm[(nf[1]+1):n]
      for (i in 1:2){
        if (distances)  
          clusterings[[i]] <- clustermethod(as.dist(xdata[indvec[[l]][[i]],
                                                    indvec[[l]][[i]]]),k,...)
        else
          clusterings[[i]] <- clustermethod(xdata[indvec[[l]][[i]],],k,...)
        jclusterings[[i]] <- rep(-1,n)
        jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$partition
        centroids <- NULL
        if(classification=="centroid"){
          if(identical(clustermethod,kmeansCBI))
            centroids <- clusterings[[i]]$result$centers
          if(identical(clustermethod,claraCBI))
            centroids <- clusterings[[i]]$result$medoids
        }
        j <- 3-i
        if (distances)
          classifications[[j]] <- classifdist(as.dist(xdata),jclusterings[[i]],
                                          method=classification,
                           centroids=centroids,nnk=nnk)[indvec[[l]][[j]]]
        else
          classifications[[j]] <- classifnp(xdata,jclusterings[[i]],
                                          method=classification,
                           centroids=centroids,nnk=nnk)[indvec[[l]][[j]]]
#        print(classifications[[j]])
      } # for i
#      browser()
      ps <- matrix(0,nrow=2,ncol=k)
      for (i in 1:2){
        ctable <- table(clusterings[[i]]$partition,classifications[[i]])
        for (kk in 1:k){
#          for (kkk in 1:k) 
          ps[i,kk] <- sum(ctable[kk,]^2-ctable[kk,]) 
          cpik <- clusterings[[i]]$partition==kk
          nik <- sum(cpik)
          if (nik>1)
            ps[i,kk] <- ps[i,kk]/(nik*(nik-1))
          else
            ps[i,kk] <- 1
 
#           cpiki <- (1:nf[i])[cpik]
#           nik <- sum(cpik)
#           if (nik>1){
# #            print(kk)
#             j11 <- 1
#             if (j11<=nik-1){
#               j1 <- cpiki[j11]
# #            for (j1 in (1:(nf[i]-1))[clusterings[[i]]$partition[1:(nf[i]-1)]==kk]){
#               print(j1)
# #        print(nf[i])
# #              for (j2 in (j1+1):nf[i])
#               j21 <- j11+1
#               if(j21<=nik){
#                   j2 <- cpiki[j21]
#                   if (clusterings[[i]]$partition[j2]==kk)
#                     ps[i,kk] <- ps[i,kk]+(classifications[[i]][j1]==
#                               classifications[[i]][j2])
#                   j21 <- j21+1
#               } # for j2
#               j11 <- j11+1
#             } # for j1
#             ps[i,kk] <- 2*ps[i,kk]/(nik*(nik-1))
#           } # if nik>0

          
         } # for kk
#      qq <- which.min(ps[i,])
#      pairs(xdata[indvec[[l]][[i]],],
#            col=(jclusterings[[i]]$clustering==qq)+
#            (jclusterings[[i]]$clustering!=qq)*
#            (jclusterings[[i]]$clustering+1),main=Gmix)
#      title(sub=min(ps[i,]))

      } # for i
      corrpred[[k]][l] <- mean(c(min(ps[1,]),min(ps[2,])))
#      browser()
#    title(sub=corrpred[l])
#    print(corrpred[l])
    } # for l
  } # for k
  mean.pred <- numeric(0)
  if (Gmin>1)
    mean.pred <- c(1)
  if (Gmin>2)
    mean.pred <- c(mean.pred,rep(NA,Gmin-2))
  for (k in Gmin:Gmax)
    mean.pred <- c(mean.pred,mean(corrpred[[k]]))
  optimalk <- max(which(mean.pred>cutoff))
  out <- list(predcorr=corrpred,mean.pred=mean.pred,optimalk=optimalk,
                cutoff=cutoff,method=clusterings[[1]]$clustermethod,
              Gmax=Gmax,M=M)
#  print(ps)
  class(out) <- "predstr"
  out
}


# prediction.strength <- function(xdata, Gmin=2, Gmax=10,M=50,
#                                 clustermethod=kmeansCBI,
#                                 classification="centroid",
#                                 cutoff=0.8,nnk=1,...){
#   require(cluster)
#   require(class)
#   xdata <- as.matrix(xdata)
#   n <- nrow(xdata)
# #  print(str(xdata))
# #  print(n)
#   p <- ncol(xdata)
#   nf <- c(floor(n/2),n-floor(n/2))
#   indvec <- clcenters <- clusterings <- jclusterings <-
#     classifications <- list()
#   prederr <- list()
#   for (k in Gmin:Gmax){
#     prederr[[k]] <- numeric(0)
#     for (l in 1:M){
#       nperm <- sample(n,n)
# #    cat("Mixpredrun ",l,"\n")
#       indvec[[l]] <- list()
#       indvec[[l]][[1]] <- nperm[1:nf[1]]
#       indvec[[l]][[2]] <- nperm[(nf[1]+1):n]
#       for (i in 1:2){
#         clusterings[[i]] <- clustermethod(xdata[indvec[[l]][[i]],],k,...)
#         jclusterings[[i]] <- rep(-1,n)
#         jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$partition
#         centroids <- NULL
#         if(classification=="centroid"){
#           if(identical(clustermethod,kmeansCBI))
#             centroids <- clusterings[[i]]$result$centers
#           if(identical(clustermethod,claraCBI))
#             centroids <- clusterings[[i]]$result$medoids
#         }
#         j <- 3-i
#         classifications[[j]] <- classifnp(xdata,jclusterings[[i]],
#                                           method=classification,
#                            centroids=centroids,nnk=nnk)[indvec[[l]][[j]]]
# #        print(classifications[[j]])
#       } # for i
#       ps <- matrix(0,nrow=2,ncol=k)
#       for (i in 1:2){
#         for (kk in 1:k){
#           nik <- sum(clusterings[[i]]$partition==kk)
#           if (nik>1){
# #            print(kk)
#             for (j1 in (1:(nf[i]-1))[clusterings[[i]]$partition[1:(nf[i]-1)]==kk]){
# #        print(j1)
# #        print(nf[i])
#               for (j2 in (j1+1):nf[i])
#                 if (clusterings[[i]]$partition[j2]==kk)
#                   ps[i,kk] <- ps[i,kk]+(classifications[[i]][j1]==
#                               classifications[[i]][j2])
#             } # for j1
#             ps[i,kk] <- 2*ps[i,kk]/(nik*(nik-1))
#           } # if nik>0
#         } # for kk
# #      qq <- which.min(ps[i,])
# #      pairs(xdata[indvec[[l]][[i]],],
# #            col=(jclusterings[[i]]$clustering==qq)+
# #            (jclusterings[[i]]$clustering!=qq)*
# #            (jclusterings[[i]]$clustering+1),main=Gmix)
# #      title(sub=min(ps[i,]))
# 
#       } # for i
#       prederr[[k]][l] <- mean(c(min(ps[1,]),min(ps[2,])))
# #      browser()
# #    title(sub=prederr[l])
# #    print(prederr[l])    
#     } # for l
#   } # for k
#   mean.pred <- numeric(0)
#   if (Gmin>1)
#     mean.pred <- c(1)
#   if (Gmin>2)
#     mean.pred <- c(mean.pred,rep(NA,Gmin-2))
#   for (k in Gmin:Gmax)
#     mean.pred <- c(mean.pred,mean(prederr[[k]]))
#   optimalk <- max(which(mean.pred>cutoff))
#   out <- list(predcorr=prederr,mean.pred=mean.pred,optimalk=optimalk,
#                 cutoff=cutoff,method=clusterings[[1]]$clustermethod,
#               Gmax=Gmax,M=M)
# #  print(ps)
#   class(out) <- "predstr"
#   out
# }

print.predstr <- function(x, ...){
  cat("Prediction strength \n")
  cat("Clustering method: ",x$method,"\n")
  cat("Maximum number of clusters: ",x$Gmax,"\n")
  cat("Resampled data sets: ",x$M,"\n")
  cat("Mean pred.str. for numbers of clusters: ",x$mean.pred,"\n")
  cat("Cutoff value: ", x$cutoff,"\n")
  cat("Largest number of clusters better than cutoff: ",x$optimalk,"\n")
}

# join.predictive <- function(xdata,mclustsummary,cutoff=0.8,M=50,...){
#   G <- mclustsummary$G
#   predvalues <- rep(1,G)
#   for (q in G:2)
#     predvalues[q] <- mixpredictive(as.matrix(xdata),G,q,M,...)$mean.pred
#   Gopt <- max((1:G)[predvalues>cutoff])
# #  print(G)
# #  print(predvalues)
# #  print(Gopt)
#   joinnormals <- mergenormals(xdata,mclustsummary,
#                               method="demp", numberstop=Gopt)
#   joinnormals$predvalues <- predvalues
#   joinnormals
# }
# 
# 
# nmergecluster <- function(xdata,G=1:9,noisek=0,modelNames=NULL,
#                           method="predictive",cutoff,...){
#   require(mclust)
#   if (noisek==0){
#     xm <- mclustBIC(xdata,G,modelNames)
#     sxm <- summary(xm,xdata)
#     if (method=="predictive")
#       out <- join.predictive(xdata,sxm,cutoff,...)
#     else
#       out <- mergenormals(xdata,sxm,method=method,cutoff=cutoff,...)
#   }
#   out$mclustsummary <- sxm
#   out
# }
# 
#   
    
mergeparameters <- function(xdata, j1, j2, probs,
                         muarray,Sigmaarray, z){
  probs[j1] <- probs[j1]+probs[j2]
  z[,j1] <- z[,j1]+z[,j2]
  covwt <- cov.wt(as.matrix(xdata),wt=z[,j1],method="ML")
  muarray[,j1] <- covwt$center
#  print(muarray)
#  print(Sigmaarray)
  Sigmaarray[,,j1] <- covwt$cov
  out <- list(probs=probs,muarray=muarray,Sigmaarray=Sigmaarray,z=z)
  out
}

# Implemented methods: "bhat", "ridge.uni",
# "ridge.ratio", "demp", "dipuni",
# "diptantrum", "predictive"
# noise is ignored if cluster 0 is noise.
# numberstop: fixed number of clusters to achieve by merging.
mergenormals <- function(xdata, mclustsummary=NULL, 
                         clustering, probs, muarray, Sigmaarray, z,
                         method=NULL, cutoff=NULL, by=0.005,
                         numberstop=NULL, renumber=TRUE, M=50, ...){
  if (is.null(cutoff)){
    cutoff <- 1
    if (method=="bhat")
      cutoff <- 0.1
    if (method=="ridge.ratio")
      cutoff <- 0.2
    if (method %in% c("dipuni","diptantrum"))
      cutoff <- 0.05
    if (method == "demp")
      cutoff <- 0.025
    if (method=="predictive")
      cutoff <- 0.75
  }    
  if (method=="predictive"){
    G <- mclustsummary$G
    predvalues <- rep(1,G)
    for (q in G:2)
      predvalues[q] <-
        mixpredictive(as.matrix(xdata)[mclustsummary$classification>0,],
                      G,q,M,...)$mean.pred
    Gopt <- max((1:G)[predvalues>cutoff])
#  print(G)
#  print(predvalues)
#  print(Gopt)
    joinnormals <- mergenormals(xdata,mclustsummary,
                              method="demp", numberstop=Gopt,
                                renumber=renumber)
    joinnormals$predvalues <- predvalues
    out <- joinnormals
    out$method <- "predictive"
    out$cutoff <- cutoff
  } # predictive
  else{
    predvalues <- NULL
    if (!is.null(mclustsummary)){
      probs <- mclustsummary$parameters$pro
      muarray <- mclustsummary$parameters$mean
      clustering <- mclustsummary$classification
      Sigmaarray <- mclustsummary$parameters$variance$sigma
      z <- mclustsummary$z
      ncomp <- mclustsummary$G
    }
    else
      ncomp <- max(clustering)
    if (!is.null(mclustsummary))
      if (mclustsummary$d==1){
        Sigmaarray <- array(0,dim=c(1,1,ncomp))
        if (mclustsummary$modelName=="E"){
          for (i in 1:ncomp)
            Sigmaarray[,,i] <- mclustsummary$parameters$variance$sigmasq
        }
        else{
          for (i in 1:ncomp)
            Sigmaarray[,,i] <- mclustsummary$parameters$variance$sigmasq[i]
        }
        muarray <- t(as.matrix(muarray))
      }
    noise <- (min(clustering)==0)
    if (noise){
      probs <- probs[c(2:(ncomp+1),1)]
      z <- z[,c(2:(ncomp+1),1)]
    }      
  #  print(mclustsummary$modelName)
  #  print(probs)
  # print(method)
  #  print(Sigmaarray)
  #  print(sum(mclustsummary$classification==0))
    valuemerged <- numeric(0)
    defunct.components <- numeric(0)
    curr.clustering <- clustering
    curr.clusters <- 1:ncomp
    mergedtonumbers <- 1:ncomp
    updatepairs <- "all"
    nmaxc <- FALSE
    if (is.numeric(numberstop)){
      cutoff <- 0
      nmaxc <- numberstop==ncomp
    }
#    print(updatepairs)
#    print("before repeat")
    repeat{
      if (length(updatepairs)!=0){
        if (method=="bhat" || method=="ridge.uni"){
  #        print(updatepairs)
  #        print(muarray)
  #        print(Sigmaarray)
          new.bhatmatrix <- bhattacharyya.matrix(muarray,Sigmaarray,
                                                 ipairs=updatepairs)
  #        print(new.bhatmatrix)
          if (!identical(updatepairs,"all")){
            m <- length(updatepairs)
            for (q in 1:m)
              bhatmatrix[updatepairs[[q]][1],updatepairs[[q]][2]] <-
                new.bhatmatrix[updatepairs[[q]][1],updatepairs[[q]][2]]
          }
          else{
            bhatmatrix <- new.bhatmatrix
  #          print("else")
  #          print(bhatmatrix)
          }
        } # if bhattacharyya || ridgeline.unimodal
        if (method=="bhat"){
          if (identical(updatepairs,"all"))
            orig.decisionmatrix <- bhatmatrix
          joinmatrix <- matrix(1,ncol=ncol(bhatmatrix),nrow=nrow(bhatmatrix))
          joinmatrix[bhatmatrix<cutoff] <- 0
          joinmatrix[defunct.components,] <- joinmatrix[,defunct.components] <- 0
          bhatmatrix[defunct.components,] <- bhatmatrix[,defunct.components] <- 0
          decisionmatrix <- bhatmatrix
  #        print("dec")
  #        print(decisionmatrix)
  #        print(joinmatrix)
        } # if bhattachryya
        if (method=="ridge.uni"){
          joinmatrix <- ridgeline.diagnosis(probs, muarray, Sigmaarray,
                                            ipairs=updatepairs,k=ncomp, 
                                            compute.ratio=FALSE,by=by,
                                  ridgelineplot="none",...)$connection.matrix
          if (identical(updatepairs,"all"))
            orig.decisionmatrix <- joinmatrix
          joinmatrix[defunct.components,] <- joinmatrix[,defunct.components] <- 0
          bhatmatrix[defunct.components,] <- bhatmatrix[,defunct.components] <- 0
          decisionmatrix <- bhatmatrix
        } # if ridgeline.unimodal
        if (method=="dipuni" || method=="diptantrum"){
#          print("updatepairs")
#          print(updatepairs)
#          exists("joinmatrix")
          new.joinmatrix <- ridgeline.diagnosis(probs, muarray, Sigmaarray,
                                            ipairs=updatepairs,k=ncomp, 
                                            compute.ratio=TRUE,by=by,
                                  ridgelineplot="none",...)$ratiomatrix
          if (!identical(updatepairs,"all") & exists("joinmatrix")){
#       print("if not")     
#       exists("joinmatrix")
           m <- length(updatepairs)
            for (q in 1:m)
              joinmatrix[updatepairs[[q]][1],updatepairs[[q]][2]] <-
                new.joinmatrix[updatepairs[[q]][1],updatepairs[[q]][2]]
          }
          else
            joinmatrix <- new.joinmatrix
#         print(defunct.components)
 #        exists("joinmatrix")
          if (identical(updatepairs,"all"))
            orig.decisionmatrix <- joinmatrix
          joinmatrix[defunct.components,] <- joinmatrix[,defunct.components] <- 0
          joinmatrix[lower.tri(joinmatrix,diag=TRUE)] <- 0
          decisionmatrix <- joinmatrix
  #        print(joinmatrix)
          jpairs <- which(joinmatrix==1,arr.ind=TRUE)
          if (nrow(jpairs)>0){
            jpair <- jpairs[1,]
  #          print(jpairs)
            if (nrow(jpairs)>1){
              involvedx <- clustering== jpairs[1,1] | clustering==jpairs[1,2]
              if (method=="dipuni")
                jbhat <- diptest.multi(as.matrix(xdata)[involvedx,],
                                       clustering[involvedx])
              if (method=="diptantrum")
                jbhat <- diptest.multi(as.matrix(xdata)[involvedx,],
                                       clustering[involvedx],
                                       pvalue="tantrum")
              for (i in 2:nrow(jpairs)){
                involvedx <- clustering== jpairs[i,1] | clustering==jpairs[i,2]
                if (method=="dipuni")
                  dm <- diptest.multi(as.matrix(xdata)[involvedx,],
                                      clustering[involvedx])
                if (method=="diptantrum")
                  dm <- diptest.multi(as.matrix(xdata[involvedx,]),
                                      clustering[involvedx],
                                      pvalue="tantrum")
                if (dm>jbhat){
                  jbhat <- dm
                  jpair <- jpairs[i,]
                }
              }
              decisionmatrix[jpair[1],jpair[2]] <- 1+jbhat
            } # if (nrow(jpairs)>1)
          } # if (nrow(jpairs)>0)
          if (max(decisionmatrix)<=1){
            jpair <- which(decisionmatrix==max(decisionmatrix),arr.ind=TRUE)[1,]
            involvedx <- clustering== jpair[1] | clustering==jpair[2]
            if (method=="dipuni")
              jbhat <- diptest.multi(as.matrix(xdata)[involvedx,],
                                     clustering[involvedx])
            if (method=="diptantrum")
              jbhat <- diptest.multi(as.matrix(xdata)[involvedx,],
                                     clustering[involvedx],
                                     pvalue="tantrum")
          }
          if (jbhat<cutoff)
            joinmatrix <- matrix(0,ncol=ncomp,nrow=ncomp)
  #        print("jbhat")
  #        print(jbhat)
          valuemerged <- c(valuemerged,jbhat)
  #        print(valuemerged)
        } # if dip
        if (method=="ridge.ratio"){
#          print(updatepairs)
          jrd <- ridgeline.diagnosis(probs, muarray, Sigmaarray,
                                            ipairs=updatepairs,k=ncomp,
                                            compute.ratio=TRUE,by=by,
                                            ratiocutoff=cutoff,
                                  ridgelineplot="none",...)
 #         print(jrd)
#          print(defunct.components)
          if (!identical(updatepairs,"all")){
            m <- length(updatepairs)
            for (q in 1:m){
              joinmatrix[updatepairs[[q]][1],updatepairs[[q]][2]] <-
                jrd$connection.matrix[updatepairs[[q]][1],updatepairs[[q]][2]]
              decisionmatrix[updatepairs[[q]][1],updatepairs[[q]][2]] <-
                jrd$ratiomatrix[updatepairs[[q]][1],updatepairs[[q]][2]]
            }
          }
          else{
            joinmatrix <- jrd$connection.matrix
            orig.decisionmatrix <- decisionmatrix <- jrd$ratiomatrix
          }
#          print(decisionmatrix)
          joinmatrix[defunct.components,] <- joinmatrix[,defunct.components] <- 0
          decisionmatrix[defunct.components,] <-
            decisionmatrix[,defunct.components] <- 0
        } # if ridgeline.densityratio
        if (method=="demp"){
 # print(updatepairs)
 # print(curr.clustering)
          if(ncomp<ncol(z)){
            zmm <- z[,1:ncomp]
            probsm <- probs[1:ncomp]
          }
          else{
            zmm <- z
            probsm <- probs
          }
          new.zmcmatrix <- zmisclassification.matrix(zmm,probs,
                                                     curr.clustering, 
                                                     ipairs=updatepairs)
 #  print(dim(new.zmcmatrix))
          if (!identical(updatepairs,"all")){
            m <- length(updatepairs)
            for (q in 1:m)
              zmcmatrix[updatepairs[[q]][1],updatepairs[[q]][2]] <-
                new.zmcmatrix[updatepairs[[q]][1],updatepairs[[q]][2]]
          }        
          else
            orig.decisionmatrix <- zmcmatrix <- new.zmcmatrix
          joinmatrix <- matrix(1,ncol=ncol(zmcmatrix),nrow=nrow(zmcmatrix))
          joinmatrix[zmcmatrix<cutoff] <- 0
          joinmatrix[defunct.components,] <- joinmatrix[,defunct.components] <- 0
          zmcmatrix[defunct.components,] <- zmcmatrix[,defunct.components] <- 0
          decisionmatrix <- zmcmatrix
  #        print(zmcmatrix)
        } # if demp
 #     print("error here")
        joinmatrix[lower.tri(joinmatrix,diag=TRUE)] <- 0
        decisionmatrix[lower.tri(decisionmatrix,diag=TRUE)] <- 0
 #     print("no")
      } # if length updatepairs!=0
  #    print(numberstop)
  #    print(updatepairs)
  #    print(joinmatrix)
  #    print(nmaxc)
      if (all(joinmatrix==0) || length(updatepairs)==0 || nmaxc){
        newclustering <- curr.clustering
        if (method=="bhat" || method=="ridge.ratio"
            || method=="demp")
          valuemerged <- c(valuemerged, max(decisionmatrix))
        break
      }
      else{
        if (method=="bhat" || method=="ridge.ratio"
            || method=="demp"
            || method=="dipuni" || method=="diptantrum"){
  #        valuemerged <- c(valuemerged, max(decisionmatrix))
          jpair <- which(decisionmatrix==max(decisionmatrix),arr.ind=TRUE)[1,]
#  cat("jpair=",jpair,"\n")
#  print(decisionmatrix)
        }
        if (method=="bhat" || method=="ridge.ratio"
            || method=="demp")
          valuemerged <- c(valuemerged, max(decisionmatrix))
  #      print(jpair)
        if (method=="ridge.uni"){
          jpairs <- which(joinmatrix==1,arr.ind=TRUE)
          jpair <- jpairs[1,]
          if (nrow(jpairs)>1){
            jbhat <- decisionmatrix[jpairs[1,1],jpairs[1,2]]
            for (i in 2:nrow(jpairs))
              if (decisionmatrix[jpairs[i,1],jpairs[i,2]]>jbhat){
                jbhat <- decisionmatrix[jpairs[i,1],jpairs[i,2]]
                jpair <- jpairs[i,]
              }
          }
        } # if ridgeline.unimodal        
        curr.clustering[curr.clustering==jpair[2]] <- jpair[1]
        defunct.components <- c(defunct.components,jpair[2])
        mergedtonumbers[mergedtonumbers==jpair[2]] <- jpair[1]
  #      curr.clusters <- (1:ncomp)[(1:ncomp)%in% curr.clustering]
        curr.clusters <- curr.clusters[curr.clusters!=jpair[2]]
        updatepairs <- list()
        k <- 1
        for (i in curr.clusters[curr.clusters<jpair[1]]){
          updatepairs[[k]] <- c(i,jpair[1])
          k <- k+1
        }
        for (i in curr.clusters[curr.clusters>jpair[1]]){
          updatepairs[[k]] <- c(jpair[1],i)
          k <- k+1
        }
  #       joinmatrix[jpair[1],jpair[2]] <- decisionmatrix[jpair[1],jpair[2]] <- 0
  #       j2index <- (1:ncomp)[joinmatrix[jpair[2],]==1 | joinmatrix[,jpair[2]]==1]
  #       j12index <- j2index[j2index>jpair[1]]
  #       j21index <- j2index[j2index<jpair[1]]
  #       joinmatrix[jpair[1],j12index] <- 1
  #       joinmatrix[j21index,jpair[1]] <- 1
  #       
  #       joinmatrix[jpair[2],] <- 0
  #       joinmatrix[,jpair[2]] <- 0
  #       decisionmatrix[jpair[2],] <- 0
  #       decisionmatrix[,jpair[2]] <- 0
  #      if (!all(c(joinmatrix[jpair[1],],joinmatrix[,jpair[1]])==0)){
  #      print(jpair)
  #        print(Sigmaarray)
        mp <- mergeparameters(xdata, jpair[1], jpair[2], probs,
                           muarray,Sigmaarray, z)
        probs <- mp$probs
        muarray <- mp$muarray
        Sigmaarray <- mp$Sigmaarray
        z <- mp$z
  #      if (method=="bhat" || method=="ridge.uni"){
  #        mb <- merge.bhat(jpair[1], jpair[2], decisionmatrix, probs,
  #                         muarray, Sigmaarray, curr.clusters)
  #        decisionmatrix <- mb
  #      }
  #      if (method!="bhat")
  #        define.merge.method()
        if (is.numeric(numberstop))
          if (numberstop==ncomp-length(defunct.components)){
  #          cat("join normals ",ncomp, " defunct ",defunct.components,"\n")
            newclustering <- curr.clustering
            break
          }
      } # else (!all(joinmatrix==0) || length(updatepairs)==0)
    } # repeat
    parameters <- NULL
    if (!is.null(mclustsummary)){
      parameters <- list()
      Gmix <- length(curr.clusters)
      Gcomp <- length(mergedtonumbers)
      for (k in 1:Gmix){
        nk <- curr.clusters[k]
        compk <- (1:Gcomp)[mergedtonumbers==nk]
        parameters[[k]] <- extract.mixturepars(mclustsummary,compk,noise)
      }
    }
    if (renumber){
      renumcl <- newclustering
      cln <- 1
      for (i in 1:max(newclustering)){
        if(sum(newclustering==i)>0){
          renumcl[newclustering==i] <- cln
          mergedtonumbers[mergedtonumbers==i] <- cln
          cln <- cln+1
        }
      }
      newclustering <- renumcl
    }
    if (length(valuemerged)==length(mergedtonumbers))
      valuemerged <- valuemerged[-length(valuemerged)]
    out <- list(clustering=newclustering,
                clusternumbers=curr.clusters,
                defunct.components=defunct.components,
                valuemerged=valuemerged,
                mergedtonumbers=mergedtonumbers,
                parameters=parameters,
                predvalues=predvalues,
                orig.decisionmatrix=orig.decisionmatrix,
                new.decisionmatrix=decisionmatrix+t(decisionmatrix)+
                 diag(diag(decisionmatrix)),
                probs=probs, muarray=muarray, Sigmaarray=Sigmaarray,
                z=z, noise=noise,
                method=method, cutoff=cutoff)
  } # method!="predictive"
  class(out) <- "mergenorm"
  out
}

summary.mergenorm <- function(object, ...){
  out <- list(clustering=object$clustering,
              onc=length(object$mergedtonumbers),
              mnc=length(object$clusternumbers),
              clusternumbers=object$clusternumbers,
                defunct.components=object$defunct.components,
                valuemerged=object$valuemerged,
                mergedtonumbers=object$mergedtonumbers,
                predvalues=object$predvalues,
                probs=object$probs[object$clusternumbers],
              muarray=object$muarray[,object$clusternumbers],
              Sigmaarray=object$Sigmaarray[,,object$clusternumbers],
                z=object$z[,object$clusternumbers],
              noise=object$noise, 
                method=object$method, cutoff=object$cutoff)
  class(out) <- "summary.mergenorm"
  out
}

print.summary.mergenorm <- function(x, ...){
  cat("* Merging Gaussian mixture components *\n\n")
  cat(" Method: ",x$method,", cutoff value: ",x$cutoff,"\n")
  cat(" Original number of components: ",x$onc,"\n")
  if (x$noise)
    cat("   (not including noise which is denoted by clustering=0)\n")
  cat(" Number of clusters after merging: ",x$mnc,"\n")
  if (x$method=="predictive"){
    cat(" Predictive strength values: \n")
    print(as.matrix(x$predvalues))
  }
  else{
    cat(" Values at which clusters were merged: \n")
    attr(x$valuemerged, "names") <- NULL
    print(cbind(((x$onc-1):(x$onc-length(x$valuemerged))),x$valuemerged))
  }
  cat(" Components assigned to clusters: \n")
  if (x$noise)
    print(as.matrix(c(0,x$mergedtonumbers)))
  else
    print(as.matrix(x$mergedtonumbers))
}
    
  
  
  


# legendposition="auto" or c(x,y)
# xdistances can be "mahalanobis", "euclidean" or "none"
# clustercol can be NULL, then allcol is used
weightplots <- function(z, clusternumbers="all", clustercol=2,
                        allcol=grey(0.2+((1:ncol(z))-1)*
                          0.6/(ncol(z)-1)),
                        lty=rep(1,ncol(z)),clusterlwd=3,
                        legendposition="none",
#                        xdistances="none",xdata=NULL,
                        weightcutoff=0.01,ask=TRUE, ...){
  k <- ncol(z)
  allnumbers <- 1:k
#  greylevels <- min(greyrange)+((1:k)-1)*(max(greyrange)-min(greyrange))/(k-1)
#  print(greylevels)
  n <- nrow(z)
#  if (xdistances=="mahalanobis"){
#    dcov <- cov(xdata)
#    dmatrix <- matrix(0,ncol=n,nrow=n)
#    for (i in 1:(n-1))
#      dmatrix[i,(i+1):n] <- dmatrix[(i+1):n,i] <-
#        mahalanobis(xdata[(i+1):n,],xdata[i,],dcov)
#  }
#  if (xdistances=="euclidean")
#    dmatrix <- as.matrix(dist(xdata))
  if (identical(clusternumbers,"all"))    
    clusternumbers <- allnumbers
  if (identical(legendposition,"auto")) legendposition=c(0,0.75)
  if (ask)
    par(ask=TRUE)
  iclustercol <- clustercol
  for (i in clusternumbers){
    oz <- order(-z[,i])
    sz <- z[oz,]
    szi <- sz[,i]>=weightcutoff
    szw <- sz[szi,]
#    if (xdistances!="none"){
#      xvals <- numeric(0)
#      xvals[1] <- 0
#      for (j in 2:n)
#        xvals[j] <- xvals[j-1]+dmatrix[oz[j],oz[j-1]]
#      print(xvals)
#      print(oz)
#      print(dmatrix[1:3,])
#      pxlab <-
#        paste("Pairwise acc. distances between observations ordered by
#               weight for cluster",i)     
#    }
#    else{
    xvals <- 1:n
    pxlab=paste("Observations ordered by posterior probability for component",
        i)
#    }
    if (is.null(clustercol))
      iclustercol <- allcol[i]
    plot(xvals[szi],szw[,i],
         col=iclustercol,ylim=c(0,1),type="l",lty=lty[i],
        xlab=pxlab,
         ylab="Estimated posterior probabilities")
    for (j in allnumbers[-i])
      points(xvals[szi],szw[,j],col=allcol[j],lty=lty[j],type="l")
    points(xvals[szi],szw[,i],
           col=iclustercol,lty=lty[i],lwd=clusterlwd,type="l")
    if (!identical(legendposition,"none")){
      leg.lwd=rep(1,ncol(z))
      leg.lwd[1] <- clusterlwd
      leg.txt <- c()
      leg.txt[1] <- paste("Cluster",i)
      k <- 2
      for (j in allnumbers[-i]){
        leg.txt[k] <- paste("Cluster",j)
        k <- k+1
      }
      leg.col <- c(iclustercol,allcol[allnumbers[-i]])
      leg.lty <- c(lty[i], lty[allnumbers[-i]])      
      legend(legendposition[1],legendposition[2],leg.txt,col=leg.col,
             lty=leg.lty,lwd=leg.lwd, ...)
    }
  }
  if (ask)
    par(ask=FALSE)
  invisible(z)
}
    
