
Inde <- function(H0=H0,pheno=pheno,envirX=envirX,clusRes=clusRes,S=S){
    cat('The individuals are independent, LM is used.\n',sep='')
    ## changed
    if(!H0){
        res <- lm(pheno~envirX+clusRes-1)
        len <- length(res$coefficients)
        alpha <- res$coefficients[1:(len-1)]
        bet <- res$coefficients[len]
    }else{
        res <- lm(pheno~envirX-1)
        alpha <- res$coefficients
        bet <- 0
    }
    zalpha <- envirX%*%alpha
    nbet <- clusRes*bet
    e <- pheno-zalpha-nbet
    sig2 <- as.numeric(var(e))
    return(list(alpha=alpha,bet=bet,sig2=sig2,e=e))
}


CNVtypeAnay <- function(pheno,pX,envirX,phi,S,FM,N,threshold, bet, alpha, sig2, sig2g,H0,thresAI,thresEM,itermax,signal,varSelection){

  seed<- sample(c(1:1000000),1)
  set.seed(seed)
  p <- rep(1/N,N)
  Ind.num <- dim(pX)[1]
  clusRes <- rep(NA,Ind.num)
  del <- list()
  Invdel <- list()
  mu <- list()
  LDdel <- rep(NA,N)
  Ddel <- rep(NA,N)
  mu1 <- rep(NA,N)
  q <- dim(pX)[2]

  ires <-  ClusProc(signal=signal,N=N,varSelection=varSelection,adjust=TRUE,thresMAF=0,scale=FALSE,thresSil=0)

  part1 <- data.frame(X=pX)
  part2 <- data.frame(Kn=ires$silWidth$adjusted$silRes.adjust)
  temp <- part2[match(rownames(part1),rownames(part2)),]
  KX <- cbind(X=part1,Kn=temp[,1])

  ## initial value for the means and covariance
  for(i in 1:N){
    mu[[i]] <- apply(as.matrix(KX[KX$Kn==(i-1),c(1:q)]),2,mean)
    del[[i]] <- cov(as.matrix(KX[KX$Kn==(i-1),c(1:q)]))
    Ddel[i] <- det(del[[i]])
    LDdel[i] <- log(Ddel[i])
    Invdel[[i]] <- solve(del[[i]])
  }
  
  ## variables preparing
  iter <- 0
  logLold <- Inf
  pin <- matrix(NA,S,N)
  P <- table(KX$Kn)/S
  clusFam <- rep(NA,FM)
  phiInv <- solve(phi)
  
  while(1){
    iter=iter+1
    cat('Iteration ',iter,':\n',sep='')
    for(i in 1:S)
      for(n in 1:N)
        pin[i,n]  <- Ddel[n]^(-1/2)*exp(-1/2*(pX[i,]-mu[[n]])%*%Invdel[[n]]%*%as.matrix(pX[i,]-mu[[n]]))*P[n]
    ## set the method to be used
    for(i in 1:S)
      clusRes[i] <- which(pin[i,]==max(pin[i,]))-1
    for(i in 1:N)
      P[i] <- mean(clusRes==(i-1))
    Sdata <- data.frame(pX=pX,clusRes=clusRes)
    ## updata the signal model
    for(i in 1:N){
      mu[[i]] <- apply(as.matrix(Sdata[Sdata$clusRes==(i-1),c(1:q)]),2,mean)
      del[[i]] <- cov(as.matrix(Sdata[Sdata$clusRes==(i-1),c(1:q)]))
      Ddel[i] <- det(del[[i]])
      LDdel[i] <- log(Ddel[i])
      Invdel[[i]] <- solve(del[[i]])
    }
    ## updata the phenotype model
    pheno[which(is.na(pheno))] <- 0

    if(all(phi[lower.tri(phi)]==0,phi[upper.tri(phi)]==0)){
        temp <- Inde(H0=H0,pheno=pheno,envirX=envirX,clusRes=clusRes,S=S)
        sig2g <- 0
        bet <- temp$bet
        sig2 <- temp$sig2
        alpha <- temp$alpha
        e <- temp$e
    }else{
        res <- doPolyGenic(envirX=envirX,snp=clusRes,pheno=pheno,phi=phi,phiInv=phiInv,thresEM=thresEM,thresAI=thresAI,H0=H0)

        if(length(res)==1) {
            temp <- Inde(H0=H0,pheno=pheno,envirX=envirX,clusRes=clusRes,S=S)
            sig2g <- 0
            bet <- temp$bet
            sig2 <- temp$sig2
            alpha <- temp$alpha
            e <- temp$e
        }else {
            cat('The individuals are correlated, LMM is used.\n',sep='')
            sig2 <- res$sig2
            sig2g <- res$sig2g
            len <- length(res$para)
            if(!H0){
                alpha <- res$para[1:(len-1)]
                bet <- res$para[len]}else{
                    alpha <- res$para[1:len]
                    bet <- 0
                }
            zalpha <- envirX%*%alpha
            nbet <- clusRes*bet
            e <- pheno-zalpha-nbet
        }        
    }

    V <- sig2*diag(S)+sig2g*phi
    ## loglikelihood calculation
    logL_sig <- 1
    for(i in 1:S)
        logL_sig <- logL_sig-1/2*log(Ddel[clusRes[i]+1])-1/2*(pX[i,]-mu[[clusRes[i]+1]])%*%solve(del[[clusRes[i]+1]])%*%as.matrix(pX[i,]-mu[[clusRes[i]+1]])
    logL <-logL_sig-1/2*determinant(V)$modulus-1/2*t(e)%*%solve(V)%*%e+sum(log(P[clusRes]))
    ##if(abs(logL-logLold)<threshold) ttt <- 1 else ttt <- 0
    #print(ttt)
    if(abs(logL-logLold)<threshold) break
    if(iter>itermax) break
    logLold <- logL
    ## print(paste("The logliklihood without of costants is now",logL,"."))
}

  return(list(clusRes=clusRes,bet=bet,alpha=alpha,sig2=sig2,sig2g=sig2g,logL1=round(logL,6)))
  
}




##' This function tests the association of CNV with continuous trait of interest. Two statistics are provided for different strategies with the intensity measurement.
##'
##' @title CNV association test procedure
##' @param signal The matrix of intensity measurements. The row names must be consistent with the Individual ID in fam file.
##' @param fam The FAM file which follows the format defined in PLINK. 
##' @param envirX The matrix of environmental variables. The intercept is automatically included and it does not need to be in this matrix.
##' @param phi The correlation matrix between individuals. It can be built with the kinship coefficient or the estimated correlation matrix with SNP data. Free software that builds this matrix is available, and one of them can be downloaded at \url{http://biostat.ac.kr/fqls/} The default is an identity matrix and it is for independent samples.
##' @param N Number of clusters one wants to fit to the data. N needs to be larger than 1 and if it is 1, error will be returned. It can be estimated with the function \code{\link{ClusProc}}.
##' @param varSelection Factor. For specifying how to handle the intensity values. It must take value on 'RAW', 'PC.9', 'PC1'and 'MEAN'. If the value is 'RAW', then the raw intensity value will be used. If it is 'PC.9', then the first several PCA scores which account for 90\% of all the variance will be used. If the value is 'PC1', then the first PCA scores will be used. If the value is 'MEAN', the mean of all the probes will be used. The default method is 'PC1'.
##' @param H0 Logicals. If it is TRUE (the default), all parameters are estimated under the assumption that there is no genetic association between CNV and phenotypes. If it is FALSE, parameters are estimated under the null or alternative hypothesis.
##' @param threshold Optional number of convergence threshold. The iteration stops if the absolute difference of log likelihood between successive iterations is less than it. The default threshold 1e-05 will be used if it's missing.
##' @param itermax Optional. The iteration stops if the times of iteration is large than this value. The default number 8 will be used if it's missing.
##' @param thresEM Optional number of convergence threshold in the EM (expectation-maximization method) procedure. The default threshold 0.005 will be used if it's missing.
##' @param thresAI Optional number of convergence threshold in the AI (average information method) procedure. The default threshold 1e-05 will be used if it's missing.
##' @return It returns object of class 'asso'. The result is obtained under the null hypothesis if H0 is TRUE, otherwise the result is obtained under null or alternative hypothesis.
##' \item{para}{The parameter estimations for the best fit.}
##' \item{clusRes}{The clustering assignment for each individual.}
##' @author Meiling Liu, Sungho Won and Weicheng Zhu
##' @examples
##' # Fit the data under the assumption that there are 3 clusters
##' fit.pc <- AssoTestProc(signal=signal,fam=fam,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
##' @export

AssoTestProc <- function(signal,fam,envirX,phi,N,varSelection=c('PC1','RAW','PC.9','MEAN'),H0=TRUE,threshold=1e-05,itermax=8,thresEM=0.005,thresAI=1e-05){  

    varSelection <- match.arg(varSelection)
    rn_signal <- row.names(signal)
    rn_envirX <- row.names(envirX)
    iid <- fam[,2]
    if(sum(!(rn_signal%in%iid))>0) {
        cat('The row name of signal data is not coincident with individual ID.\n')
        return(0)
    }
    if(sum(!(rn_envirX%in%iid))>0) {
        cat('The row name of covariant data is not coincident with individual ID.\n')
        return(0)
    }

    envirX <- cbind(1,envirX)
    len <- sum(!duplicated(fam[,1]))
    pheno <- as.numeric(fam[,6])
    S <- length(pheno)
    if(is.null(phi)) phi <- diag(S)
    
    p <- rep(1/N,N)
    INsnp <- rbinom(S,2,p)
    X  <- cbind(envirX,INsnp)
    y  <- matrix(pheno,ncol=1)
    yy <- t(y)%*%y
    Xy <- t(X)%*%y
    XX <- t(X)%*%X
    phiInv  <- solve(phi)
    trphiInv<- sum(diag(phiInv))
    q  <- ncol(X)
    residual <- y - X%*%matrix(solve(XX)%*%Xy,ncol=1)
    
    inis     <- getInitial(residual,phi,S)

    in_sig2g <- inis[2]
    in_sig2 <- inis[1]
    pca <- princomp(signal)
    pX1 <- signal%*%loadings(pca)[,1]
    cut <- 1

    if(varSelection=='PC.9'){
        prop <- 0.9
        pca <- princomp(signal) ## PCA
        sds <- pca$sdev
        vars <- sds^2
        varprop  <- vars/sum(vars)
        cumvars <- as.vector(cumsum(varprop))
        while(cumvars[cut]<prop)  cut=cut+1
        ## print(paste0("we use Comp.1 to Comp.",cut))
        Invcov <- matrix(0,nrow=cut,ncol=cut)
        diag(Invcov) <- 1/vars[1:cut]
        comptable <- data.frame(sdev=pca$sdev,vars=vars,cumu=cumvars)
        ## print("Variable Selection")
        ## print(t(comptable[1:(cut+1),]))
        coef <- pca$loadings[,1:cut]
        pX <- signal%*%coef
  }
  if(varSelection=='RAW') {
    pX <- signal
    cut <- ncol(pX)
  }
  if(varSelection=='MEAN')    pX <- as.matrix(apply(signal,1,mean))
  if(varSelection=='PC1'){
    pca <- princomp(signal)
    pX <- signal%*%loadings(pca)[,1]
  }

    ## under all
    if(!H0){
        fit <- lm(pheno~pX1+envirX-1)
        bet <- fit$coefficient[1]
        alpha <- fit$coefficient[-1]
        res_H1<- CNVtypeAnay(pX=pX,pheno=pheno,envirX=envirX,phi=phi,S=S,FM=len,N=N,threshold=threshold,bet=bet,alpha=alpha,sig2=in_sig2,sig2g=in_sig2g,H0=FALSE,thresEM=thresEM,thresAI=thresAI,itermax=itermax,signal=signal,varSelection=varSelection)
     } else {
         fit <- lm(pheno~envirX-1)
         bet <- 0
         alpha <- fit$coefficient
         res_H0 <- CNVtypeAnay(pX=pX,pheno=pheno,envirX=envirX,phi=phi,S=S,FM=len,N=N,threshold=threshold,bet=bet,alpha=alpha,sig2=in_sig2,sig2g=in_sig2g,H0=TRUE,thresEM=thresEM,thresAI=thresAI,itermax=itermax,signal=signal,varSelection=varSelection)
}

    if(H0){
        resfinal <- list(para=list(bet=res_H0$bet,alpha=res_H0$alpha,sig2=res_H0$sig2,sig2g=res_H0$sig2g),clusRes=res_H0$clusRes,H0=H0)
    } else
        {
      resfinal <- list(para=list(bet=res_H1$bet,alpha=res_H1$alpha,sig2=res_H1$sig2,sig2g=res_H1$sig2g),clusRes=res_H1$clusRes,H0=H0)
  }

    class(resfinal) <- 'asso'
    return(resfinal)
}







##' Prints formatted results from the association study returned by \code{\link{AssoTestProc}}.
##'
##' @title Prints association study results
##' @param x The association study results obtained from the \code{\link{AssoTestProc}}.
##' @param ... Usual arguments passed to the print function.
##' @author Meiling Liu
##' @method print asso
##' @examples
##' # Fit the data under the assumption that there are 3 clusters
##' asso.fit <- AssoTestProc(signal=signal,fam=fam,envirX=envirX,phi=phi,N=3,varSelection='PC.9')
##' print(asso.fit)
##' @export
print.asso <- function(x, ...){
    if(x$H0) cat('\n\n Under H0:\n') else 
    cat('\n\n Under H0 and H1:\n')
    cat('The coefficent:\n')
    temp <- c(bet=0,alpha=round(x$para$alpha,4))
    print(temp,quote=FALSE,row.names=FALSE)
    cat('The estimated variance:\n')
    temp <- data.frame(sig2=x$para$sig2,sig2g=x$para$sig2g)
    print(temp,quote=FALSE,row.names=FALSE)
}

