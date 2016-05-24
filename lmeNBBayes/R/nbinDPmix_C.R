## ===== Subroutines for multMeta =========
cov2corNA <- function(x)
{
  cv2cr <- x
  cv2cr[!is.na(x)] <- cov2cor(matrix(x[!is.na(x)],
                                     nrow=sum(!is.na(x)[1,]),
                                     ncol=sum(!is.na(x)[,1])))
  ##diag(cv2cr) <- 1
  return(cv2cr)
}



intfun <-  function(D,meanlog,sdlog,N)  D*log(1+N/D)*dlnorm(D,meanlog,sdlog)
E.KN <- function(N,meanlog,sdlog,width=3)
  {
    ## compute E(K|N) v.s. mD under D ~ dlnorm(meanlog,sdlog)
    l <- integrate(f=intfun,lower=0,
                   upper=exp(meanlog - width * sdlog),
                   meanlog=meanlog,sdlog=sdlog,N=N)$value
    m <- integrate(f=intfun,
                   lower=exp(meanlog - width * sdlog),
                   upper=exp(meanlog + width * sdlog),
                   meanlog=meanlog,sdlog=sdlog,N=N)$value
    u <- integrate(f=intfun,
                   lower=exp(meanlog + width * sdlog), upper=Inf,
                   meanlog=meanlog,sdlog=sdlog,N=N)$value
    return(l + m + u)
  }



cor2cov <- function(cor.mat,sd.vec)
{
  cov.mat <- cor.mat
  temp <- cov.mat*sd.vec ## each row of cov.mat is multiplied by each entry of sd.vec
  temp <- t(t(temp)*sd.vec) 
  return(temp)
}


adjustPosDef <- function (mat,zero=1e-15) 
{
    eg <- eigen(mat)
    if (sum(eg$values < zero) > 0) {
        note <- "The non positive semi-definite-ness is adjusted by the truncation"
        adjust.mat <- 0
        for (i in 1:length(eg$values)) adjust.mat <- adjust.mat + 
            max(c(eg$values[i], zero)) * outer(eg$vector[, i], 
                eg$vector[, i])
        if (!is.null(colnames(mat))) {
            colnames(adjust.mat) <- colnames(mat)
            rownames(adjust.mat) <- rownames(mat)
        }
    }
    else {
        adjust.mat <- mat
        note <- "The matrix is already positive semi-definite"
    }
    return(list(adjust.mat = adjust.mat, note = note))
}


logL.G <- function(Y,ID,
                 X,vec.beta,vec.gPre
                 )
  {
    ## focused likelihood Pr(Yij | beta, G)
    ## gPre, gNre is a vector of length N
    ## if (is.null(vec.gNew)) labelnp <- rep(0,length(Y))

    if (is.vector(X))
      X <- matrix(X,ncol=1)
    logL <- 0
    uniID <- unique(ID)
    for (ipat in 1 : length(uniID) )
      {
        patID <- uniID[ipat]
        iY <- Y[ID==patID]
        iX <- X[ID==patID,,drop=FALSE]
        i.gPre <- vec.gPre[ipat]
        logL <- logL + sum(dnbinom(iY,
                                   size=exp(iX%*%vec.beta),prob=i.gPre,
                                   log=TRUE))
      }
    return (logL)
  }

logL.aGh.rGh <- function(Y,ID,
                  X,vec.beta,
                  vec.aGs,vec.rGs ## LENGTH OF N
                  )
  {
    ## focused likelihood Pr(Yij | beta, aGs[h], rGs[h])
    ## gPre, gNre is a vector of length N
    ## if (is.null(vec.gNew)) labelnp <- rep(0,length(Y))
    
    if (is.vector(X))
      X <- matrix(X,ncol=1)
    logL <- 0
    uniID <- unique(ID)
    for (ipat in 1 : length(uniID) )
      {
        patID <- uniID[ipat]
        iY <- Y[ID==patID]
        iX <- X[ID==patID,,drop=FALSE]
        aGh <- vec.aGs[ipat]
        rGh <- vec.rGs[ipat]
        Lchoose <- log.choose(a=exp(iX%*%vec.beta),n=iY)
        logL <- logL + sum(Lchoose)+
          lbeta(sum(exp(iX%*%vec.beta))+aGh,sum(iY)+rGh)-lbeta(aGh,rGh)
      }
    return (logL)
  }


  
log.choose <- function(a,n)
  {
    ## Compute choose(a+n-1,n) for non-integer a
    lgamma(a+n)-lgamma(a)-lfactorial(n)
  }
llk.FG_i <- function(ys,rs,aGs,bGs,ps)
{
  ## Pr( {Yij}_j=1^ni | beta, F_G )
  ##     = prod_{j=1}^ni choose(rij+yij-1,yij)*sum_ih^M ps[ih]*beta(r_ip+ +aGs[ih], y_ip+ +bGs[ih])/beta(aGs[ih],bGs[ih])
  ## ys: length ni
  ## rs: length ni
  ## aGs: length M
  ## rGs: length M
  ## ps: length M
  term1 <- sum(log.choose(a=rs,n=ys))
  ratio <- exp(lbeta(sum(rs)+aGs,sum(ys)+bGs) - lbeta(aGs,bGs))
  term2 <- log(sum(ps*ratio))
  return(term1+term2)
}
llk.FG <- function(Ys,Xs,ID,mat.beta,mat.aGs,mat.bGs,mat.ps)
  {
    ## mat.aGs: B* by M where B* = length(useSample)
    ## mat.bGs: B* by M
    ## mat.ps:  B* by M
    uniID <- unique(ID)
    llkeach <- rep(NA,nrow(mat.aGs))
    for (iB in 1 : nrow(mat.aGs))
      {
        beta <- mat.beta[iB,]
        aGs <- mat.aGs[iB,]
        bGs <- mat.bGs[iB,]
        ps <- mat.ps[iB,]

        temp.llk <- 0
        for (ipat in 1 : length(uniID))
          {
            Yi <- Ys[ID==uniID[ipat]]
            Xi <- Xs[ID==uniID[ipat],,drop=FALSE]
            
            rs <- exp(Xi%*%beta) ## length ni
            temp.llk <- temp.llk + llk.FG_i(ys=Yi,rs=rs,aGs=aGs,bGs=bGs,ps=ps)
            if (!is.finite(temp.llk)) stop("!!")
          }

        llkeach[iB] <- temp.llk
      }
    return(llkeach)
  }

getDIC.FG <- function(olmeNBB, data, ID, useSample, lower.alpha=0.0001,upper.alpha=0.99999,inc.alpha=0.0005)
  {
    ## focused likelihood = P( Yi | beta, F_G )
    alphas <- seq(lower.alpha,upper.alpha,inc.alpha)
    ## focused likelihood integrates out the component labels
    ## P( Yij, j=1,..,ni | beta, RE dist)
    formula <- olmeNBB$para$formula
    X <- model.matrix(object = formula, data = data)
    Y <- model.response(model.frame(formula = formula, data = data))

    if (olmeNBB$para$DP == 0) stop("Use getDIC for parametric Bayesian procedure")
    if (is.vector(olmeNBB$beta)) 
        olmeNBB$beta <- matrix(olmeNBB$beta, ncol = 1)
    
    D.bar <- -2*mean(llk.FG(Ys=Y,Xs=X,ID=ID,
                            mat.beta = olmeNBB$beta[useSample,, drop = FALSE],
                            mat.aGs  = olmeNBB$aGs[useSample,, drop = FALSE],
                            mat.bGs = olmeNBB$rGs[useSample,, drop = FALSE],
                            mat.ps = olmeNBB$weightH1[useSample,,drop=FALSE]))
    ## compute the Pr( {Yij}_{j=1}^ni | bar.beta, bar.F_G)
    ## return a length(alphas) by length(B) matrix
    qM <- dqmix(
                weightH1 = olmeNBB$weightH1[useSample,],
                aGs = olmeNBB$aGs[useSample,],
                rGs = olmeNBB$rGs[useSample,],
                alphas = alphas, ## the grid of points at which the density is evaluated
                dens=TRUE
                )[[1]]
    ## bar.F_G is computed by averaging the estimated density at each grid of points
    aggregate.dens <- colMeans(qM)
    ## bar.beta is average of beta^{(b)} after thinning 
    aggregate.beta <- colMeans(olmeNBB$beta[useSample,,drop=FALSE])

    llk.pat <- 0
    uniID <- unique(ID)
    for (ipat in 1 : length(uniID))
      {
        ys <- Y[ID==uniID[ipat]]
        xs <- X[ID==uniID[ipat],,drop=FALSE]
        rs <- exp(xs%*%aggregate.beta)
        ## Pr( {Yij}_j=1^ni | bar.beta, bar.F_G) = int_0^1 prod_{j=1}^ni Pr( Yij | bar.beta, g) bar.f_G(g) dg
        vec.dens <- rep(NA,length(alphas))
        for (ialpha in 1 : length(alphas)){
          vec.dens[ialpha]<- prod(dnbinom(ys,size=rs,prob=alphas[ialpha]))*aggregate.dens[ialpha]
        }
        llk.pat <- llk.pat + log(sum(vec.dens*inc.alpha))
      }
    D.hat <- -2 *llk.pat
    effect.para <- D.bar - D.hat
    DIC <- effect.para + D.bar
    return(list(DIC = DIC, effect.para = effect.para))
  }


getDIC_aGh.rGh.G <- function(olmeNBB,data,ID,useSample)
{
  formula <- olmeNBB$para$formula
  X <- model.matrix(object=formula,data=data)
  Y <- model.response(model.frame(formula=formula,data=data))
  ## The effective number of parameters of the model is computed
  ## Dbar: posterior mean of the deviance bar(-2*logLlik)
  D.bar <- -2*mean(olmeNBB$logL[useSample]) ## + constant
  ## Dhat: -2*logLik(theta.bar)
  if (is.vector(olmeNBB$beta)) olmeNBB$beta <- matrix(olmeNBB$beta,ncol=1)
  ##if (is.vector(olmeNBB$gNew)) olmeNBB$gNew <- matrix(olmeNBB$gNew,ncol=1)
  ##if (olmeNBB$para$Reduce){
    if (olmeNBB$para$DP==0)
      {
        ## parametric procedure:
        olmeNBB$aGs_pat <- matrix(olmeNBB$aG,nrow=length(olmeNBB$aG),ncol=length(unique(ID)))
        olmeNBB$rGs_pat <- matrix(olmeNBB$rG,nrow=length(olmeNBB$rG),ncol=length(unique(ID)))
      }
    ## Pr( { Yij }_j=1^ni | beta, aGh, rGh)
    D.hat <- -2*logL.aGh.rGh(Y=Y,ID=ID,X=X,
                             vec.beta = colMeans(olmeNBB$beta[useSample,,drop=FALSE]), 
                             vec.aGs = colMeans(olmeNBB$aGs_pat[useSample,,drop=FALSE]),
                             vec.rGs = colMeans(olmeNBB$rGs_pat[useSample,,drop=FALSE])
                             )
  ## }else{
  ##   ## Pr( { Yij }_j=1^ni | beta, G)
  ##   D.hat <- -2*logL.G(Y=Y,ID=ID,X=X,
  ##                    vec.beta = colMeans(olmeNBB$beta[useSample,,drop=FALSE]),
  ##                    vec.gPre = colMeans(olmeNBB$g1s[useSample,,drop=FALSE])
  ##                    )
  ## }
  effect.para <- D.bar - D.hat
  DIC <- effect.para + D.bar
  return (list(DIC=DIC,effect.para=effect.para))
}


getDIC <- function(olmeNBB, data,
                  ID, useSample=NULL,focus = c("FG","G","aGh.rGh","para"), 
                  lower.alpha=0.0001,upper.alpha=0.99999,inc.alpha=0.0005){
  if (length(focus)==1) focus <- focus[1]
  if (is.null(useSample)) useSample <- 1 : olmeNBB$para$B
  ## Pr(Yij | beta, aGs[h], rGs[h], pis[ih], ih=1,...,M ) works only for semiparametric procedure
  if (focus== "FG")
    re <- getDIC.FG(olmeNBB=olmeNBB, data=data,
                    ID=ID, useSample=useSample,
                    lower.alpha=lower.alpha,upper.alpha=upper.alpha,inc.alpha=inc.alpha)
  ## Pr(Yij | beta, aGs[h], rGs[h], pis[ih])
  ## Pr(Yij | beta, G)
  if (focus %in% c("G","aGh.rGh","para"))
    re <- getDIC_aGh.rGh.G(olmeNBB=olmeNBB,data=data,
                           ID=ID,useSample=useSample)
  return(re)
}



tempCdens <- function(mu,Sigma,Random){
  eS <- eigen(Sigma)
  ev <- eS$values;evec <- eS$vector
  return( mu + evec %*% diag(ev) %*% (Random))
}
index.batch.Bayes <- function(data,labelnp,ID,olmeNBB,thin=NULL,printFreq=10^5,unExpIncrease=TRUE)
  {
    burnin <- olmeNBB$para$burnin
    if (is.null(thin))
      {
        if (is.null(olmeNBB$para$thin))  thin <- 1
        else  thin <- olmeNBB$para$thin
      }
    formula <- olmeNBB$para$formula
    X <- model.matrix(object=formula,data=data)
    Y <- model.response(model.frame(formula=formula,data=data))

    if (is.vector(X)) X <- matrix(X,ncol=1)
    NtotAll <- length(Y)
    if (nrow(X)!= NtotAll) stop ("nrow(X) != length(Y)")
    if (length(ID)!= NtotAll)  stop ("length(ID)!= length(Y)")
    if (length(labelnp)!= NtotAll)  stop ("labelnp!= length(Y)")
  
    if (olmeNBB$para$DP==0)
      {
        olmeNBB$aGs <- matrix(olmeNBB$aG,ncol=1)
        olmeNBB$rGs <- matrix(olmeNBB$rG,ncol=1)
        olmeNBB$weightH1 <- matrix(1,ncol=1,nrow=length(olmeNBB$aGs))
      }
    
    useSample <- useSamp(thin=thin, burnin=burnin, B=nrow(olmeNBB$beta))
    maxni <- max(tapply(rep(1,length(ID)),ID,sum))


    ## change the index of ID to numeric from 1 to # patients
    temID <- ID  
    Npat <- length(unique(temID))
    uniID <- unique(temID)
    ID <- rep(NA,length(temID))
    for (i in 1 : length(uniID))
      {
        ID[temID == uniID[i]] <- i
      }
    ## ID starts from zero
    mID <- ID-1

    
    mat_betas <- olmeNBB$beta[useSample,,drop=FALSE]
    mat_aGs <- olmeNBB$aGs[useSample,,drop=FALSE]
    mat_rGs <- olmeNBB$rGs[useSample,,drop=FALSE]
    mat_weightH1 <- olmeNBB$weightH1[useSample,,drop=FALSE]
    B <- nrow(mat_betas) 
    M <- ncol(olmeNBB$aGs)
    
    if (unExpIncrease){
      
      re <- .Call("index_b_Bayes_UnexpIncrease",
                  as.numeric(Y), as.numeric(c(X)),
                  as.numeric(labelnp), as.integer(mID),
                  as.integer(B),as.integer(maxni),as.integer(M),
                  as.integer(Npat),
                  as.numeric(c(t(mat_betas))), 
                  as.numeric(c(t(mat_aGs))),
                  as.numeric(c(t(mat_rGs))),
                  as.numeric(c(t(mat_weightH1))),
                  as.integer(printFreq),
                  package ="lmeNBBayes")
    }else{
      re <- .Call("index_b_Bayes_UnexpDecrease",
                  as.numeric(Y), as.numeric(c(X)),
                  as.numeric(labelnp), as.integer(mID),
                  as.integer(B),as.integer(maxni),as.integer(M),
                  as.integer(Npat),
                  as.numeric(c(t(mat_betas))), 
                  as.numeric(c(t(mat_aGs))),
                  as.numeric(c(t(mat_rGs))),
                  as.numeric(c(t(mat_weightH1))),
                  as.integer(printFreq),
                  package ="")
    }
    
    re <- matrix(re[[1]],nrow=B,ncol=Npat,byrow=TRUE)
    colnames(re) <- unique(temID)
    rownames(re) <- paste("sim",1:B,sep="")
    ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
    for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0
        
    patwoNorO <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
    if (length(patwoNorO)==0) patwoNorO <- NULL;

    re[,patwoNorO] <- NA

    res <- list()
    res$condProb <- re
    res$condProbSummary <- condProbCI(ID=temID,condProb=re)
    
    res$para$labelnp <- labelnp
    res$para$ID <-ID
    res$para$CEL <- Y
    res$para$thin <- thin
    res$para$B <- B
    res$para$burnin <- burnin
    class(res) <- "IndexBatch"
    return(res)
  }


condProbCI <- function(ID,condProb)
  {
    uniID <- unique(ID)
    
    condProb <- cbind(
                      colMeans(condProb),
                      apply(condProb,2,sd),
                      t(apply(condProb,2,
                              quantile,na.rm=TRUE,
                              prob=c(0.025,0.975)))
                      )
    colnames(condProb) <- c("CondProb","SE","2.5%","97.5%")
    rownames(condProb) <- uniID
    return(condProb)
  }



newCat <- function(label1,label2)
  {
    ## The length of label1 and label2 are the same
    newCat <- rep(0,length(label1))
    ## label1 and label2 must be at the same length
    for (ilabel1 in unique(label1))
      {
        for (ilabel2 in unique(label2))
          {
            newCat[label1==ilabel1 & label2 == ilabel2] <- paste(ilabel1,ilabel2,collapse=":",sep=":")
          }
      }
    return (as.factor(newCat))
  }


repeatAsID <-function(values,ID)
  {
    ## ID does not have to be numerics
    ## values is a vector of length length(unique(ID))
    ivalue <- 1
    output <- rep(NA,length(ID))
    for (iID in unique(ID))
      {
        output[ID==iID] <- values[ivalue]
        ivalue <- ivalue + 1
      }
    return (output)
  }


slim <- function(vec,ID)
  {

    uniID <- unique(ID)
    re <- rep(NA,length(uniID))
    for (i in 1 : length(uniID))
      {
        re[i] <- vec[ID==uniID[i]][1]
      }
    
    return(re)
  }




lnpara <- function(EX,SDX){
  ## log(D) ~ normal(mean=mulog,sd=sdlog) 
  ## D ~ log-normal(meanlog=meanlog,sdlog=sdlog)
  ## Given the desired E(D) and SD(D), find meanlog = E(log(D)) and sdlog = SD(log(D)) parameters of log-normal

  ## E(logX) = exp(EX + SDX^2/2)
  ## Var(logX) = ( exp(SDX^2) - 1 )*exp(2*EX + SDX^2 ) = ( exp(SDX^2) - 1 )*( E(logX)^2 )
  temp <- (SDX^2)/(EX^2) 
  meanlog <- log(EX)-log(temp+1)/2
  sdlog <- sqrt(log(1+temp))
  return(list(meanlog=meanlog,sdlog=sdlog))
}

invlnpara <- function(meanlog,sdlog){
  ## log(D) ~ normal(mean=mulog,sd=sdlog) 
  ## D ~ log-normal(meanlog=meanlog,sdlog=sdlog)
  ## Given the desired E(D) and SD(D), find meanlog = E(log(D)) and sdlog = SD(log(D)) parameters of log-normal
  mean <- exp(meanlog+(sdlog^2)/2)
  return(list(mean=mean,
              sdlog=(exp(sdlog^2)-1)*(mean^2)))
}



lnpara <- function(EX,SDX){
  ## log(D) ~ normal(mean=mulog,sd=sdlog) 
  ## D ~ log-normal(meanlog=meanlog,sdlog=sdlog)
  ## Given the desired E(D) and SD(D), find meanlog = E(log(D)) and sdlog = SD(log(D)) parameters of log-normal

  temp <- (SDX^2)/(EX^2) 
  meanlog <- log(EX)-log(temp+1)/2
  sdlog <- sqrt(log(1+temp))
  return(list(meanlog=meanlog,sdlog=sdlog))
}

int.M <- function(D,M,mean.norm,sd.norm){
  if (M == 1)
    1/(D+1)*dlnorm(D,meanlog=mean.norm,sdlog=sd.norm)
  else
    (1-1/(D+1))^(M-1)*dlnorm(D,meanlog=mean.norm,sdlog=sd.norm)
}



piM <- function(mean.norm=0, sd.norm=1, width=2,M) {
  integral.l <- integrate(f=int.M, lower=0, upper=exp(mean.norm - width * sd.norm),
                          mean.norm=mean.norm, sd.norm=sd.norm,M=M)$value
  integral.m <- integrate(f=int.M, lower=exp(mean.norm - width * sd.norm),
                          upper=exp(mean.norm + width * sd.norm),
                          mean.norm=mean.norm, sd.norm=sd.norm,M=M)$value
  integral.u <- integrate(f=int.M, lower=exp(mean.norm + width * sd.norm), upper=Inf,
                          mean.norm=mean.norm, sd.norm=sd.norm,M=M)$value
  return(integral.l + integral.m + integral.u)
}





lmeNBBayes <- function(formula,          ##   A vector of length sum ni, containing responses
                       data,
                       ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
                       ID,         ##   A Vector of length sum ni, indicating patients
                       B = 105000, ##     A scalar, the number of Gibbs iteration 
                       burnin = 5000,  
                       printFreq = B,
                       M = NULL,
                       probIndex = FALSE,
                       thin =1, ## optional
                       labelnp=NULL, ## necessary if probIndex ==1
                       epsilonM = 1e-4,## nonpara
                       para = list(mu_beta = NULL,Sigma_beta = NULL,max_aG=30,mu_lnD=NULL,sd_lnD=NULL),
                       DP=TRUE,
                       thinned.sample=FALSE,
                       proposalSD = NULL
                       ## Does not matter if DP=FALSE
                       )
  {
    Xmodmat <- model.matrix(object=formula,data=data)
    covariatesNames<- colnames(Xmodmat)
    Y <- model.response(model.frame(formula=formula,data=data))
    ## If ID is a character vector of length sum ni,
    ## it is modified to an integer vector, indicating the first appearing patient
    ## as 1, the second one as 2, and so on..

    ## This code generate samples from NBRE model with constant random effect ~ DP mixture of Beta
    if (is.vector(Xmodmat)) Xmodmat <- matrix(Xmodmat,ncol=1)
    NtotAll <- length(Y)
    if (nrow(Xmodmat)!= NtotAll) stop ("nrow(Xmodmat) != length(Y)")
    if (length(ID)!= NtotAll)  stop ("length(ID)!= length(Y)")
    if (!is.null(labelnp) & length(labelnp)!= NtotAll)  stop ("labelnp!= length(Y)")
    if (thinned.sample & (!is.numeric(thin)) & (!is.numeric(burnin)))
      stop("If you only want thinned samples, you must give the thinning and burnin parameters")
    
    dims <- dim(Xmodmat)
    Ntot <- dims[1]
    pCov <- dims[2]
    ## proposalSD
    if (is.null(proposalSD)) proposalSD <- list() 
    if (is.null(proposalSD$min$aG))   proposalSD$min$aG <- 0.05
    if (is.null(proposalSD$min$rG))   proposalSD$min$rG <- 0.05
    if (is.null(proposalSD$min$beta)) proposalSD$min$beta <- 1e-4
    if (is.null(proposalSD$max$aG))   proposalSD$max$aG <- 5
    if (is.null(proposalSD$max$rG))   proposalSD$max$rG <- 5
    if (is.null(proposalSD$max$beta)) proposalSD$max$beta <- 10
    ## prior para
    if (is.null(para$mu_beta)) para$mu_beta <- rep(0,pCov)
    if (is.null(para$Sigma_beta))para$Sigma_beta <- diag(10,pCov)
    if (is.null(para$max_aG)) para$max_aG <- 30
    ## Checker
    if (length(para$mu_beta)!=ncol(Xmodmat))
      stop("The dimension of the fixed effect hyperparameter is wrong!")
    if (!is.matrix(para$Sigma_beta) & length(para$Sigma_beta)==1) para$Sigma_beta <- matrix(para$Sigma_beta,1,1)
    if (nrow(para$Sigma_beta) != ncol(Xmodmat) || ncol(para$Sigma_beta) != ncol(Xmodmat) ) stop()
    
    if (DP)
      {
        ## Additional arguments:
        ## para$mu_lnD, para$sd_lnD, their corresponding proposal variance's upper and lower bounds and M
        EX <- 0.5
        SDX <- 0.5
        logDpara <- lnpara(EX=EX,SDX=SDX)
        if (is.null(para$mu_lnD)){
          para$mu_lnD <- logDpara$meanlog
        }
        if (is.null(para$sd_lnD)){
          para$sd_lnD <- logDpara$sdlog
        }
        if (length(para$mu_lnD) != 1 ||length(para$sd_lnD) != 1) stop()
        atlnD <- length(para$mu_beta)
        if (is.null(proposalSD$min$lnD)) proposalSD$min$lnD <- 1e-6
        if (is.null(proposalSD$max$lnD)) proposalSD$max$lnD <- 1
            
        
        if (is.null(M)){
          for (M in 1 : 1000)
            {
              EpiM <- piM(M=M,
                          mean.norm=para$mu_lnD,
                          sd.norm=para$sd_lnD)
              
            if (EpiM < epsilonM ) break
            }
        }
      }
    
    eigenTemp <- eigen(para$Sigma_beta)
    evalue_sigma_beta <- eigenTemp$values
    evec_sigma_beta <- c(eigenTemp$vector)
    if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
    Inv_sigma_beta <- c( solve(para$Sigma_beta) )

    X <- c(Xmodmat) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }

    ## change the index of ID to numeric from 1 to # patients
    temID <- ID  
    N <- length(unique(temID))
    uniID <- unique(temID)
    ID <- rep(NA,length(temID))
    for (i in 1 : length(uniID))
      {
        ID[temID == uniID[i]] <- i
      }
    
    mID <- ID-1
    ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    
    maxni <- max(tapply(rep(1,length(ID)),ID,sum))
    Npat <- length(unique(ID))
    
    if (probIndex)
      {
        ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
        patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
        for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0
        
        patwoNorO <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
        if (length(patwoNorO)==0) patwoNorO <- -1000;
      }else{
        
        labelnp <- rep(0,length(Y))
      }
    
    Btilde <- B
    if (B %% thin != 0 ) stop("B %% thin !=0")
    if (burnin %% thin !=0) stop("burnin %% thin !=0")
    min_proposalSD <- c(aG=proposalSD$min$aG,rG=proposalSD$min$rG,beta=proposalSD$min$beta);
    max_proposalSD <- c(aG=proposalSD$max$aG,rG=proposalSD$max$rG,beta=proposalSD$max$beta);
    if (DP)
      {
        cat("\n M:",M)
        min_proposalSD <- c( min_proposalSD,D=proposalSD$min$lnD)
        max_proposalSD <- c( max_proposalSD,D=proposalSD$max$lnD)
        if (length(min_proposalSD) != 4 ||length(max_proposalSD) != 4) stop()

        re <- .Call("ReduceDmvn",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(para$max_aG),
                    as.numeric(para$mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.numeric(evec_sigma_beta),

                    as.numeric(para$mu_lnD),  ## REAL
                    as.numeric(para$sd_lnD),
                    
                    as.integer(M),
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(probIndex),
                    as.integer(thin),
                    as.integer(thinned.sample),
                    as.double(min_proposalSD),
                    as.double(max_proposalSD),
                    package = "lmeNBBayes"
                    )
      
        if (thinned.sample){
          B <- (B - burnin)/thin
        }
        for ( i in 13 : 14 ) re[[i]] <- matrix(re[[i]],nrow=B,ncol=Npat,byrow=TRUE)
        names(re) <- c("aGs","rGs","vs","weightH1",
                       "condProb","h1s","g1s",
                       "beta",
                       "D","logL",
                       "AR","prp","aGs_pat","rGs_pat")
        
        ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
        for ( i in 1 : 4 ) re[[i]] <- matrix(re[[i]],nrow=B,ncol=M,byrow=TRUE)
        re[[5]]  <- matrix(re[[5]],nrow=(Btilde-burnin)/thin,ncol=Npat,byrow=TRUE)
        for ( i in 6 : 7 ) re[[i]] <- matrix(re[[i]],nrow=B,ncol=Npat,byrow=TRUE)
        re[[8]] <- matrix(re[[8]],nrow=B,ncol= pCov,byrow=TRUE)## ncol= pCov or pCov + 1
        if (probIndex)
          {
            ## patients with no new scans
            if (sum(patwoNorO < 0) > 1) patwoNorO <- NULL
            re$condProb[,patwoNorO] <- NA
          }else{
            re$condProb <- NULL
          }
        names(re$prp) <- names(re$AR) <-c("aG", "rG","beta","D")

      }else{ ## NOT DP

   
        re <- .Call("Beta1reduce",
                    as.numeric(Y),           ## REAL
                    as.numeric(X),           ## REAL
                    as.integer(mID),         ## INTEGER
                    as.integer(B),           ## INTEGER
                    as.integer(maxni),       ## INTEGER
                    as.integer(Npat),        ## INTEGER
                    as.numeric(labelnp),     ## REAL
                    as.numeric(para$max_aG),
                    as.numeric(para$mu_beta),     ## REAL
                    as.numeric(evalue_sigma_beta),  ## REAL
                    as.numeric(Inv_sigma_beta),  ## REAL
                    as.numeric(evec_sigma_beta),
                    as.integer(burnin),      ## INTEGER
                    as.integer(printFreq),
                    as.integer(probIndex),
                    as.integer(thin),
                    as.integer(thinned.sample),
                    as.double(min_proposalSD),
                    as.double(max_proposalSD),
                    package = "lmeNBBayes"
                    )
        if (thinned.sample){
          B <- (B - burnin)/thin
        }
       
        re[[3]] <- matrix(re[[3]],nrow=B,ncol= Npat, byrow=TRUE) ## 
        re[[4]] <- matrix(re[[4]],nrow=B,ncol=pCov,byrow=TRUE)
        re[[7]]  <- matrix(re[[7]],nrow=(Btilde-burnin)/thin,ncol=Npat,byrow=TRUE)
        names(re) <- c("aG","rG",
                       "g1s",
                       "beta",
                       "AR","prp",
                       "condProb",
                       "logL"
                       )
        if (probIndex)
          {
            ## patients with no new scans
            if (sum(patwoNorO < 0) > 1) patwoNorO <- NULL
            re$condProb[,patwoNorO] <- NA
          }
        names(re$prp) <- names(re$AR) <-c("aG", "rG","beta")
      }
    if (thinned.sample){
      thin <- 1
      burnin <- 0
    }
    
    if (probIndex)
      {
        re$para$labelnp <- labelnp
        re$condProbSummary <- condProbCI(ID,re$condProb)
        rownames(re$condProbSummary) <- uniID
      }
    re$para <- para
    names(re$para$mu_beta) <- rownames(re$para$Sigma_beta)<-
      colnames(re$para$Sigma_beta ) <- colnames(re$beta) <- covariatesNames
    re$para$CEL <- Y
    re$para$ID <- temID ## return original IDs
    re$para$X <- Xmodmat
    if (DP) re$para$M <- M
    re$para$B <- B
    re$para$burnin <- burnin
    re$para$thin <- thin
    re$para$probIndex <- probIndex
    re$para$burnin <- burnin
    re$para$DP <- DP
    re$para$thinned.sample <- thinned.sample
    re$para$formula <- formula
    re$para$min_proposalSD <-  min_proposalSD
    re$para$max_proposalSD <-  max_proposalSD
    re$para$proposalSD <-  proposalSD
    class(re) <- "LinearMixedEffectNBBayes"
    return (re)
  }









print.LinearMixedEffectNBBayes <- function(x,...){

    para <- x$para
    if (is.vector(x$beta)) pCov <- 1 else pCov <- ncol(x$beta)
    cat("\n -----Negative binomial mixed effect regression----- ")
    if (para$DP){
      cat("\n ---random effect is from DP mixture of beta distributions---")
    }else{
      cat("\n ---random effect is from a single beta distribution---")
    }

    cat("\n====INPUT==== ")

    cat("\n formula: ");print(para$formula)
    if (is.null(para$thin))
      para$thin <- 1
    cat(paste(" MCMC parameters: B=",para$B,", burnin=",para$burnin,", thin=",para$thin,sep=""))
    ##if (para$Reduce)
    ##  {
        cat("\n The target distribution in McMC is a partially marginalized posterior distribution")
        cat("\n (random effects are integrated out).")
    ##  }
    cat("\n ------------------")
    cat("\n [[Hyperparameters]] ")
    cat("\n [Fixed Effect]       ")
    cat("\n beta ~ multivariate normal with mean and covariance matrix:\n  ")
    muSigma <- data.frame(para$mu_beta,"   |   ",para$Sigma_beta);
    colnames(muSigma) <- c("mean","   |   ","covariance",rep(" ",pCov-1))
    rownames(muSigma) <- names(para$mu_beta)
    print(muSigma)
    
    cat(paste(" [Random effect dist'n] aG,rG ~ Unif(min=",0.5,", max=",para$max_aG,")",sep=""))
    if (para$DP){
      cat(paste("\n [Precision parameter]      D ~ Unif(min=",para$a_D,", max=",para$ib_D,")",sep=""))
      cat(paste("\n [Truncation parameter]     M = ",para$M,sep=""))
    }

    
    if (!is.null(para$thin))
      {
        useSample <- useSamp(B=para$B,thin=para$thin,burnin=para$burnin)
        useSampleAll <- FALSE
      }
 
    cat("\n====OUTPUT====")

    cat("\n [Fixed effect]")
    if (para$DP){
       cat(" and [Precision parameter] \n")
       bb <- cbind(x$beta[useSample,],x$D[useSample])
       row.names =  c( names(para$mu_beta),"D")
    }else{
      cat("\n")
      bb <- cbind(x$beta[useSample,])
       row.names = names(para$mu_beta)
    }
    bsummary <- data.frame(colMeans(bb),apply(bb,2,sd),t(apply(bb,2,quantile,prob=c(0.025,0.975))),
                           row.names =  row.names)
    colnames(bsummary) <- c("Estimate","Std. Error","lower CI","upper CI")
    print(round(bsummary,4))
    cat("-----------------")
    cat("\n Reported estimates are posterior means computed after discarding burn-in")
    if (para$thin==1){
      cat(". No Thinning.")
    }else{
      cat(paste(" and thinning at every ",para$thin," samples.",sep=""))
    }
    cat("\n Posterior mean of the logLikelihood",mean(x$logL[useSample]))
    cat("\n Acceptance Rate in Metropolice Hasting rates")
    cat("\n ",paste(names(x$AR),round(x$AR,4),sep=":"))

    if (para$probIndex)
      {
        cat("\n-----------------\n")
        cat("\n ===== Conditional probability index ======\n")
        prIndex(x)
      }
  }




print.IndexBatch <- prIndex <- function(x,...)
  {
    condProb <- x$condProbSummary
    ID <- x$para$ID
    uniID <- unique(ID)
    CEL <- x$para$CEL
    labelnp <- x$para$labelnp
    ## condProb:: A length(uniID) by 3 matrix, first column: a point estimate,
    ## the second column lower CI, the third column upper CI
    ## index of patient, reorder the patients in the increasing ways of probability index1 
    Suborderps <- order(condProb[,1])
    ## ps[order(ps)[1]] == min(ps)
    ## contains the index (1,2,...,Npat) of pat in decreasing way of index1
    orderps <- uniID[Suborderps]
    ## contains the ID (integer but the increment is not always 1) of pat in increasing way
    r <- NULL
    patwoNO <- patwNew0 <- rep(NA,length(uniID))
    ## m1G_1 must be > 0
    ## \hat{ E(1/G) } - 1
    ii <- rep(NA,length(uniID))
    for (i in 1 : length(uniID) )
      {
        i.orderp <- orderps[i]
        pickedIndex <- ID==i.orderp
        ## vector of length ID, containing TRUE if the position agree with the observations of ID orderps[i]
        pickedpt <- which(uniID==i.orderp)
        ## single element indicating the position of the orderps[i] in uniID
        CELpat <- CEL[pickedIndex]
        labelnppat <- labelnp[pickedIndex]
        CELpat0 <- CELpat[labelnppat==0]
        CELpat1 <- CELpat[labelnppat==1]
        if (length(CELpat0)==0) CELpat0 <- "--"
        if (length(CELpat1)==0)
          {
            CELpat1 <- "--"
            patwNew0[i] <- FALSE
          }else
        patwNew0[i] <- sum(CELpat1,na.rm=TRUE)>0 ## TRUE if patients have sum of CEL counts on new scan greater than zero
        ## Patients with no new/old scans are omitted at the end
        patwoNO[i] <- !is.na(condProb[pickedpt,1]) ## TRUE if patients have a new scan and old scan  
        r <- rbind(r,
                   c(
                     "pre-Scan"=paste(CELpat0,collapse="/"),
                     "new-Scan"=paste(CELpat1,collapse="/")
                     )
                   )
        ## ======== printable format =========;
        ii[i] <- paste(sprintf("%1.3f",condProb[pickedpt,1]),
                         " (",sprintf("%1.3f",condProb[pickedpt,3]),
                         ",",sprintf("%1.3f",condProb[pickedpt,4]),")",sep="")
      }
    
    outp <- cbind(ii,r)
    colnames(outp) <- c("conditional probability","pre-scan","new-scan")
    rownames(outp) <- orderps

    cat("\n Estimated conditional probability index and its 95 % CI.")
    if (!is.null(x$para$qsum))
      cat("\n The scalar function to summarize new scans:",x$para$qsum)
    cat("\n-----------------\n")
    if (sum(patwoNO&patwNew0) > 0) print(data.frame(outp[patwoNO&patwNew0,,drop=FALSE]))
    else cat("\n No patient has CPI < 1")
    cat("\n-----------------")
    cat("\n Patients are ordered by decreasing conditional probability.")
    cat("\n Patients with no pre-scans or no new scans are not reported.")
    cat("\n Patients whose new scans are all zero are not reported as conditional probability of such patients are always one.\n")
  }
    


UsedPara <- function(olmeNBB,getSample=FALSE)
  {
    ## If getSample=TRUE, then outputDPFit must be the output of getSample with mod=4/3
    if (!getSample){
      o <- olmeNBB$para
    }else{
      o <- olmeNBB
      o$DP <- TRUE
    }

    UsedPara <- c(
                  "$B$"=o$B,burnin=o$burnin,
                  maxD=o$max_aG,
                  "$\\mu_{\\beta}$"=paste("(",paste(sprintf("%1.2f",c(o$mu_beta[1],o$mu_beta[2])),collapse=","),"...)",sep="",collapse=""),
                  "$\\Sigma_{\\beta}$"=paste("(",paste(sprintf("%1.2f",c(o$Sigma_beta[1,1],o$Sigma_beta[2,1])),collapse=","),
                    "...)")
                  )
    
      
    if (o$DP){
      UsedPara<-c(
                  UsedPara,"M"=o$M,
                  "$ib_D$"=sprintf("%1.2f",o$ib_D),"$a_D$"=sprintf("%1.2f",o$a_D),M=o$M
                  )
    }
    
    return (UsedPara)
    
  }




plotbeta <- function(shape1,shape2,poi=FALSE,col="red",lwd=1,lty=2)
  {
    xs <- seq(0,1,0.001)
    if (!poi)plot(xs,dbeta(xs,shape1,shape2),type="l",col=col,lwd=lwd,lty=lty)
    else points(xs,dbeta(xs,shape1,shape2),type="l",col=col,lwd=lwd,lty=lty)
  }
plotlnorm <- function(mulog,sdlog,poi=FALSE,col="red")
  {
    xs <- seq(0,10,0.001)
    if (!poi) plot(xs,dlnorm(xs,mulog,sdlog),type="l",col=col)
    else points(xs,dlnorm(xs,mulog,sdlog),type="l",col=col)
  }


useSamp <- function(thin,burnin,B)
  {
    temp <- (1:B)*thin
    unadj <- temp[burnin < temp & temp <= B]
    start <- unadj[1] - burnin -1
    return (unadj-start)
  }

dqmix <- function(weightH1,aGs,rGs,alphas=seq(0,0.99,0.01),dens=TRUE)
  {
    ## return a B by length(alphas) matrix, containing the quantile estimates at each corresponding alphas  
    ## weightH1s: B by M matrix
    ## aGs: B by M matrix
    ## rGs: B by M matrix
    B <- nrow(aGs)
    M <- ncol(aGs)
    if(nrow(weightH1) != B) stop()
    if(ncol(weightH1) != M)stop()
    if(nrow(rGs) != B) stop()
    if(ncol(rGs) != M) stop()
    pis <- c(t(weightH1))
    aGs <- c(t(aGs))
    rGs <- c(t(rGs))
    
    if (dens)
      {
        re <- .Call("mixdDist",
                    as.numeric(pis),
                    as.numeric(aGs),
                    as.numeric(rGs),as.numeric(alphas),
                    as.integer(B),as.integer(M),
                    package ="lmeNBBayes")
      }else{
        re <- .Call("mixQDist",
                    as.numeric(pis),
                    as.numeric(aGs),
                    as.numeric(rGs),as.numeric(alphas),
                    as.integer(B),as.integer(M),
                    package ="lmeNBBayes")
      }
    for (i in 1 : length(re)) re[[i]] <- matrix(re[[i]],nrow=B,ncol=length(alphas),byrow=TRUE)
    return (re)
  }

######## subroutines for getS ##############

## Functions to compute E(Yij) and Var(Yij) assuming regression coefficients are random
momBeta <- function(aG,rG,mu.alpha=1.298,sigma.alpha=0.262,dig="%1.2f")
  {
    ## Y | G, alpha ~ NB(size=exp(alpha),prob=G)
    ## G ~ Beta(aG,rG)
    ## alpha ~ N(mu.alpha,sigma.alpha)

    ## E(Y) = E(E(E(Y|G,alpha))) = E(E( [1/G - 1]exp(alpha))) = [ E(1/G)-1 ]*E(exp(alpha))
    ## E(1/G) = ...some calculations... = (aG+rG-1)/(aG-1)
    ## exp(alpha) ~ logN(mu.alpha,sigma.alpha)

    ## E(Y) = ( (aG+rG-1)/(aG-1)-1 )*exp(mu.alpha + sigma.alpha^2/2 )
    E1G <- (aG+rG-1)/(aG-1)
    EexpAlpha <- exp(mu.alpha + sigma.alpha^2/2 )
    EY <- ( E1G -1 )*EexpAlpha

    ## E(1/G^2) = (aG+rG-1)*(aG+rG-2)/((aG-1)*(aG-2))
    E1G2 <- (aG+rG-1)*(aG+rG-2)/( (aG-1)*(aG-2))
    ## Var( exp(alpha)) = (exp(sigma.alpha^2) - 1)*exp(2*mu.alpha + sigma.alpha^2)
    VexpAlpha <- (exp(sigma.alpha^2) - 1)*exp(2*mu.alpha + sigma.alpha^2)
    ## Var(Y) = Var_G( E(Y|G) ) +  E_G( Var(Y|G) )
    ##        = Var_G ( E_alpha ( E(Y|G,alpha))) + E_G( Var_alpha(E(Y|G,alpha)) + E_G( E_alpha( Var(Y|G,alpha))  )
    ##        = Var_G ( 1/G - 1)*E( exp(alpha) )^2 + E_G(Var_alpha( [1/G - 1]exp(alpha) )) + E_G( E_alpha( exp(alpha)*([1-G]/ G^2 )   ))
    ##        = Var_G ( 1/G )*E( exp(alpha) )^2 + E_G([1/G - 1]^2)*Var_alpha( exp(alpha)) + E_G([1-G]/(G^2)) *E_alpha( exp(alpha) )
    ##        = (E(1/G^2)-E(1/G)^2)*E( exp(alpha) )^2
    ##        + [E(1/G^2) - 2*E(1/G) + 1 ]*Var( exp(alpha))
    ##        + [E(1/G^2)-E(1/G)]*E( exp(alpha) )

    VarY <- (E1G2-E1G^2)*((EexpAlpha)^2) +
      (E1G2 - 2*E1G + 1)*VexpAlpha +
        (E1G2 - E1G)*EexpAlpha

    return(list(EY=sprintf(dig,EY),VarY=sprintf(dig,VarY),SEY=sprintf(dig,sqrt(VarY))))
  }


momGL <- function(EG=1.820,VG=0.303,mu.alpha=1.298,sigma.alpha=0.262,dig="%1.2f")
  {
    ## Y|G,alpha ~ NB(size=exp(alpha),prob=1/(G+1))
    ## G ~ dist(mean=EG,var=VG)
    ## alpha ~ N(mu.alpha,sigma.alpha)

    EG2 <- VG + EG^2
    EexpAlpha <- exp(mu.alpha + sigma.alpha^2/2 )
    VexpAlpha <- (exp(sigma.alpha^2) - 1)*exp(2*mu.alpha + sigma.alpha^2)
    ## E(Y) = E(E(E(Y|G,alpha))) = E(E( G*exp(alpha))) = E(G)*E(exp(alpha))
    EY <- EG*EexpAlpha
    ## Var(Y) = Var_G( E(Y|G) ) +  E_G( Var(Y|G) )
    ##        = Var_G ( E_alpha ( E(Y|G,alpha))) + E_G( Var_alpha(E(Y|G,alpha)) + E_G( E_alpha( Var(Y|G,alpha))  )
    ##        = Var_G ( G*E( exp(alpha) ) ) + E_G( Var_alpha( G*exp(alpha) )) + E_G( E_alpha( exp(alpha)*G*(G+1)   ))
    ##        = Var_G ( G )*E(exp(alpha))^2 + E_G(G^2)*Var( exp(alpha) ) + E_G( G*(G+1))*E( exp(alpha) )
    ##        = Var(G)*E(exp(alpha))^2    + E(G^2)*Var( exp(alpha) ) + ( E(G^2) + E(G) )*E( exp(alpha) )
    VarY <-     VG*(EexpAlpha^2) +    EG2*VexpAlpha + (EG2 + EG)*EexpAlpha
  
    return(list(EY=sprintf(dig,EY),VarY=sprintf(dig,VarY),SEY=sprintf(dig,sqrt(VarY))))
  }

##

F_inv_beta2 <- function(ps,aG1,rG1,aG2,rG2,pi)
  {
    quant <- rep(NA,length(ps))
    for (i in 1 : length(ps))
      {
        p <- ps[i]
        G = function (t) pi*pbeta(t,shape1=aG1,shape2=rG1)+(1-pi)*pbeta(t,shape1=aG2,shape2=rG2) - p
        quant[i] <- uniroot(G,c(0,1000))$root
      }   
    return(quant)
  }


index.b.each <- function(Y,ID,labelnp,X,betas,aGs,rGs,pis)
  {
    M <- length(aGs)
    ## betas can be vector and X can be a matrix
    if (is.vector(X)) X <- matrix(X,ncol=1);
    ## compute the conditional probability of
    ## observing Y_{i,new+} >= y_{i,new+} given observing Y_{i,pre+}=y_{i,pew+}
    uniID <- unique(ID)
    Npat <- length(uniID)
    probs <- rep(NA,Npat)
    for (ipat in 1 : Npat)
      {
        Yi <- Y[ID==uniID[ipat]]
        Xi <- X[ID==uniID[ipat],,drop=FALSE]
        labelnpi <- labelnp[ID==uniID[ipat]]
        ## step 1: compute y_{i,pre+} and y_{i,new+} for all the patients:
        Ypresum <- sum(Yi[labelnpi==0])
        Ynewsum <- sum(Yi[labelnpi==1])
        if (sum(labelnpi==0)==0 | sum(labelnpi==1)==0 )
          {
            ## no new scans or no pre scans; no way to calculate the prob index!
            probs[ipat] <- NA
            next;
          }
        if (Ynewsum==0) ## If y_{i,new+} = 0 then the conditional prob of interest is one for that patient
          {
            probs[ipat] <- 1;
            next;
          }
        ## step 2: compute X_{ij}^T beta for all i and j
        rnewsum <- 0
        rpresum <- 0
        for (ivec in 1 : sum(ID==uniID[ipat]))
          {
            rij <- exp(sum(Xi[ivec,]*betas)) ## exp(Xij%*%beta)
            if (labelnpi[ivec]==0) rpresum <- rpresum + rij     ## sum_j exp(Xij%*%beta)
            else if (labelnpi[ivec]==1) rnewsum <- rnewsum + rij
          }
        num  <- 0
        ## step 3: compute the numerator sum_{k=0}^{y_{i,new+}-1} k
        for (k in 0:(Ynewsum-1))
          {
            BoverB <- 0
            for (ih in 1 : M )
              {
                if (pis[ih] == 0) next
                BoverB <- BoverB + pis[ih]*exp(lbeta(rpresum+rnewsum+aGs[ih],k+Ypresum+rGs[ih])-lbeta(aGs[ih],rGs[ih]))
              }
            
            num <- num + exp(log.choose(rnewsum,n=k))*BoverB
            ##num <- num + gamma(rnewsum+k)/(gamma(rnewsum)*gamma(k+1))*BoverB
          }
        ## step 4: compute denominator 
        den <- 0
        for (ih in 1 : M )
          {
            if (pis[ih] == 0) next
            den <- den + pis[ih]*exp( lbeta(rpresum+aGs[ih],Ypresum+rGs[ih]) - lbeta(aGs[ih],rGs[ih]) )
          }
        probs[ipat] <- 1- num/den
      }
    return(probs)
  }


numfun.norm <- function(g,t,Rnp,Rpp,Ypp,ph,ah,rh)
  {
    dnbinom(t,size=Rnp,prob=1/(1+g))*dnbinom(Ypp,size=Rpp,prob=1/(1+g))*ph*dgamma(g,shape=ah,scale=rh)
  }

denfun.norm <- function(g,Ypp,Rpp,ph,ah,rh)
  {
    dnbinom(Ypp,size=Rpp,prob=1/(1+g))*ph*dgamma(g,shape=ah,scale=rh)
  }

index.gammas <- function(Y,ID,labelnp,X,betas,
                         shapes, ## shape
                         scales, ## scale
                         pis)
  {
    M <- length(shapes)
    ## betas can be vector and X can be a matrix
    if (is.vector(X)) X <- matrix(X,ncol=1);
    ## compute the conditional probability of
    ## observing Y_{i,new+} >= y_{i,new+} given observing Y_{i,pre+}=y_{i,pew+}
    uniID <- unique(ID)
    Npat <- length(uniID)
    probs <- rep(NA,Npat);
    for (ipat in 1 : Npat)
      {
        pick <- ID==uniID[ipat]
        Yi <- Y[pick]
        Xi <- X[pick,,drop=FALSE]
        labelnpi <- labelnp[pick]
        ## step 1: compute y_{i,pre+} and y_{i,new+} for all the patients:
        Ypresum <- sum(Yi[labelnpi==0])
        Ynewsum <- sum(Yi[labelnpi==1])
        if (sum(labelnpi==0)==0 || sum(labelnpi==1)==0 )
          {
            ## no new scans or no pre scans; no way to calculate the prob index!
            probs[ipat] <- NA
            next;
          }
        if (Ynewsum==0) ## If y_{i,new+} = 0 then the conditional prob of interest is one for that patient
          {
            probs[ipat] <- 1;
            next;
          }
        ## step 2: compute X_{ij}^T beta for all i and j
        rnewsum <- 0
        rpresum <- 0
        for (ivec in 1 : sum(ID==ipat))
          {
            rij <- exp(sum(Xi[ivec,]*betas))
            if (labelnpi[ivec]==0) rpresum <- rpresum + rij
            else if (labelnpi[ivec]==1) rnewsum <- rnewsum + rij
          }
        num  <- 0
        ## step 3: compute the numerator sum_{k=0}^{y_{i,new+}-1} k
        for (k in 0:(Ynewsum-1))
          {
            for (ih in 1 : M )
              {
                if (pis[ih] == 0) next
                num <- num + integrate(f=numfun.norm, lower=0, upper=Inf,
                                       t=k,Rnp=rnewsum,Rpp=rpresum,Ypp=Ypresum,ph=pis[ih],
                                       ah=shapes[ih],rh=scales[ih])$value


              }
          }
        ## step 4: compute denominator 
        den <- 0
        for (ih in 1 : M )
          {
            if (pis[ih] == 0) next
            den <- den +integrate(f=denfun.norm, lower=0, upper=Inf,
                                  Rpp=rpresum,Ypp=Ypresum,ph=pis[ih],
                                  ah=shapes[ih],rh=scales[ih])$value


          }
        probs[ipat] <- 1 - num/den
      }
    return(probs)
  }




numfun.YZ <- function(g,t,Rnp,Rpp,Ypp,ph,mu,sd)
  {
    dnbinom(t,size=Rnp,prob=1/(1+g))*
      dnbinom(Ypp,size=Rpp,prob=1/(1+g))*
        ph*dnorm(g,mean=mu,sd=sd)
  }

denfun.YZ <- function(g,Ypp,Rpp,ph,mu,sd)
  {
    dnbinom(Ypp,size=Rpp,prob=1/(1+g))*ph*dnorm(g,mean=mu,sd=sd)
  }

index.YZ <- function(Y,ID,labelnp,X,betas,
                     shape, ## shape
                     scale, ## scale
                     mu,sd,pi)
  {
    ## betas can be vector and X can be a matrix
    if (is.vector(X)) X <- matrix(X,ncol=1);
    ## compute the conditional probability of
    ## observing Y_{i,new+} >= y_{i,new+} given observing Y_{i,pre+}=y_{i,pew+}
    uniID <- unique(ID)
    Npat <- length(uniID)
    probs <- rep(NA,Npat);
    for (ipat in 1 : Npat)
      {
        Yi <- Y[ID==uniID[ipat]]
        Xi <- X[ID==uniID[ipat],,drop=FALSE]
        labelnpi <- labelnp[ID==uniID[ipat]]
        ## step 1: compute y_{i,pre+} and y_{i,new+} for all the patients:
        Ypresum <- sum(Yi[labelnpi==0])
        Ynewsum <- sum(Yi[labelnpi==1])
        if (sum(labelnpi==0)==0 || sum(labelnpi==1)==0 )
          {
            ## no new scans or no pre scans; no way to calculate the prob index!
            probs[ipat] <- NA
            next;
          }
        if (Ynewsum==0) ## If y_{i,new+} = 0 then the conditional prob of interest is one for that patient
          {
            probs[ipat] <- 1;
            next;
          }
        ## step 2: compute X_{ij}^T beta for all i and j
        rnewsum <- 0
        rpresum <- 0
        for (ivec in 1 : sum(ID==ipat))
          {
            rij <- exp(sum(Xi[ivec,]*betas))
            if (labelnpi[ivec]==0) rpresum <- rpresum + rij
            else if (labelnpi[ivec]==1) rnewsum <- rnewsum + rij
          }
        num  <- 0
        ## step 3: compute the numerator sum_{k=0}^{y_{i,new+}-1} k
        for (k in 0:(Ynewsum-1))
          {
            num <- num + integrate(f=numfun.norm, lower=0, upper=Inf,
                                   t=k,Rnp=rnewsum,Rpp=rpresum,Ypp=Ypresum,ph=pi,
                                   ah=shape,rh=scale)$value +
                                     integrate(f=numfun.YZ,lower=0,upper=Inf,
                                               t=k,Rnp=rnewsum,Rpp=rpresum,Ypp=Ypresum,ph=1-pi,
                                               mu=mu,sd=sd)$value
          }
        
        ## step 4: compute denominator 
        den <- 0
        
        den <- integrate(f=denfun.norm, lower=0, upper=Inf,
                         Rpp=rpresum,Ypp=Ypresum,ph=pi,
                         ah=shape,rh=scale)$value + integrate(f=denfun.YZ, lower=0, upper=Inf,
                                    Rpp=rpresum,Ypp=Ypresum,ph=1-pi,
                                    mu=mu,sd=sd)$value
        probs[ipat] <- 1 - num/den
      }
    return(probs)
  }


rmvnorm <- 
  function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
            method = c("eigen", "svd", "chol"), pre0.9_9994 = FALSE) 
{
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  sigma1 <- sigma
  dimnames(sigma1) <- NULL
  if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
    warning("sigma is numerically not symmetric")
  }
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
      t(ev$vectors)
  }
  else if (method == "svd") {
    sigsvd <- svd(sigma)
    if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  }
  else if (method == "chol") {
    retval <- chol(sigma, pivot = TRUE)
    o <- order(attr(retval, "pivot"))
    retval <- retval[, o]
  }
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% 
    retval
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}


getS.StatInMed <- function(iseed = "random",
                           rev = 4,
                           dist = "b",
                           mod = 0,
                           probs = seq(0,0.99,0.01),
                           ts = seq(0.001,0.99,0.001),
                           trueCPI = FALSE,
                           full=FALSE,
                           Scenario = "SPMS"
                           )
  {
    ## mod = 0: generate sample
    ## mod = 1: quantiles of the true populations at given probs
    ## mod = 2: densities of the true populations
    ## mod = 3: parameters of the simulation model
    ## dist = "b","b2","YZ" 
    piTRT <- 0.6736842
    ## if trueCPI == TRUE then returns the most precise conditional probability computed based on the
    ## treatment assignment 
    Npat <- 180; ## upto review 4, Npat=160 in total this number must be divisible by NpatEnterPerMonth

    if (iseed=="random") set.seed(sample(1e+6,1)) else  set.seed(iseed)

    ## The prior for regression coefficients beta

    ## (Intercept) trt010:timeInt1 trt011:timeInt1 trt010:timeInt2 trt011:timeInt2
    ## Estimated based on the 5 IFN-beta trials

    if (Scenario == "full"){
      ## para.full$mu_beta
      mu_beta <- round(c(1.3091326,  0.0150677, -0.7759003, -0.1715055, -1.0394682
                         ## StatInMedFirstSubmission
                         ## 1.3103072,  0.0145406, -0.7758956, -0.1718227, -1.0401282
                         ),6)
      
      ## para.full$Sigma_beta
      Sigma_beta <- matrix(round(c(
                                   0.06093404,    -0.02348381,    -0.05339049,    -0.02575333,    -0.07351022,
                                  -0.02348381,     0.01698531,     0.01337642,     0.01842227,     0.02813052,
                                  -0.05339049,     0.01337642,     0.49999949,     0.04797195,     0.74974055,
                                  -0.02575333,     0.01842227,     0.04797195,     0.02716636,     0.07697036,
                                  -0.07351022,     0.02813052,     0.74974055,     0.07697036,     1.15908970
                                   ## StatInMedFirstSubmission
                                   ##   0.06001687,    -0.02345988,    -0.05568177,    -0.02575713,    -0.07730916,
                                   ## -0.02345988,     0.01706351,     0.01419726,     0.01859172,     0.02936096,
                                   ## -0.05568177,     0.01419726,     0.50037189,     0.04904017,     0.75136313,
                                   ## -0.02575713,     0.01859172,     0.04904017,     0.02754559,     0.07884069,
                                   ## -0.07730916,     0.02936096,     0.75136313,     0.07884069,     1.16294508
                                   ),6)
                           ,nrow=length(mu_beta),byrow=TRUE)
      ## para.full$mu_lnD
      mu_lnD <- -0.7731277 ##StatInMedFirstSubmission -0.7565044 
      ## para.full$sd_lnD
      sd_lnD <-0.2559623   ##StatInMedFirstSubmission 0.2608056   
    }else if (Scenario == "SPMS"){
      ##para.SPMS$mu_beta
      mu_beta <- round(c(1.25548248,  0.09307017, -0.89594757, -0.01817658, -1.05446830
                         ## StatInMedFirstSubmission
                         ## 1.25680296,  0.09376186, -0.89599224, -0.01692496, -1.05552658
                         ),6)
      ## para.SPMS$Sigma_beta
      Sigma_beta <- matrix(round(c(0.017934818,    0.006469281,    -0.03534636,   -0.008756420,    -0.03344118,
                                   0.006469281,    0.017793152,    -0.08635236,   -0.009483972,    -0.11631098,
                                  -0.035346364,   -0.086352363,     0.56018192,    0.070573793,     0.77596421,
                                  -0.008756420,   -0.009483972,     0.07057379,    0.013904651,     0.09328106,
                                  -0.033441184,   -0.116310975,     0.77596421,    0.093281064,     1.10824255
                                   ## StatInMedFirstSubmission
                                   ## 0.017715774,    0.006411256,    -0.03561045,   -0.008621458,    -0.03419703,
                                   ##  0.006411256,    0.017846706,    -0.08644089,  -0.009355102,    -0.11660896,
                                   ## -0.035610451,   -0.086440894,     0.56035742,    0.070343557,     0.77681689,
                                   ## -0.008621458,   -0.009355102,     0.07034356,    0.013857689,     0.09321587,
                                   ## -0.034197026,  -0.116608960,     0.77681689,    0.093215871,     1.10993034
                                   ),6),
                           nrow=length(mu_beta),byrow=TRUE)
      ## para.SPMS$mu_lnD
      mu_lnD <- -0.4823671  ## StatInMedFirstSubmission-0.5175115
      ## para.SPMS$sd_lnD
      sd_lnD <- 0.3556396 ## StatInMedFirstSubmission 0.3617086 
    }
    if (mod == 3){
      ## return parameter information of selected model
      outputMod3 <- list(mu_beta=mu_beta, Sigma_beta = Sigma_beta,mu_lnD=mu_lnD,sd_lnD=sd_lnD,Scenario=Scenario)
      ## Compute the E(Yij) at baseline for patients from each RE cluster 
      mu.alpha <- mu_beta[1]
      sigma.alpha <- Sigma_beta[1,1]
    }
    
    if (dist == "b"){
      aG1 <- 3  
      rG1 <- 0.8 
      if (mod %in% c(0,4:6)){
        gs <- rbeta(Npat,aG1,rG1)
        hs <- rep(1,Npat)
      }else if (mod==1){
        ## mod = 1: quantiles of the true populations at given probs
        return(cbind(probs=probs,
                     quantile=qbeta(probs,shape1=aG1,shape2=rG1)))
      }else if (mod==2){
        ## mod = 2: densities of the true populations
        return (cbind(ts=ts,
                      dens=dbeta(ts,shape1=aG1,shape2=rG1)) )
      }else if (mod == 3){
        ## mod = 3: parameters of the simulation model
        c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
        outputMod3 <- c(outputMod3,list(K=1,c1=c1,aGs=c(aG1=aG1), rGs=c(rG1=rG1)))
      }

    }else if (dist == "b2"){
      ## mixture of two beta distributions
      ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
      ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
      pi <- 0.3
      
      aG1 <- 10
      rG1 <- 10
      aG2 <- 20
      rG2 <- 1
      ## generate the initial random effect values of everyone
      if (mod %in% c(0,4:6)){
        hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
        Npat_dist1 <- sum(hs==1)
        Npat_dist2 <- sum(hs==2)
        
        gs <- rep(NA,Npat)
        gs[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
        gs[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
      }else if (mod==1){
        return (cbind(probs=probs,
                      quantile=F_inv_beta2(ps=probs,aG1,rG1,aG2,rG2,pi=pi)))
      }else if (mod==2){
        return (cbind(ts=ts,
                      dens=(pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2))))
      }else if (mod == 3){

        c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
        c2 <- momBeta(aG2,rG2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
        outputMod3 <- c(outputMod3,list(
                                        K=2,c1=c1,c2=c2,
                                        aGs=c(aG1=aG1,aG2=aG2),
                                        rGs=c(rG1=rG1,rG2=rG2),
                                        pi=pi))
      }
    }else if ( dist == "YZ"){
      
      pi <- 0.85

      alpha <- exp(-0.5)
      ## a bimodal distribution with 85 % of Gi from a gamma distribution with mean 0.647 and variance 2.374
      scale <- 2.374/0.647*alpha 
      shape <- 0.647^2/2.374
      mu <- 3*alpha
      sd <- sqrt(0.25)*alpha
      
      ## generate the initial random effect values of everyone
      if (mod %in% c(0,4:6) ){
        hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
        Npat_dist1 <- sum(hs==1)
        Npat_dist2 <- sum(hs==2)
        
        gs <- rep(NA,Npat)
        gs[hs==1] <- rgamma(Npat_dist1,shape=shape,scale=scale);
        gs[hs==2] <- rnorm(Npat_dist2,mean=mu,sd=sd);
        gs[gs < 0 ] <- 0
        gs <- 1/(1+gs)
        
      }else if (mod==1){
        return ("this option is currently not available")
      }else if (mod == 2){
        ts.trans <-1/ts-1
        return (cbind(ts=ts,
                      (pi*dgamma(ts.trans,shape=shape,scale=scale)+(1-pi)*dnorm(ts.trans,mean=mu,sd=sd) )*(1/ts^2)))
        
      }else if (mod == 3){
        c1 <- momGL(EG=shape*scale,VG=shape*(scale^2),mu.alpha=mu.alpha,sigma.alpha=sigma.alpha) ## expectation and variance of Y
        c2 <- momGL(EG=mu,VG=sd^2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
        outputMod3 <- outputMod3 <- c(outputMod3,list(c1=c1,c2=c2,scale=scale,shape=shape,mu=mu,sd=sd,pi=pi,K=2))
      }
    }
    ## return parameter information of selected model
    if (mod == 3) return (outputMod3)
    
    ## Generate regression coefficients of this dataset
    betFull <- matrix(rmvnorm(n=1,mean=mu_beta,
                                sigma=Sigma_beta),ncol=1)
    trtAss.pat <- sample(0:1,Npat,c(1-piTRT,piTRT),replace=TRUE)
    ## Sequential samples
    ni <- 10
    NpatEnterPerMonth <- 15
    DSMBVisitInterval <- 4 ## months
    
    ## === necessary for mod= 0, 3 and mod = 4
    ## d contains a full dataset
    days <- NULL
    for (ipatGroup in 1 : (Npat/NpatEnterPerMonth))
      {
        ScandaysForSingleGroup <- ipatGroup:(ipatGroup+ni-1)
        days <- c(days,rep(ScandaysForSingleGroup,NpatEnterPerMonth))
      }
    Y <- XforFit <- XforTCPI <- NULL
    ## timeInt for all patients (used in model fit)
    Xagg <- cbind(rep(1,ni),
                  c(rep(0,2),rep(1,4),rep(0,ni-6)),
                  c(rep(0,6),rep(1,ni-6)))
    zeros <- rep(0,ni)
    Xplcb <- cbind(Xagg[,1], ## ni = 9
                   Xagg[,2],
                   zeros, 
                   Xagg[,3],
                   zeros)
    Xtrt <- cbind(Xagg[,1], ## ni = 9
                  zeros, 
                  Xagg[,2],
                  zeros,
                  Xagg[,3])
    for ( ipat in 1 : Npat)
      {
        ## placebo
        if (trtAss.pat[ipat]==0){
          X <- Xplcb 
        }else if (trtAss.pat[ipat]==1) X <- Xtrt ## trt
        ## the number of repeated measures are the same
        ## we assume that the time effects occurs once after the treatments are in effect
        got <- rnbinom(ni,size = exp(X%*%betFull), prob = gs[ipat])         
        Y <- c(Y,got)
        ## XforFit and XforTCPI must be the same if piTRT = 0
        XforFit <- rbind(XforFit,Xagg)
        if (trueCPI) XforTCPI <- rbind(XforTCPI,X)
      }
    colnames(XforFit) <- c("Intercept","timeInt1","timeInt2")

    trtlabel <- factor(rep(trtAss.pat,each=ni),labels=c("plcb","trt"))

    d <- data.frame(Y=Y,
                    XforFit,
                    ID=rep(1:Npat,each=ni),
                    gs = rep(gs,each=ni),
                    scan = rep(-1:(ni-2),Npat),
                    ## day contains the day when the scan was taken
                    ## 10 patients enter a trial every month
                    days = days,
                    hs = rep(hs,each=ni),
                    trtAss = trtlabel)

    if (trueCPI){
      d$XforTCPI <- XforTCPI
    }
    if (full) return (d) 
    ## DSMB visit is assumed to be every 4 months
    d <- subset(d,subset= days <= DSMBVisitInterval*rev)
    d$labelnp <- rep(0,nrow(d))
    d$labelnp[ DSMBVisitInterval*(rev-1) < d$days ] <- 1
    ## The first two scans (screening and base-line scans are treated as pre-scans)
    d$labelnp[ d$scan <= 0 ] <- 0
    betPlcb <- betFull[c(1,2,4)]
    ## cat("betFull",betFull)
    if (mod == 0){
      XX <- cbind(d$Intercept,d$timeInt1, d$timeInt2)
      if (dist=="b")
        {
          temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
                               betas=betPlcb,aGs=aG1,rGs=rG1,pis=1)
          if (trueCPI){
            tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
                                 betas=betFull,aGs=aG1,rGs=rG1,pis=1)
          }
        }else if (dist=="b2"){
          temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
                               betas=betPlcb,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))  
          if (trueCPI){
            tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
                                 betas=betFull,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))
          }
        }else if (dist == "YZ"){
          
          temp <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
                           betas=betPlcb,
                           shape=shape, ## shape
                           scale=scale, ## scale
                           mu=mu,sd=sd,pi=pi)
          if (trueCPI){
            tCPI <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
                             betas=betFull,
                             shape=shape, ## shape
                             scale=scale, ## scale
                             mu=mu,sd=sd,pi=pi)
          }
        }
      
      d$betPlcb <- c(betPlcb,rep(NA,nrow(d)-length(betPlcb)))
      d$betFull <- c(betFull,rep(NA,nrow(d)-length(betFull)))
      d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp) ))
      if (trueCPI){
        d$probIndexTRUE <- c(tCPI,rep(NA,nrow(d)-length(tCPI) ))
      }
      return(d)
    }
  }





## Used in the first submission to Stat In Med
## getSNormUnblind <- function(iseed = "random",
##                             rev = 4,
##                             dist = "b",
##                             mod = 0,
##                             probs = seq(0,0.99,0.01),
##                             ts = seq(0.001,0.99,0.001),
##                             trueCPI = FALSE,
##                             full=FALSE,
##                             Scenario = "SPMS"
##                             )
##   {
##     ## mod = 0: generate sample
##     ## mod = 1: quantiles of the true populations at given probs
##     ## mod = 2: densities of the true populations
##     ## mod = 3: parameters of the simulation model
##     ## dist = "b","b2","YZ" 
##     piTRT <- 0.6736842
##     ## if trueCPI == TRUE then returns the most precise conditional probability computed based on the
##     ## treatment assignment 
##     Npat <- 180; ## upto review 4, Npat=160 in total this number must be divisible by NpatEnterPerMonth

##     if (iseed=="random") set.seed(sample(1e+6,1)) else  set.seed(iseed)

##     ## The prior for regression coefficients beta

##     ## (Intercept) trt010:timeInt1 trt011:timeInt1 trt010:timeInt2 trt011:timeInt2
##     ## Estimated based on the 5 IFN-beta trials

##     if (Scenario == "full"){
##       ## para.full$mu_beta
##       mu_beta <- round(c(1.3103072,  0.0145406, -0.7758956, -0.1718227, -1.0401282
##                          ##1.31044123,  0.01501563, -0.77578748, -0.17174723, -1.04176877
##                          ),4)
      
##       ## para.full$Sigma_beta
##       Sigma_beta <- matrix(round(c(
##                                    0.06001687,    -0.02345988,    -0.05568177,    -0.02575713,    -0.07730916,
##                                  -0.02345988,     0.01706351,     0.01419726,     0.01859172,     0.02936096,
##                                  -0.05568177,     0.01419726,     0.50037189,     0.04904017,     0.75136313,
##                                  -0.02575713,     0.01859172,     0.04904017,     0.02754559,     0.07884069,
##                                  -0.07730916,     0.02936096,     0.75136313,     0.07884069,     1.16294508),6)
##                            ,nrow=length(mu_beta),byrow=TRUE)
                           
##                            ## round(c( 0.05220635, -0.02079985, -0.05010230, -0.02297177, -0.06964878,
##                            ##         -0.02079985,  0.01395803,  0.01255353,  0.01580000,  0.02625352,
##                            ##         -0.05010230,  0.01255353,  0.44997693,  0.04404106,  0.67839680,
##                            ##         -0.02297177,  0.01580000,  0.04404106,  0.02295083,  0.07111575,
##                            ##         -0.06964878,  0.02625352,  0.67839680,  0.07111575,  1.05203182),4)
##                            ## ,nrow=length(mu_beta),byrow=TRUE)
##       ## para.full$mu_lnD
##       mu_lnD <- -0.7565044 ##-0.7539186
##       ## para.full$sd_lnD
##       sd_lnD <- 0.2608056   #9.318458e-05
##     }else if (Scenario == "SPMS"){
##       ## para.SPMS$full
##       mu_beta <- round(c(1.25680296,  0.09376186, -0.89599224, -0.01692496, -1.05552658),6)
##       ##1.25549440, 0.09461850, -0.89751584, -0.01645465, -1.05669044),4)
##       Sigma_beta <- matrix(round(c( 0.017715774,    0.006411256,    -0.03561045,   -0.008621458,    -0.03419703,
##                                     0.006411256,    0.017846706,    -0.08644089,  -0.009355102,    -0.11660896,
##                                    -0.035610451,   -0.086440894,     0.56035742,    0.070343557,     0.77681689,
##                                    -0.008621458,   -0.009355102,     0.07034356,    0.013857689,     0.09321587,
##                                    -0.034197026,  -0.116608960,     0.77681689,    0.093215871,     1.10993034),6),
##                            nrow=length(mu_beta),byrow=TRUE)
##       mu_lnD <-  -0.5175115 ##-0.5117792 
##       sd_lnD <-  0.3617086 ## 5.39951e-06
##     }
##     if (mod == 3){
##       ## return parameter information of selected model
##       outputMod3 <- list(mu_beta=mu_beta, Sigma_beta = Sigma_beta,mu_lnD=mu_lnD,sd_lnD=sd_lnD,Scenario=Scenario)
##       ## Compute the E(Yij) at baseline for patients from each RE cluster 
##       mu.alpha <- mu_beta[1]
##       sigma.alpha <- Sigma_beta[1,1]
##     }
    
##     if (dist == "b"){
##       aG1 <- 3  
##       rG1 <- 0.8 ##0.8
##       if (mod %in% c(0,4:6)){
##         gs <- rbeta(Npat,aG1,rG1)
##         hs <- rep(1,Npat)
##       }else if (mod==1){
##         ## mod = 1: quantiles of the true populations at given probs
##         return(cbind(probs=probs,
##                      quantile=qbeta(probs,shape1=aG1,shape2=rG1)))
##       }else if (mod==2){
##         ## mod = 2: densities of the true populations
##         return (cbind(ts=ts,
##                       dens=dbeta(ts,shape1=aG1,shape2=rG1)) )
##       }else if (mod == 3){
##         ## mod = 3: parameters of the simulation model
##         c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- c(outputMod3,list(K=1,c1=c1,aGs=c(aG1=aG1), rGs=c(rG1=rG1)))
##       }

##     }else if (dist == "b2"){
##       ## mixture of two beta distributions
##       ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##       ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##       pi <- 0.3
      
##       aG1 <- 10
##       rG1 <- 10
##       aG2 <- 20
##       rG2 <- 1
##       ## generate the initial random effect values of everyone
##       if (mod %in% c(0,4:6)){
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gs <- rep(NA,Npat)
##         gs[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
##         gs[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
##       }else if (mod==1){
##         return (cbind(probs=probs,
##                       quantile=F_inv_beta2(ps=probs,aG1,rG1,aG2,rG2,pi=pi)))
##       }else if (mod==2){
##         return (cbind(ts=ts,
##                       dens=(pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2))))
##       }else if (mod == 3){

##         c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         c2 <- momBeta(aG2,rG2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- c(outputMod3,list(
##                                         K=2,c1=c1,c2=c2,
##                                         aGs=c(aG1=aG1,aG2=aG2),
##                                         rGs=c(rG1=rG1,rG2=rG2),
##                                         pi=pi))
##       }
##     }else if ( dist == "YZ"){
      
##       pi <- 0.85

##       alpha <- exp(-0.5)
##       ## a bimodal distribution with 85 % of Gi from a gamma distribution with mean 0.647 and variance 2.374
##       scale <- 2.374/0.647*alpha 
##       shape <- 0.647^2/2.374
##       mu <- 3*alpha
##       sd <- sqrt(0.25)*alpha
      
##       ## generate the initial random effect values of everyone
##       if (mod %in% c(0,4:6) ){
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gs <- rep(NA,Npat)
##         gs[hs==1] <- rgamma(Npat_dist1,shape=shape,scale=scale);
##         gs[hs==2] <- rnorm(Npat_dist2,mean=mu,sd=sd);
##         gs[gs < 0 ] <- 0
##         gs <- 1/(1+gs)
        
##       }else if (mod==1){
##         return (NULL)
##       }else if (mod == 2){
##         ts.trans <-1/ts-1
##         return (cbind(ts=ts,
##                       (pi*dgamma(ts.trans,shape=shape,scale=scale)+(1-pi)*dnorm(ts.trans,mean=mu,sd=sd) )*(1/ts^2)))
        
##       }else if (mod == 3){
##         c1 <- momGL(EG=shape*scale,VG=shape*(scale^2),mu.alpha=mu.alpha,sigma.alpha=sigma.alpha) ## expectation and variance of Y
##         c2 <- momGL(EG=mu,VG=sd^2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- outputMod3 <- c(outputMod3,list(c1=c1,c2=c2,scale=scale,shape=shape,mu=mu,sd=sd,pi=pi,K=2))
##       }
##     }
##     ## return parameter information of selected model
##     if (mod == 3) return (outputMod3)
    
##     ## Generate regression coefficients of this dataset
##     betFull <- matrix(rmvnorm(n=1,mean=mu_beta,
##                                 sigma=Sigma_beta),ncol=1)
##     trtAss.pat <- sample(0:1,Npat,c(1-piTRT,piTRT),replace=TRUE)
##     ## Sequential samples
##     ni <- 10
##     NpatEnterPerMonth <- 15
##     DSMBVisitInterval <- 4 ## months
    
##     ## === necessary for mod= 0, 3 and mod = 4
##     ## d contains a full dataset
##     days <- NULL
##     for (ipatGroup in 1 : (Npat/NpatEnterPerMonth))
##       {
##         ScandaysForSingleGroup <- ipatGroup:(ipatGroup+ni-1)
##         days <- c(days,rep(ScandaysForSingleGroup,NpatEnterPerMonth))
##       }
##     Y <- XforFit <- XforTCPI <- NULL
##     ## timeInt for all patients (used in model fit)
##     Xagg <- cbind(rep(1,ni),
##                   c(rep(0,2),rep(1,4),rep(0,ni-6)),
##                   c(rep(0,6),rep(1,ni-6)))
##     zeros <- rep(0,ni)
##     Xplcb <- cbind(Xagg[,1], ## ni = 9
##                    Xagg[,2],
##                    zeros, 
##                    Xagg[,3],
##                    zeros)
##     Xtrt <- cbind(Xagg[,1], ## ni = 9
##                   zeros, 
##                   Xagg[,2],
##                   zeros,
##                   Xagg[,3])
##     for ( ipat in 1 : Npat)
##       {
##         ## placebo
##         if (trtAss.pat[ipat]==0){
##           X <- Xplcb 
##         }else if (trtAss.pat[ipat]==1) X <- Xtrt ## trt
##         ## the number of repeated measures are the same
##         ## we assume that the time effects occurs once after the treatments are in effect
##         got <- rnbinom(ni,size = exp(X%*%betFull), prob = gs[ipat])         
##         Y <- c(Y,got)
##         ## XforFit and XforTCPI must be the same if piTRT = 0
##         XforFit <- rbind(XforFit,Xagg)
##         if (trueCPI) XforTCPI <- rbind(XforTCPI,X)
##       }
##     colnames(XforFit) <- c("Intercept","timeInt1","timeInt2")

##     trtlabel <- factor(rep(trtAss.pat,each=ni),labels=c("plcb","trt"))

##     d <- data.frame(Y=Y,
##                     XforFit,
##                     ID=rep(1:Npat,each=ni),
##                     gs = rep(gs,each=ni),
##                     scan = rep(-1:(ni-2),Npat),
##                     ## day contains the day when the scan was taken
##                     ## 10 patients enter a trial every month
##                     days = days,
##                     hs = rep(hs,each=ni),
##                     trtAss = trtlabel)

##     if (trueCPI){
##       d$XforTCPI <- XforTCPI
##     }
##     if (full) return (d) 
##     ## DSMB visit is assumed to be every 4 months
##     d <- subset(d,subset= days <= DSMBVisitInterval*rev)
##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ DSMBVisitInterval*(rev-1) < d$days ] <- 1
##     ## The first two scans (screening and base-line scans are treated as pre-scans)
##     d$labelnp[ d$scan <= 0 ] <- 0
##     betPlcb <- betFull[c(1,2,4)]
    
##     if (mod == 0){
##       XX <- cbind(d$Intercept,d$timeInt1, d$timeInt2)
##       if (dist=="b")
##         {
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                                betas=betPlcb,aGs=aG1,rGs=rG1,pis=1)
##           if (trueCPI){
##             tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                                  betas=betFull,aGs=aG1,rGs=rG1,pis=1)
##           }
##         }else if (dist=="b2"){
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                                betas=betPlcb,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))  
##           if (trueCPI){
##             tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                                  betas=betFull,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))
##           }
##         }else if (dist == "YZ"){
          
##           temp <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                            betas=betPlcb,
##                            shape=shape, ## shape
##                            scale=scale, ## scale
##                            mu=mu,sd=sd,pi=pi)
##           if (trueCPI){
##             tCPI <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                              betas=betFull,
##                              shape=shape, ## shape
##                              scale=scale, ## scale
##                              mu=mu,sd=sd,pi=pi)
##           }
##         }
      
##       d$betPlcb <- c(betPlcb,rep(NA,nrow(d)-length(betPlcb)))
##       d$betFull <- c(betFull,rep(NA,nrow(d)-length(betFull)))
##       d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp) ))
##       if (trueCPI){
##         d$probIndexTRUE <- c(tCPI,rep(NA,nrow(d)-length(tCPI) ))
##       }
##       return(d)
##     }
##   }



getTrueProb <- function(aG1=0.5,rG1=0.5,aG2=3,rG2=0.8,dist="gamma")
  {
    
    if (dist=="beta")
      {
        ## compute P(G1 < G2) = int_{y1=0}^1 (int_{y2=y1}^1 beta(y2;a2,b2) dy2 )*beta(y1;a1,b1) dy1
        insideIntBeta <- function(y1,a2,b2,a1,b1)
          integrate(dbeta,shape1=a2,shape2=b2,lower=y1, upper=1)$value*dbeta(y1,shape1=a1,shape2=b1)
        
        return (integrate(insideIntBeta,a2=aG2,b2=rG2,a1=aG1,b1=rG1,lower=0,upper=1)$value)
        
      }else if (dist=="gamma")
        {
          
          ## compute P(G1 < G2) = int_{y2=0}^Inf (int_{y1=0}^y1 gamma(y1;a1,b1) dy1 gamma(y2;a2,b2) dy2
          insideIntGam <- function(y2,a2,b2,a1,b1)
            {
              integrate(dgamma,shape=a1,scale=b1,lower=0, upper=y2)$value*dgamma(y2,shape=a2,scale=b2)
            }
          return  (integrate(insideIntGam,
                             a2=aG2,b2=rG2,a1=aG1,b1=rG1,lower=0,upper=Inf)$value)
        }
  }


mixgamma <- function(sp1,sc1,sp2,sc2,xs = seq(0,15,0.001),ylim=c(0,2),points=FALSE)
  {
    dg1 <- dgamma(xs,shape=sp1,scale=sc1)
    dg2 <- dgamma(xs,shape=sp2,scale=sc2)
    ys <- (1/2)*dg1+(1/2)*dg2
    if (!points){
      plot(xs,ys,ylim=ylim,type="l",
           main=paste("The mixture of gamma distributions:\n 1/2*dgamma(shape=",
             sp1,",scale=",sc1,")+1/2*dgamma(shape=",sp2,",scale=",sc2,")",sep=""))
      points(xs,dg1,type="l",col="red")
      points(xs,dg2,type="l",col="blue")
    }else{
      points(xs,ys,type="l",col="blue")
    }
  }


getM <- function(epsilonM,a_D,r_D,prob)
  {
    ## returns the appropriate truncation number M
    Dbig <- qgamma(prob,shape=a_D,scale=1/r_D)
    return(round(1+ log(epsilonM)/(log(Dbig/(1+Dbig))),0))
  }

Nuniq <-function(test,sampleindex=NULL){
  if (is.null(sampleindex)){
    tt <- apply(test$h1,1,unique)
  }else{
    tt <- apply(test$h1[sampleindex,],1,unique)
  }
  return( unlist(lapply(tt,length)))
}

colmeansd <- function (mat, name = NULL, sigDig = 2, sigDigSD = NULL,space=" ",na.rm=FALSE) 
{
    if (is.vector(mat)) 
        mat <- matrix(mat, ncol = 1)
    if (length(sigDig) == 1) 
        sigDig <- rep(sigDig, ncol(mat))
    else if (length(sigDig) != ncol(mat)) 
        stop("length(sigDig) != ncol(mat)")
    
    if (length(sigDigSD) == 0)
      sigDigSD <- sigDig
    else if (length(sigDigSD) == 1) 
      sigDigSD <- rep(sigDigSD, ncol(mat))
    else if (length(sigDigSD) != ncol(mat)) 
      stop("length(sigDigSD) != ncol(mat)")
    
    mea <- colMeans(mat,na.rm=na.rm)
    sdd <- apply(mat, 2, sd,na.rm=na.rm)
    re <- NULL
    MysigDig <- paste("%1.", sigDig, "f", sep = "")
    MysigDigSD <- paste("%1.", sigDigSD, "f", sep = "")
    for (i in 1:ncol(mat)) {
        if (is.na(mea[i])) {
            re <- c(re, "------")
        }
        else {
            re <- c(re, paste(sprintf(MysigDig[i], mea[i]),space, "(", 
                sprintf(MysigDigSD[i], sdd[i], 3), ")", sep = ""))
        }
    }
    if (!is.null(name)) 
        names(re) <- name
    if (is.null(name) & !is.null(colnames(mat))) 
        names(re) <- colnames(mat)
    return(re)
}


pointsgamma<-function(shape,scale,xmin=0,xmax=5,main="")
  {
    ## plotting the density of gamma dist'n with given shape and scale
    xs <-seq(xmin,xmax,length.out=1000)
    ys <- dgamma(xs,shape=shape,scale=scale)
    points(xs,ys,main=main,type="l",col="blue")
  }

plotgamma<-function(shape,scale,xmin=0,xmax=5,main="")
  {
    ## plotting the density of gamma dist'n with given shape and scale
    xs <-seq(xmin,xmax,length.out=1000)
    ys <- dgamma(xs,shape=shape,scale=scale)
    plot(xs,ys,main=main,type="l")
  }


plotnbinom <- function(size,prob,xmin=0,xmax=30,main="",size2=NULL,prob2=NULL,xlab="")
  {
    ## plotting the density of gamma dist'n with given size and prob
    xs <-xmin:xmax
    ys <- dnbinom(xs,size=size,prob=prob)
    plot(xs,ys,main=main,type="b",xlab=xlab)
    if (!is.null(size2))
      {
        ys <- dnbinom(xs,size=size2,prob=prob2)
        points(xs,ys,col="red",type="b")
      }
  }

plotHistAcf<- function(postsample,up.burnin=3000,mainplot,mainhist,xlim=NULL,col="black",everythin=1)
  {
    plot(postsample,type="l",main=mainplot,col=col)
    samp <- postsample[-(1:up.burnin)]
    samp <- samp[(1:(length(samp)/everythin))*everythin]
    med <- median(samp)
    if (is.null(xlim))
      {
        hist(samp,main=paste(mainhist," median:",round(med,3),sep=""))
      }else{
        hist(samp,main=paste(mainhist," median:",round(med,3),sep=""),
             xlim=xlim)##,breaks=0:1000)
      }
    ci <- round(quantile(samp,prob=c(0.025,0.975)),2)
    points(x=ci,y=c(0,0),type="b",col="blue",lwd=3)
    
    abline(v=med,col="blue")
    acf(samp,main="")
  }




plotGs <- function(vec,main,xlim=c(0,8),upto=10,length.out=1000,breaks=seq(0,15,0.001),ylim)
  {
    med <- round(median(vec),3)
    mea <- round(mean(vec),3)
    hist(vec,probability=TRUE,
         main=paste(main," mean: ",mea," median: ",med,sep=""),
         xlim=xlim,breaks=c(breaks,100000),ylim=ylim)
    ## density estimates are obtained at seq(0,upto,length.out=length.out)
    estden <- density(vec,kernel = "gaussian",from=0,to=upto,n=length.out)
    points(estden$x,estden$y,type="l")
  }





## ##=== Meta analysis function === ##
## rigamma <- function (n, shape, scale) 
##   {
##     if (shape > 0 & scale > 0) 
##       1/rgamma(n = n, shape=shape, scale=1/scale)
##     else stop("rigamma: invalid parameters\n")
##   }

## metaUnivGibb <- function(ys,
## ### a vector of length Nstudy, containing the observed estimates
##                          sigmas2,
## ### a vector of length Nstudy, containing the number of patients (or should it be the number of MS scan?) involved to estimate ys
##                          a=1,
##                          b=1,
##                          B=100000)
##   {
##     Nstudy <- length(ys)
    
##     ## Model:
##     ## ys[istudy] ~ N(mean=theta[istudy],sd=sqrt(sigmas2[istudy]))
##     ## theta_i ~ N(mu,tau2)
##     ## mu ~ unif(-1000,1000)
##     ## tau2 ~ inv-gamma(shape=a,scale=b)
##     ## where i = 1,...,Nstudy and Nstudy is the number of dataset

##     ## In our MS clinical analysis context,
##     ## y_i = hat.beta_i (i=1,..,p) or hat.logaG_i
    
##     ## informative prior construction
##     ## Y ~ inv-gamma(shape=a,scale=b) then E(Y)=b/(a-1) for a > 1, Var(Y)=b^2/((a-1)^2(a-2)) for a > 2
##     ## a <= 2 => infinite variance 
##     postThetas <- matrix(NA,B,Nstudy)
##     postMu <- rep(NA,B)
##     postTau2 <- rep(NA,B)
##     ## Initialization of unknown parameters
##     thetas <- rep(0,Nstudy)
##     mu <- mean(ys)
##     tau2 <- mean(sigmas2)
    
##     for (iB in 1 : B)
##       {
##         if (iB %% 5000 == 0) cat("\n",iB, "iterations are done" )
##         ## Step 1: update theta_i
##         for (istudy in 1 : Nstudy)
##           {
##             vari <- 1/(1/sigmas2[istudy]+1/tau2)
##             thetas[istudy]<- rnorm(1,
##                                    mean=(ys[istudy]/sigmas2[istudy] + mu/tau2)*vari,
##                                    sd=sqrt(vari)
##                                    )
##           }
        
##         mu <- rnorm(1,
##                     mean=sum(thetas)/Nstudy,
##                     sd=sqrt(tau2/Nstudy))
        
##         tau2 <- rigamma(
##                         1,
##                         shape=(Nstudy/2+a),
##                         scale=sum((thetas-mu)^2)/2+b
##                         )
        
##         postThetas[iB,] <- thetas
##         postMu[iB] <- mu
##         postTau2[iB] <- tau2
##       }
##     return(list(theta=postThetas,mu=postMu,tau2=postTau2))
##   }





## listMean <- function(listobj,useSample)
##   {
##     if (is.null(useSample))
##       {
##         useSample <- 1: length(listobj)
##       }
##     output <- listobj[[useSample[1]]]
##     for (i in useSample[-1])
##       {
##         output <- output +listobj[[i]]
##       }
##     return (output/length(useSample))
##   }



## metaMultGibb <- function(hat.betas,
##                          ## p by N.study matrix, where p = # covariates 
##                          sigmas2,
##                          ## a list of covariance matrix 
##                          v0="uninfo",
##                          S0=NULL,
##                          B=10000)
##   {
##     ## TheWishart distribution is a multivariate analogue of the gamma distribution
##     Nstudy <- ncol(hat.betas)
##     p <- nrow(hat.betas)
##     sigmas2Inv <- lapply(sigmas2,solve)
##     ## Model:
##     ## hat.betas[[istudy]] ~ mvnorm(mean=betas[,istudy], var=sigmas2[[istudy]])
##     ## betas ~ mvnorm(mean=mu, var=tau2)
##     ## :::: prior :::: 
##     ## mu[i] ~ unif(-1000,1000), i = 1,...,p
##     ## tau2 ~ inv-wishart(v0,S0) <=> solve(tau2) ~ wishart(v0,solve(S0))
##     ## where i = 1,...,Nstudy and Nstudy is the number of dataset

##     ## In our MS clinical analysis context,
##     ## y_i = hat.beta_i (i=1,..,p) or hat.logaG_i


    
##     ## The default values of hyperparameters
##     ## inverse-wisrt is centered at diag(p) with large variance
##     ## p109 Hogg
##     if (v0=="info"){
##       v0 <- 10
##       S0 <- diag(p)*(v0-p-1) ## Identity matrix
      
##     }else if (v0=="uninfo"){
##       v0 <- p + 2
##       S0 <- diag(c(0.5,0.10)) ## Identity matrix
##     }

    
##     ## storages
##     postBetas <- postTau2 <- list()
##     postMu <- matrix(NA,p,B)

##     ## Initialization of unknown parameters
##     betas <- matrix(0,nrow=p,ncol=Nstudy)
##     mu <- rowMeans(hat.betas)
##     tau2 <- listSum(sigmas2)/Nstudy
##     tau2Inv <- solve( tau2 )
    
##     for (iB in 1 : B)
##       {
##         if (iB %% 5000 == 0) cat("\n",iB, "iterations are done" )
##         ## Step 1: update betas[,istudy]
##         for (istudy in 1 : Nstudy)
##           {
##             vari <- solve(tau2Inv+sigmas2Inv[[istudy]] )
##             betas[,istudy]<- mvrnorm(1,
##                                      mu=vari%*%( tau2Inv%*%mu+sigmas2Inv[[istudy]]%*%hat.betas[,istudy]),
##                                      Sigma=vari
##                                      )
##           }

##         mu <- mvrnorm(1,
##                       mu=rowSums(betas)/Nstudy, 
##                       Sigma=tau2/Nstudy
##                       )


##         ## temp <- 0; for (istudy in 1 : Nstudy) temp <- temp + (betas[,istudy]-mu)%*%t(betas[,istudy]-mu)
##         ## generate sample from Inverse-Wishart distribution
##         Sn <- S0 + (betas-mu)%*%t(betas-mu) ## outer product sum_{istudy} (betas[,istudy]-mu)%*%t(betas[,istudy]-mu)
##         tau2Inv <- matrix(rWishart(n=1,
##                                    df=v0+Nstudy,
##                                    Sigma=solve(Sn)
##                                    ),p,p)

##         postMu[,iB] <- mu
##         postTau2[[iB]] <- tau2 <- solve(tau2Inv)
##         postBetas[[iB]] <- betas
##       }
##     return(list(betas=postBetas,mu=postMu,tau2=postTau2))
##   }









## nbinREDPmix <- function(Y,   ##     A vector of length sum ni, containing responses 
##                         X,   ##     A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                         ID,  ##     A Vector of length sum ni, indicating patients
##                         B = 20000,   ##     A scalar, the number of Gibbs iteration 
##                         M = NULL,    ##     A scalar, the truncation value
##                         labelnp,
##                         dist = 1,      ##     dist == 1 means the base distribution of the random effects are assumed to be 
##                         max_aG = 3.0,    ##     from P0 = gamma(shape=aG,scale=1/rG)
##                         a_r = 3.0,       ##      where aG ~ Unif(0,max_aG ) and 
##                         r_r = 0.2,    ##     rG ~ gamma(shape=a_r,scale=1/r_r)
##                         a_D = 0.5,    ##     small shape a.D and small scale ib.D put weights on both small and large values of D
##                         r_D = 2,      ##      D ~ gamma(shape=a.D,scale=1/r.D)
##                         a_qG = 20.0,     ##     qG ~ beta(shape1=a_qG,scale=r_qG)
##                         r_qG = 1.0,
##                         mu_beta = 0,   ##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta)
##                         sigma_beta = 5,##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta) 
##                         burnin = round(B/3),
##                         epsilonM = 0.001,
##                         prob = 0.999,
##                         printfFreq = B
##                         )
##   {

##     ## ==== Description of the parameters in code ====
##     ## g1s:            A vector of length N, containing the random effect 1 values for each pat
##     ## g2s:            A vector of length N, containing the random effect 2 values for each pat
##     ## vs:             A vector of length M, containing the parameters to construct pis
##     ## weightH1(vs):   A vector of length M, containing the probabilities of the categorical distribution of Hi, function of vs
##     ## beta:           A vector of length p, containing the coefficients
##     ## h1s:            A vector of length N, containing the cluster labels of RE1 of each observation, cluster index is between 1 and M+1
##     ## h2s:            A vector of length N, containing the cluster labels of RE2 of each observation, cluster index is between 1 and M+1
##     ## aGs:            A vector of length M, containing the shape parameters of the gamma distribution
##     ## rGs:            A vector of length M, containing the scale parameters of the gamma distribution

##     ## ## === Model ===
##     ## // 1) y_ij | Gij = gij,beta, ~  NB(r_ij,p_i)
##     ## //    where r_ij = exp(t(x_ij)*beta), p_{ij} = 1/(g_ij+1)
##     ## //    gij  = g1i   if j is in pre scan
##     ## //         = g2i   if j is in new scan
##     ## //    and the distribution of g2i is specified as;
##     ## //    g2i  = g1i   with prob qG2
##     ## //         = g*    with prob 1-qG2   (The default value of p.g2 = 0)
##     ## //
##     ## // 2) Gi|ih ~ gamma(shape=aGs[ih],scale=1/rGs[ih])
##     ## //  aGs[ih],ih=1,...,M ~ Unif(0,max_aG )
##     ## //  rGs[ih],ih=1,...,M ~ gamma(shape=a_r,scale=1/r_r)
##     ## //
##     ## // 3) beta[i] ~ normal(0,sigma_beta) i=1,..,p
##     ## // 4) D ~ gamma(shape=a.D,scale=b.D)
##     if (is.null(M)) M = getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)

##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }
##     ## The patients with labelnp = 1 (no old scans) for all repeated measures
##     ## are treated as lack of new scans and used to estimate beta only.
##     ## skip the computation of H2
##     ## All patients have old scans

##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

##     patwoONscan <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1) 
##     if (length(patwoONscan)==0) patwoONscan <- -999;
##     p <- ncol(X)
##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
##     mID <- ID-1
##     mpatwoONscan <- patwoONscan - 1
##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))

##     re <- .Call("res",
##                 as.numeric(Y),           ## REAL
##                 as.numeric(X),           ## REAL
##                 as.integer(mID),         ## INTEGER
##                 as.integer(B),           ## INTEGER
##                 as.integer(maxni),       ## INTEGER
##                 as.integer(M),           ## INTEGER
##                 as.integer(Npat),        ## INTEGER
##                 as.numeric(labelnp),     ## REAL
##                 as.numeric(max_aG),      ## REAL
##                 as.numeric(a_r),         ## REAL
##                 as.numeric(r_r),         ## REAL
##                 as.numeric(a_D),         ## REAL
##                 as.numeric(r_D),         ## REAL
##                 as.numeric(a_qG),        ## REAL 
##                 as.numeric(r_qG),        ## REAL
##                 as.numeric(mu_beta),     ## REAL
##                 as.numeric(sigma_beta),  ## REAL
##                 as.integer(mpatwoONscan),## INTEGER
##                 as.integer(dist),        ## INTEGER
##                 as.integer(burnin),      ## INTEGER
##                 as.integer(printfFreq),
##                 package = "lmeNBBayes"
##                 )

##     for ( i in 1 : 5) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)

##     for ( i in 6 : 9) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##     re[[10]] <- matrix(re[[10]],B,p,byrow=TRUE)

##     names(re) <- c("h1","h2","g1","g2","j",
##                    "aG","rG","v","weightH1",
##                    "beta","D","qG","AR","prp")

##     para <- list(mu_beta = mu_beta,
##                  sigma_beta = sigma_beta,
##                  M = M,
##                  a_D = a_D,
##                  r_D = r_D,
##                  dist = dist,
##                  burnin = burnin,
##                  a_r = a_r,
##                  r_r = r_r,
##                  a_qG = a_qG,
##                  r_qG = r_qG,
##                  max_aG = max_aG,
##                  patwoONscan=(patwoONscan+1))

##     re$para <- para
##     names(re$AR) <-c("gs", "aG", paste("beta",1:p))
##     re$hs <- re$h1 + 1  
##     re$h2s <- re$h2 + 1 
##     return(re)  
##   }



## nbinREDPmix2 <- function(Y,   ##     A vector of length sum ni, containing responses 
##                          X,   ##     A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                          ID,  ##     A Vector of length sum ni, indicating patients
##                          B = 20000,   ##     A scalar, the number of Gibbs iteration 
##                          M = NULL,    ##     A scalar, the truncation value
##                          labelnp,
##                          dist = 1,    ##     dist == 1 means the base distribution of the random effects are assumed to be 
##                          max_aG = 3.0,    ##     from P0 = gamma(shape=aG,scale=1/rG)
##                          ##a_r = 3.0,       ##      where aG ~ Unif(0,max_aG ) and 
##                          ##r_r = 0.2,    ##     rG ~ gamma(shape=a_r,scale=1/r_r)
##                          a_D = 0.5,    ##     small shape a.D and small scale ib.D put weights on both small and large values of D
##                          r_D = 2,      ##      D ~ gamma(shape=a.D,scale=1/r.D)
##                          a_qG = 20.0,     ##     qG ~ beta(shape1=a_qG,scale=r_qG)
##                          r_qG = 1.0,
##                          mu_beta = 0,   ##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta)
##                          sigma_beta = 5,##     hyperpara of beta: beta ~ norm(mu_beta,sd=sigma_beta) 
##                          burnin = round(B/3),
##                          epsilonM = 0.001,
##                          prob = 0.999,
##                          printfFreq = B
##                          )
##   {

##     ## ==== Description of the parameters in code ====
##     ## g1s:            A vector of length N, containing the random effect 1 values for each pat
##     ## g2s:            A vector of length N, containing the random effect 2 values for each pat
##     ## vs:             A vector of length M, containing the parameters to construct pis
##     ## weightH1(vs):   A vector of length M, containing the probabilities of the categorical distribution of Hi, function of vs
##     ## beta:           A vector of length p, containing the coefficients
##     ## h1s:            A vector of length N, containing the cluster labels of RE1 of each observation, cluster index is between 1 and M+1
##     ## h2s:            A vector of length N, containing the cluster labels of RE2 of each observation, cluster index is between 1 and M+1
##     ## aGs:            A vector of length M, containing the shape parameters of the gamma distribution
##     ## rGs:            A vector of length M, containing the scale parameters of the gamma distribution

##     ## ## === Model ===
##     ## // 1) y_ij | Gij = gij,beta, ~  NB(r_ij,p_i)
##     ## //    where r_ij = exp(t(x_ij)*beta), p_{ij} = 1/(g_ij+1)
##     ## //    gij  = g1i   if j is in pre scan
##     ## //         = g2i   if j is in new scan
##     ## //    and the distribution of g2i is specified as;
##     ## //    g2i  = g1i   with prob qG2
##     ## //         = g*    with prob 1-qG2   (The default value of p.g2 = 0)
##     ## //
##     ## // 2) Gi|ih ~ gamma(shape=aGs[ih],scale=1/rGs[ih])
##     ## //  aGs[ih],ih=1,...,M ~ Unif(0,max_aG )
##     ## //  rGs[ih],ih=1,...,M ~ gamma(shape=a_r,scale=1/r_r)
##     ## //
##     ## // 3) beta[i] ~ normal(0,sigma_beta) i=1,..,p
##     ## // 4) D ~ gamma(shape=a.D,scale=b.D)
##     if (is.null(M)) M = getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)

##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }
##     ## The patients with labelnp = 1 (no old scans) for all repeated measures
##     ## are treated as lack of new scans and used to estimate beta only.
##     ## skip the computation of H2
##     ## All patients have old scans

##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

##     patwoONscan <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1) 
##     if (length(patwoONscan)==0) patwoONscan <- -999;
##     p <- ncol(X)
##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
##     mID <- ID-1
##     mpatwoONscan <- patwoONscan - 1
##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))

##     re <- .Call("resBeta",
##                 as.numeric(Y),           ## REAL
##                 as.numeric(X),           ## REAL
##                 as.integer(mID),         ## INTEGER
##                 as.integer(B),           ## INTEGER
##                 as.integer(maxni),       ## INTEGER
##                 as.integer(M),           ## INTEGER
##                 as.integer(Npat),        ## INTEGER
##                 as.numeric(labelnp),     ## REAL
##                 as.numeric(max_aG),      ## REAL
##                 ##as.numeric(a_r),         ## REAL
##                 ##as.numeric(r_r),         ## REAL
##                 as.numeric(a_D),         ## REAL
##                 as.numeric(r_D),         ## REAL
##                 as.numeric(a_qG),        ## REAL 
##                 as.numeric(r_qG),        ## REAL
##                 as.numeric(mu_beta),     ## REAL
##                 as.numeric(sigma_beta),  ## REAL
##                 as.integer(mpatwoONscan),## INTEGER
##                 as.integer(dist),        ## INTEGER
##                 as.integer(burnin),      ## INTEGER
##                 as.integer(printfFreq),
##                 package = "lmeNBBayes"
##                 )

##     for ( i in 1 : 5) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)

##     for ( i in 6 : 9) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##     re[[10]] <- matrix(re[[10]],B,p,byrow=TRUE)

##     names(re) <- c("h1","h2","g1","g2","j",
##                    "aG","rG","v","weightH1",
##                    "beta","D","qG","AR","prp")

##     para <- list(mu_beta = mu_beta,
##                  sigma_beta = sigma_beta,
##                  M = M,
##                  a_D = a_D,
##                  r_D = r_D,
##                  dist = dist,
##                  burnin = burnin,
##                  a_qG = a_qG,
##                  r_qG = r_qG,
##                  max_aG = max_aG,
##                  patwoONscan=(patwoONscan+1))

##     re$para <- para
##     names(re$AR) <-c("gs", "aG", paste("beta",1:p))
##     re$hs <- re$h1 + 1  
##     re$h2s <- re$h2 + 1 
##     return(re)  
##   }





## getSampleOld <- function(
##                       iseed=911,
##                       rev=4,
##                       dist="b",
##                       mod=0,
##                       beta1=0,
##                       probs=seq(0,0.99,0.01),
##                       ts = seq(0,0.99,0.01)
##                       )
## {
##   ## mod = 0: generate sample
##   ## mod = 1: quantiles of the true populations
##   ## mod = 2: densities of the true populations
##   ni <- 11
##   Npat <- 200;
##   ni_rev <- c(5,8,11)[rev-1]
##   set.seed(iseed)
##   ## generate unobservable latent variables gs:
##   if (dist=="g")
##     {
##       ## r.e. mixture of two gamma distributions
##       ## cluster1: E(Y)=exp(0.5)*E(G) = exp(0.5)*sp1*sc1 = 0.4121803
##       ## cluster2: E(Y)=exp(0.5)*E(G) = exp(0.5)*sp2*sc2 = 3.956931
##       pi <- 0.7
##       temp <- rmultinom(1,Npat,c(pi,1-pi))
##       Npat_dist1 <- temp[1,1]
##       Npat_dist2 <- temp[2,1]
##       sp1 <- 0.5
##       sc1 <- 0.5
##       sp2 <- 3
##       sc2 <- 0.8
##       beta0 <- 0.5
##       sampldist1 <- rgamma(Npat_dist1,shape=sp1,scale=sc1)
##       sampldist2 <- rgamma(Npat_dist2,shape=sp2,scale=sc2)
##       gsBASE <- c(sampldist1,sampldist2)
##       gsBASE = 1/(1+gsBASE)
##       if (mod==1)
##         {
##           F_inv_gam <- function(p,sp1,sc1,sp2,sc2,pi)
##             {
##               G = function (t)
##                 pi*pgamma(t,shape=sp1,scale=sc1)+(1-pi)*pgamma(t,shape=sp2,scale=sc2) - p
##               return(uniroot(G,c(0,100))$root)
##             }
##           lq <- length(probs);
##           quan <- rep(NA,lq);
##           for (i in 1 : lq) quan[i] <- F_inv_gam(p=probs[i],sp1=sp1,sc1=sc1,sp2=sp2,sc2=sc2,pi=pi)
##           quan = sort(1/(1+quan))
##           return (rbind(probs=probs,quantile=quan))
##         }else if (mod == 2){

##             ts.trans <-1/ts-1
##             return ((pi*dgamma(ts.trans,shape=sp1,scale=sc1)+(1-pi)*dgamma(ts.trans,shape=sp2,scale=sc2) )*(1/ts^2))
##           }else if (mod == 3){
##             return (list(beta0=beta0,K=2))
##           }
##     }else if (dist== "l")
##       {
##         ## single log-normal distribution
##         ## E(Y)=exp(beta)*E(G)=exp(1)*exp(meanlog+sdlog^2/2)= 5.078419
##         beta0 <- 1
##         meanlog <- -2.5
##         sdlog <- 2.5
##         gsBASE <- rlnorm(Npat,meanlog=meanlog,sdlog=sdlog)
##         gsBASE = 1/(1+gsBASE)
##         quan <- qlnorm(probs,meanlog=meanlog,sdlog=sdlog)
##         if (mod==1){
##           return( 
##                  rbind(probs=probs,
##                        quantile=sort(1/(1+quan))
##                        )
##                  )
##         }else if (mod==2){

##             ts.trans <-1/ts-1
##             return(dlnorm(ts.trans,meanlog=meanlog,sdlog=sdlog)*(1/ts^2))
##           }else if (mod == 3){
##             return (list(beta0=beta0,K=1))
##           }
##       }
##     else if (dist == "b")
##       {
##         ## single beta distribution
##         ## E(1/G) = ...some calculations... = (aG+rG-1)/(aG-1)
##         ## E(Y)=exp(beta0)*mu_{1/G}=exp(1)*((aG1+rG1-1)/(aG1-1)+1)= 4.238999
##         beta0 <- 1
##         aG1 <- 15.3
##         rG1 <- 8
##         gsBASE <- rbeta(Npat,aG1,rG1)
##         if (mod==1){
##           return( 
##                  rbind(probs=probs,
##                        quantile=qbeta(probs,shape1=aG1,shape2=rG1)
##                        )
##                  )
##         }else if (mod==2){

##             return (dbeta(ts,shape1=aG1,shape2=rG1))

##           }else if (mod == 3){

##             return (list(beta0=beta0,K=1))

##           }
##     }else if (dist == "b2")
##       {
##         ## mixture of two beta distributions
##         ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##         ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##         pi <- 0.3
##         temp <- rmultinom(1,Npat,c(pi,1-pi))
##         Npat_dist1 <- temp[1,1]
##         Npat_dist2 <- temp[2,1]
##         beta0 <- 2
##         aG1 <- 13
##         rG1 <- 18
##         aG2 <- 20
##         rG2 <- 1
##         sampldist1 <- rbeta(Npat_dist1,aG1,rG1);
##         sampldist2 <- rbeta(Npat_dist2,aG2,rG2);
##         gsBASE <- c(sampldist1,sampldist2)
##         if (mod==1){
##           F_inv_beta2 <- function(p,aG1,rG1,aG2,rG2,pi)
##             {
##               G = function (t) pi*pbeta(t,shape1=aG1,shape2=rG1)+(1-pi)*pbeta(t,shape1=aG2,shape2=rG2) - p
##               return(uniroot(G,c(0,1000))$root)
##             }
##           lq <- length(probs);
##           quant <- rep(NA,lq);
##           for (i in 1 : lq) quant[i] <- F_inv_beta2(p=probs[i],aG1,rG1,aG2,rG2,pi=pi)
##           return (rbind(probs=probs,quantile=quant))
##         }
##         else if (mod==2){

##           return ((pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2)) )
##         }else if (mod == 3){
##           return(list(beta0=beta0,K=2))
##         }
##       }
##   else if (dist == "b3")
##     {
##       ## mixture of two beta distributions
##       ## cluster 1: E(Y)=exp(beta)*mu_G=exp(0.5)*((aG1+rG1-1)/(aG1-1)+1)=11.33496
##       ## cluster 2: E(Y)=exp(beta)*mu_G=exp(0.5)*((aG2+rG2-1)/(aG2-1)+1)=5.035284
##       ## cluster 3: E(Y)=exp(beta)*mu_G=exp(0.5)*((aG3+rG3-1)/(aG3-1)+1)=2.000417
##       beta0 <- 1
##       p1 <- 0.4
##       p2 <- 0.1
##       temp <- rmultinom(1,200,c(p1,p2,1-p1-p2))
##       Npat_dist1 <- temp[1,1]
##       Npat_dist2 <- temp[2,1]
##       Npat_dist3 <- temp[3,1]
##       aG1 <- 5
##       rG1 <- 19.5
##       aG2 <- 19.5
##       rG2 <- 19.5
##       aG3 <- 25
##       rG3 <- 0.01
##       sampldist1 <- rbeta(Npat_dist1,aG1,rG1);
##       sampldist2 <- rbeta(Npat_dist2,aG2,rG2);
##       sampldist3 <- rbeta(Npat_dist3,aG3,rG3);
##       gsBASE <- c(sampldist1,sampldist2,sampldist3)
##       if (mod==1)
##         {
##           F_inv_beta3 <- function(p,aG1,rG1,aG2,rG2,aG3,rG3,p1,p2)
##             {
##               G = function (t)
##                 p1*pbeta(t,shape1=aG1,shape2=rG1)+p2*pbeta(t,shape1=aG2,shape2=rG2)+(1-p1-p2)*pbeta(t,shape1=aG3,shape2=rG3)-p
##               return(uniroot(G,c(0,1000))$root)
##             }
##           lq <- length(probs);
##           quant <- rep(NA,lq);
##         for (i in 1 : lq) quant[i] <- F_inv_beta3(p=probs[i],aG1,rG1,aG2,rG2,aG3,rG3,
##                                                   p1=p1,
##                                                   p2=p2)
##           return (rbind(probs=probs,quantile=quant))
##         }else if (mod==2)
##           {
##             return (p1*dbeta(ts,shape1=aG1,shape2=rG1)+p2*dbeta(ts,shape1=aG2,shape2=rG2)+(1-p1-p2)*dbeta(ts,shape1=aG3,shape2=rG3 ))
##           }else if (mod == 3){
##             return( list(beta0=beta0,K=3))
##           }
##     }

##   Y <- NULL
##   for ( ipat in 1 : Npat)
##     {
##       ## the number of repeated measures are the same
##       got <- rnbinom(ni,size = exp(beta0+beta1*(1:ni)), prob = gsBASE[ipat])         
##       Y <- c(Y,got)
##     }
##   cat("beta0",beta0,"\n")

##   d <- data.frame(Y=Y,
##                   X=rep(1,Npat*ni),
##                   ID=rep(1:Npat,each=ni),
##                   gsBASE1 = rep(gsBASE,each=ni),
##                   scan = rep(1:ni,Npat)
##                   )

##   d <- subset(d,subset=scan <= ni_rev)
##   d$labelnp <- rep(c(rep(0,ni_rev-3),rep(1,3)),Npat)

##   return(d)
## }




## DPfitold <- function(Y,   ##   A vector of length sum ni, containing responses 
##                   X,   ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                   ID,  ##   A Vector of length sum ni, indicating patients
##                   B = 105000,     ##     A scalar, the number of Gibbs iteration 
##                   max_aG = 20.0,  ##     from P0 = gamma(shape=aG,scale=1/rG)
##                   mu_beta = NULL,    ##     hyperpara of beta: beta[ibeta] ~ norm(mu_beta[ibeta],sd=sigma_beta[ibeta])
##                   Sigma_beta = NULL, ## Beta24
##                   sigma_beta = NULL, ##    
##                   burnin = 5000,  
##                   printFreq = B,
##                   model=8,
##                   M=NULL,
##                   a_D=1,
##                   r_D=1,
##                   epsilonM=0.01,  ## nonpara
##                   prob=0.9,       ## nonpara
##                   labelnp=NULL, 
##                   sd_eta=0.3,      ## Beta13 and Beta23
##                   initial_beta=0.0,## Beta23
##                   a_qG=1,          ## Beta14 and Beta24
##                   r_qG=1,          ## Beta14 and Beta24
##                   mu_aG = 1.5,        ## Beta24
##                   sd_aG = 0.75,      ## Beta24
##                   mu_rG = 1.5,        ## Beta24
##                   sd_rG = 0.75      ## Beta24
##                   )
##   {
##     Ntot <- length(Y)

##     if (is.vector(X)) X <- matrix(X,ncol=1)
##     pCov <- ncol(X)
##     if (nrow(X)!= Ntot) stop ("nrow(X) != length(Y)")
##     if (length(ID)!= Ntot)  stop ("length(ID)!= length(Y)")

##     if (is.null(sigma_beta)) sigma_beta <- rep(5,pCov)
##     if (is.null(mu_beta)) mu_beta <- rep(0,pCov)
##     if (is.null(labelnp)) labelnp <- rep(0,length(Y))

##     ## if (is.null(M)) M = getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##     ## change the index of ID to numeric from 1 to # patients
##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }

##     p <- pCov
##     if( is.vector(X)) p <- 1 
##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
##     mID <- ID-1

##     ## The patients with labelnp = 1 (no old scans) for all repeated measures
##     ## are treated as lack of new scans and used to estimate beta only.
##     ## skip the computation of H2
##     ## All patients have old scans

##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

##     patwoNS <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
##     if (length(patwoNS)==0) patwoNS <- -999;
##     patwoNS <- patwoNS - 1

##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))
##     if (model==1)
##       {
##         ##parametric model: gi is constant over time 
##         labelnp <- rep(0,length(Y))    
##         re <- .Call("Beta1",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(mu_aG),
##                     as.numeric(sd_aG),
##                     as.numeric(mu_rG),
##                     as.numeric(sd_rG),
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )

##         re[[3]] <- matrix(re[[3]],B,N,byrow=TRUE)
##         re[[4]] <- matrix(re[[4]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","gPre","beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model==2)
##       {
##         ## parametric model: gi is allowed to change between new scan g1 and old scan g2 period
##         ##                   g1 and g2 are assumed to be independent
##         re <- .Call("Beta2",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )

##         for ( i in 3:4 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[5]] <- matrix(re[[5]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","g1s","g2s","beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model==3)
##       {
##         ## Nonparametric model: gi is constant over time 
##         labelnp <- rep(0,length(Y))
##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##         re <- .Call("Beta12",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.numeric(a_D),  ## REAL
##                     as.numeric(r_D),  ## REAL
##                     as.integer(M),
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 5:6 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[7]] <- matrix(re[[7]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","vs","weightH1",
##                        "h1s","g1s",
##                        "beta","D","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG",paste("beta",1:p))
##         return(re)
##       }
##     else if (model==4)
##       {
##         ## nonparametric model: gi is allowed to change between new scan g1 and old scan g2 period
##         ##                      g1 and g2 are assumed to be independent 
##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##         re <- .Call("Beta22",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.numeric(a_D),         ## REAL
##                     as.numeric(r_D),         ## REAL
##                     as.integer(M),           ## INTEGER
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 5:8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[9]] <- matrix(re[[9]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","vs","weightH1",
##                        "h1s","h2s","g1s","g2s",
##                        "beta","D","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG,
##                      M=M,
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model == 5)
##       {
##         min_prec_eta = 500
##         max_prec_eta = 501
##         re <- .Call("Beta13",
##                     as.numeric(Y),           ## 1REAL
##                     as.numeric(X),           ## 2REAL
##                     as.integer(mID),         ## 3INTEGER
##                     as.integer(B),           ## 4INTEGER
##                     as.integer(maxni),       ## 5INTEGER
##                     as.integer(Npat),        ## 6INTEGER
##                     as.numeric(labelnp),     ## 7REAL
##                     as.numeric(min_prec_eta),   ## 8REAL
##                     as.numeric(max_prec_eta),   ## 9REAL
##                     as.numeric(max_aG),      ## 10REAL
##                     as.numeric(mu_beta),     ## 11REAL
##                     as.numeric(sigma_beta),  ## 12REAL
##                     as.integer(burnin),      ## 13INTEGER
##                     as.integer(printFreq),   ## 14INTEGER
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 4 : 6 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[7]] <- matrix(re[[7]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","prec_eta",
##                        "g1s","etas","g2s",
##                        "beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG
##                      )
##         re$para <- para
##         names(re$AR) <-c("g1s","etas","aG", "rG", "prec_eta",paste("beta",1:p))
##         return(re)

##       }
##     else if (model==6)
##       {
##         ##  nonparametric model,  
##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)

##         re <- .Call("Beta23",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.integer(M),           ## INTEGER
##                     as.numeric(a_D),         ## REAL
##                     as.numeric(r_D),         ## REAL
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(sd_eta),      ## REAL
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),   ## INTEGER
##                     as.numeric(initial_beta),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         for ( i in 1 : 4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 5 : 8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[5]] <- re[[5]]+1 ## "h1s"
##         re[[9]] <- matrix(re[[9]],B,p,byrow=TRUE)
##         names(re) <- c("aGs","rGs","vs","weightH1",
##                        "h1s","g1s","g2s","etas"
##                        ,"betas","D",##"prec_eta",
##                        "AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      max_aG = max_aG,
##                      M=M,
##                      ID=ID,
##                      labelnp=labelnp,
##                      patWOnew=patwonew
##                      ) 
##         re$para <- para
##         names(re$AR) <- names(re$prp) <- c("g1s","etas",paste("aG",1:M,sep=""),
##                                            paste("rG",1:M,sep=""),
##                                            paste("beta",1:p,sep=""))##"prec_eta")
##         re$MeanWH <- colMeans(re$weightH1);
##         names(re$MeanWH) <- paste("cluster",1:M)
##         return(re)                
##       }
##     else if (model==7)
##       {

##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ beta(aG,rG)
##         ## // beta ~ rnorm(mu_beta,sigma_beta)
##         ## // aG, rG ~ unif(0,max_aG)

##         re <- .Call("Beta14",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),
##                     as.numeric(r_qG),
##                     as.numeric(max_aG),      ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(sigma_beta),  ## REAL
##                     ##as.numeric(a_D),         ## REAL
##                     ##as.numeric(r_D),         ## REAL
##                     ##as.integer(M),           ## INTEGER
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 4 : 7 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[8]] <- matrix(re[[8]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","qG",
##                        "gPre","g2s","gNew","js",
##                        "beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      sigma_beta = sigma_beta,
##                      burnin = burnin,
##                      B=B,
##                      max_aG = max_aG,
##                      patwoNS= patwoNS+1
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", paste("beta",1:p))
##         return(re)
##       }
##     else if (model==8)
##       {
##         ## nonparametric model
##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ sum_{h=1}^M pi_h beta(aG_h,rG_h)
##         ## // beta ~ rnorm(mu_beta,sigma_beta)
##         ## // aG, rG ~ unif(0,max_aG)
##         ##  nonparametric model,

##         if (is.null(M)) M  <- getM(epsilonM=epsilonM,a_D=a_D,r_D=r_D,prob=prob)
##         if (is.null(Sigma_beta)) Sigma_beta <-  diag(5,pCov)
##         evalue_sigma_beta <- eigen(Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
##         if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
##         Inv_sigma_beta <- c( solve(Sigma_beta) )

##         re <- .Call("Beta24",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.integer(M),           ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),        ## REAL
##                     as.numeric(r_qG),        ## REAL
##                     as.numeric(mu_aG),       ## REAL
##                     as.numeric(sd_aG),       ## REAL
##                     as.numeric(mu_rG),       ## REAL
##                     as.numeric(sd_rG),       ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(evalue_sigma_beta),  ## REAL
##                     as.numeric(Inv_sigma_beta),  ## REAL
##                     as.numeric(a_D),         ## REAL
##                     as.numeric(r_D),         ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     package = "lmeNBBayes"
##                     )
##         for ( i in 2 : 7 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         for ( i in 8 : 11 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         re[[12]] <- matrix(re[[12]],B,p,byrow=TRUE)
##         names(re) <- c("qG",
##                        "gPre","g2s","gNew","js","h1s","h2s",
##                        "weightH1","vs","aGs","rGs",
##                        "beta","AR","prp")
##         re$para <- list(mu_beta = mu_beta,
##                         sigma_beta = (Sigma_beta),
##                         burnin = burnin,
##                         M=M,
##                         B=B,
##                         model=model,
##                         mu_aG = mu_aG,
##                         sd_aG = sd_aG,
##                         mu_rG = mu_rG,
##                         sd_rG = sd_rG,
##                         a_D=a_D,
##                         r_D=r_D,
##                         a_qG=a_qG,
##                         r_qG=r_qG,
##                         Npat=N,
##                         Ntot=length(Y),
##                         patwoNS = patwoNS+1,
##                         CEL=Y,
##                         ID=ID
##                         )
##         names(re$AR) <- names(re$prp) <- c(paste("aG",1:M,sep=""),
##                                            paste("rG",1:M,sep=""),
##                                            paste("beta",sep="")
##                                            )
##         re$h1s <- re$h1s + 1
##         re$h2s <- re$h2s + 1
##         re$MeanWH <- colMeans(re$weightH1);
##         names(re$MeanWH) <- paste("cluster",1:M)
##         return(re)

##       }
##   }


## getSample2 <- function(iseed=1,
##                        mod=1, ## mod==2 returns the simulation parameter 
##                        aG1=10,rG1=10,aG2=20,rG2=1)
## {
##   ## The random effect changes its value over time
##   ## 20 % of patients switch their random effect values between pre and new scan period
##   ## By review 2, 80 patients are recruited to the study 
##   Npat <- 100; ## Npat must be divisible by NpatEnterPerMonth

##   beta0 <- 1.5
##   beta1 <- 0.1
##   set.seed(iseed)

##   cat("\n\n beta0",beta0,"beta1",beta1,"\n\n")
##   ## The number of patients entering the study per month
##   DSMBVisitInterval <- 4 ## months

##   ##
##   ## I_{ij} follows a Markov Chain with two states 0/1 with transition probability p_{i0} = 0.9 and p_{i1}=0.1 for i = 0 or 1
##   ## That is I_{ij} is more likely to stay at state 0.
##   ##
##   ## mixture of two beta distributions
##   ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(1)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##   ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(1)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##   Ys <- Gs <- Ispat <- ID <- Is <- NULL


##     cat("\n\n Random effect changes at time 2, at time 3 (DSMB review) and at time 4 \n\n")

##     Ys <- Gs <-Is <- NULL
##     ni <- 6
##     timecov <- c(0,0,1:(ni-2))
##     NpatEach <- 5
##     NpatEach1 <- 2
##     Ispat <- matrix(c(
##                       rep(c(0,0,0,1,1,1),NpatEach),
##                       rep(c(1,1,1,0,0,0),NpatEach),
##                       rep(c(0,0,1,1,1,1),NpatEach1), ## NpatEach 
##                       rep(c(0,0,0,0,1,1),NpatEach1),
##                       rep(c(1,1,0,0,0,0),NpatEach1),
##                       rep(c(1,1,1,1,0,0),NpatEach1),
##                       rep(rep(1,ni),40), ## 40 
##                       rep(rep(0,ni),130) ## 130 reasonabl
##                       ),ncol=ni,byrow=TRUE)
##     Npat <-nrow(Ispat)
##     for (ipat in 1 : Npat)
##       {
##         gs <- ys <- NULL
##         ## random effect corresponding to Ispat[ipat,ivec]==0
##         g0 <- rbeta(1,shape1=aG2,shape2=rG2);
##         g1 <- rbeta(1,shape1=aG1,shape2=rG1);
##         for (ivec in 1 : ni)
##           {
##             if (Ispat[ipat,ivec] == 0)
##               {
##                 ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g0)
##                 gs <- c(gs,g0)
##               }else if (Ispat[ipat,ivec] == 1)
##                 {
##                   ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g1)
##                   gs <- c(gs,g1)
##                 }
##           }
##         Ys <- c(Ys,ys)
##         Gs <- c(Gs,gs)
##         Is <- c(Is,Ispat[ipat,])
##       }
##     d <- data.frame(Y=Ys,gs=Gs,Is=Is,
##                     ID=rep(1:Npat,each=ni),
##                     days=rep(1:ni,Npat),
##                     X1=rep(1,Npat*ni),
##                     X2=rep(timecov,Npat))

##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ 3 < d$days ] <- 1 

##     if (mod== 1)return(d)
##     else if (mod == 2){
##       ## return the parameter values:
##       c1 <- momentsBeta(aG1=aG1,rG1=rG1,beta=beta0)
##       c2 <- momentsBeta(aG1=aG2,rG1=rG2,beta=beta0)
##       REclass <- REclass(d)
##       return (list(c1=c1,c2=c2,beta0=beta0,beta1=beta1,aG1=aG1,rG1=rG1,aG2=aG2,rG2=rG2,REclass=REclass))
##     }

## }



## REclass <- function(d)
##   {
##     d <- data.frame(d)
##     Npat <- length(unique(d$ID))
##     REclass <- rep(NA,Npat)
##     for (ipat in 1 : Npat)
##       {
##         Is_pat <- d$Is[d$ID==ipat]
##         if (identical(Is_pat,c(1,1,1,1,1,1))){REclass[ipat] <- 1}
##         if (identical(Is_pat,c(0,0,0,0,0,0))){REclass[ipat] <- 2}
##         if (identical(Is_pat,c(1,1,0,0,0,0))){REclass[ipat] <- 3}
##         if (identical(Is_pat,c(0,0,0,0,1,1))){REclass[ipat] <- 4}
##         if (identical(Is_pat,c(1,1,1,1,0,0))){REclass[ipat] <- 5}
##         if (identical(Is_pat,c(0,0,1,1,1,1))){REclass[ipat] <- 6}
##         if (identical(Is_pat,c(0,0,0,1,1,1))){REclass[ipat] <- 7}
##         if (identical(Is_pat,c(1,1,1,0,0,0))){REclass[ipat] <- 8}
##       }
##     return(REclass)
##   }


## getSample2 <- function(iseed=1,
##                        mod=1, ## mod==2 returns the simulation parameter 
##                        aG1=10,rG1=10,aG2=20,rG2=1)
## {
##   ## The random effect changes its value over time
##   ## 20n % of patients switch their random effect values between pre and new scan period
##   ## By review 2, 80 patients are recruited to the study 
##   Npat <- 100; ## Npat must be divisible by NpatEnterPerMonth

##   beta0 <- 1.5
##   beta1 <- 0.1
##   set.seed(iseed)

##   cat("\n\n beta0",beta0,"beta1",beta1,"\n\n")
##   ## The number of patients entering the study per month
##   DSMBVisitInterval <- 4 ## months

##   ##
##   ## I_{ij} follows a Markov Chain with two states 0/1 with transition probability p_{i0} = 0.9 and p_{i1}=0.1 for i = 0 or 1
##   ## That is I_{ij} is more likely to stay at state 0.
##   ##
##   ## mixture of two beta distributions
##   ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(1)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##   ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(1)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##   Ys <- Gs <- Ispat <- ID <- Is <- NULL


##     cat("\n\n Random effect changes at time 2, at time 3 (DSMB review) and at time 4 \n\n")

##     Ys <- Gs <-Is <- NULL
##     ni <- 6
##     timecov <- c(0,0,1:(ni-2))
##     NpatEach <- 5
##     NpatEach1 <- 2
##     Ispat <- matrix(c(
##                       rep(c(0,0,0,1,1,1),NpatEach),
##                       rep(c(1,1,1,0,0,0),NpatEach),
##                       rep(c(0,0,1,1,1,1),NpatEach1), ## NpatEach 
##                       rep(c(0,0,0,0,1,1),NpatEach1),
##                       rep(c(1,1,0,0,0,0),NpatEach1),
##                       rep(c(1,1,1,1,0,0),NpatEach1),
##                       rep(rep(1,ni),40), ## 40 
##                       rep(rep(0,ni),130) ## 130 reasonabl
##                       ),ncol=ni,byrow=TRUE)
##     Npat <-nrow(Ispat)
##     for (ipat in 1 : Npat)
##       {
##         gs <- ys <- NULL
##         ## random effect corresponding to Ispat[ipat,ivec]==0
##         g0 <- rbeta(1,shape1=aG2,shape2=rG2);
##         g1 <- rbeta(1,shape1=aG1,shape2=rG1);
##         for (ivec in 1 : ni)
##           {
##             if (Ispat[ipat,ivec] == 0)
##               {
##                 ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g0)
##                 gs <- c(gs,g0)
##               }else if (Ispat[ipat,ivec] == 1)
##                 {
##                   ys[ivec] <- rnbinom(1,size = exp(beta0+beta1*timecov[ivec]),prob = g1)
##                   gs <- c(gs,g1)
##                 }
##           }
##         Ys <- c(Ys,ys)
##         Gs <- c(Gs,gs)
##         Is <- c(Is,Ispat[ipat,])
##       }
##     d <- data.frame(Y=Ys,gs=Gs,Is=Is,
##                     ID=rep(1:Npat,each=ni),
##                     days=rep(1:ni,Npat),
##                     X1=rep(1,Npat*ni),
##                     X2=rep(timecov,Npat))

##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ 3 < d$days ] <- 1 

##     if (mod== 1)return(d)
##     else if (mod == 2){
##       ## return the parameter values:
##       c1 <- momentsBeta(aG1=aG1,rG1=rG1,beta=beta0)
##       c2 <- momentsBeta(aG1=aG2,rG1=rG2,beta=beta0)
##       REclass <- REclass(d)
##       return (list(c1=c1,c2=c2,beta0=beta0,beta1=beta1,aG1=aG1,rG1=rG1,aG2=aG2,rG2=rG2,REclass=REclass))
##     }

## }




## getSampleS <- function(iseed=1,propTreat=2/3,RedFactor=0.5,detail=FALSE)
##   {
##     ## The subset of patients in treatment arm experience the treatment effect
##     ## Those with treatment effect decrease the expectation of CEL counts by RedFactor at every time point.
##     ## RedFactor is constant over time.

##     ## propTreat is the proportion of the patients in treatment arm with treatment effect
    
##     ## if full = TRUE then full dataset is returned and rev is ignored
##     ni <- 8
##     Npat <- 120;
##     i.trt <- c(
##                rep(1,round(Npat*propTreat*0.5)), ## 0.5 because half of the patients receive treatments
##                rep(0,Npat-round(Npat*propTreat*0.5))
##                )
##     pi <- 0.3
    
##     beta0 <- 1.5
##     beta1 <- -0.05
    
##     aG1 <- 10
##     rG1 <- 10
##     aG2 <- 20
##     rG2 <- 1

##     ## generate the initial random effect values of everyone
##     hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##     Npat_dist1 <- sum(hs==1)
##     Npat_dist2 <- sum(hs==2)
    
##     gsBASE <- rep(NA,Npat)
##     gsBASE[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
##     gsBASE[hs==2] <- rbeta(Npat_dist2,aG2,rG2);

##     ## generate samples for dist = "b","b3", "g", and "l" (except "b2")
##     Y <- NULL
##     timecov <- c(0,0,(1:(ni-2)))
##     for ( ipat in 1 : Npat)
##       {
##         ## the number of repeated measures are the same
##         ## we assume that the treatment effects occurs once after the screening and baseline
##         got <- rnbinom(ni,size = exp(beta0+beta1*timecov)*RedFactor^(i.trt[ipat]*(timecov>0)), prob = gsBASE[ipat])         
##         Y <- c(Y,got)
##       }
##     ## cat("beta0",beta0,"\n")
##     d <- data.frame(Y=Y,
##                     X1=rep(1,Npat*ni),
##                     X2=rep(timecov,Npat),
##                     labelnp=rep(c(rep(0,2),rep(1,ni-2)),Npat), ## the first two scans are pre-scans 
##                     ID=rep(1:Npat,each=ni),
##                     gsBASE = rep(gsBASE,each=ni),
##                     Treat = rep(1:0,each=(Npat/2)*ni), ## the first half of the patients receive treatment
##                     effectTreat = rep(i.trt,each=ni),
##                     hs = rep(hs,each=ni)
##                     )

##     if (detail)
##       {
##         E1G.1 <- momentsBeta(aG1,rG1,beta0)[3] ## E1G
##         E1G.2 <- momentsBeta(aG2,rG2,beta0)[3]
        
##         E1G <- E1G.1*pi+E1G.2*(1-pi)

##         E1s<-c(E1G,E1G.1,E1G.2)
##         names(E1s) <- c("Overall","H=1","H=2")
##         for (i in 1 : length(E1s))
##           {
##             noTreatEffect <- (E1s[i]-1)*exp(beta0+beta1*timecov)*RedFactor^(0*(timecov>0))
##             TreatEffect <- (E1s[i]-1)*exp(beta0+beta1*timecov)*RedFactor^(1*(timecov>0))
            
##             rr <- rbind(noTreatEffect,TreatEffect)
##             colnames(rr) <- timecov
##             cat("\n The expected CEL counts E(Yt)",names(E1s)[i],"\n")
##             print(rr)
##           }
##       }
    
##     return(d)
##   }

## ## E(K|N,D) digamma(D+N)-digamma(D) ~= D*log(1+N/D) D*log(1+N/D)
## K_N.D <- function(D,N=200) D*log(1+N/D)

## ## E(K|N) = E(E(K|N,D)) where D ~ gamma(shape=aD,scale1/rD)
## ## E(L|N) depends on aD and rD and N
## K_N <- function(aD,rD,N)
##   {
##     f <- function(D,N,rD,aD) log(1+N/D)*exp(-rD*D)*D^aD
##     k <- integrate(f=f,N=N,rD=rD,aD=aD,lower=0,upper=Inf)
##     return(k$value*rD^aD/gamma(aD))
##   }


## K_NisK <- function(aD,rD,N,K)
##   {
##     return(K_N(aD,rD,N)-K)
##   }





## nbinDPREchange <- function(Y,          ##   A vector of length sum ni, containing responses 
##                            X,          ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                            ID,         ##   A Vector of length sum ni, indicating patients
##                            B = 105000, ##     A scalar, the number of Gibbs iteration 
##                            burnin = 5000,  
##                            printFreq = B,
##                            M = NULL,
##                            labelnp, ## necessary if probIndex ==1
##                            epsilonM = 0.01,## nonpara
##                            para = list(mu_beta = NULL,Sigma_beta = NULL,a_D = 0.01, ib_D = 5,max_aG=30)
##                            )
##   {

##     ## This code generate samples from NBRE model with constant random effect ~ DP mixture of Beta
##     if (is.vector(X)) X <- matrix(X,ncol=1)
##     NtotAll <- length(Y)
##     if (nrow(X)!= NtotAll) stop("nrow(X) != length(Y)")
##     if (length(ID)!= NtotAll)  stop("length(ID)!= length(Y)")
##     if (length(labelnp)!= NtotAll)  stop ("labelnp!= length(Y)")

    
##     dims <- dim(X)
##     Ntot <- dims[1]
##     pCov <- dims[2]
##     ## cat("pCov",pCov)
##     if (is.null(para$mu_beta))
##       {
##         para$mu_beta <- rep(0,pCov)
##       }
##     mu_beta <- para$mu_beta

##     if (is.null(M)) M  <- round(1 + log(epsilonM)/log(para$ib_D/(1+para$ib_D)))
##     if (is.null(para$Sigma_beta)) {
##       para$Sigma_beta <-  diag(5,pCov)
##     }
##     Sigma_beta = para$Sigma_beta
    
##     evalue_sigma_beta <- eigen(Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
##     if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
##     Inv_sigma_beta <- c( solve(Sigma_beta) )

##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }

##     ## change the index of ID to numeric from 1 to # patients
##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }
    
##     mID <- ID-1
##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    
##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))


##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0

##     re <- .Call("gibbsREchange",
##                 as.numeric(Y),           ## REAL
##                 as.numeric(X),           ## REAL
##                 as.integer(mID),         ## INTEGER
##                 as.integer(B),           ## INTEGER
##                 as.integer(maxni),       ## INTEGER
##                 as.integer(Npat),        ## INTEGER
##                 as.numeric(labelnp),     ## REAL
##                 as.numeric(para$max_aG),
##                 as.numeric(para$mu_beta),     ## REAL
##                 as.numeric(evalue_sigma_beta),  ## REAL
##                 as.numeric(Inv_sigma_beta),  ## REAL
##                 as.numeric(para$a_D),
##                 as.numeric(para$ib_D),
##                 as.integer(M),
##                 as.integer(burnin),      ## INTEGER
##                 as.integer(printFreq),
##                 package = "lmeNBBayes"
##                 )


##     ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##     ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##     for ( i in 1 : 4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##     for ( i in 5 : 8 ) re[[i]] <- matrix(re[[i]],B,Npat,byrow=TRUE)
##     re[[9]] <- matrix(re[[9]],B,pCov,byrow=TRUE)
##     names(re) <- c("aGs","rGs","vs","weightH1",
##                    "h1s","h2s","g1s","g2s",
##                    "beta",
##                    "logL","D",
##                    "AR","prp")

##     re$para <- para
##     re$para$M <- M
##     names(re$AR) <-names(re$prp) <- c("aG", "rG","beta","D")
    
    
##     return (re)
##   }


## DPfit <- function(Y,          ##   A vector of length sum ni, containing responses 
##                   X,          ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                   ID,         ##   A Vector of length sum ni, indicating patients
##                   B = 105000, ##     A scalar, the number of Gibbs iteration 
##                   burnin = 5000,  
##                   printFreq = B,
##                   model = "nonpara-unif",
##                   M = NULL,
##                   epsilonM = 0.01,## nonpara
##                   prob=0.9,        ## nonpara
##                   labelnp=NULL,
##                   para = list(mu_beta = NULL,Sigma_beta = NULL,a_D = 1, r_D = 1, a_qG = 1, r_qG = 1,
##                     mu_aG = 0,sd_aG = 2, mu_rG = 0,sd_rG = 2, Nclust=NULL),
##                   initBeta=NULL
##                   )
##   {

##     Ntot <- length(Y)


##     ## === prior input check =====## 
##     if (is.vector(X)) X <- matrix(X,ncol=1)
    


##     if (nrow(X)!= Ntot) stop ("nrow(X) != length(Y)")
##     if (length(ID)!= Ntot)  stop ("length(ID)!= length(Y)")
    
##     if (is.null(para$a_qG)){
##       para$a_qG <- 1
##       cat("\n a_qG needs to be specified!! set to 1")
##     }
##     a_qG = para$a_qG
    
##     if (is.null(para$r_qG) ) {
##       para$r_qG <- 1
##       cat("\n r_qG needs to be specified!! set to 1")
##     }
##     r_qG = para$r_qG

##     if (is.vector(X)) X <- matrix(X,ncol=1)
##     pCov <- ncol(X)

##     if (is.null(para$mu_beta))
##       {
##         para$mu_beta <- rep(0,pCov)
##       }
    
##     mu_beta = para$mu_beta
    
##     if (is.null(para$Sigma_beta)) {
##       para$Sigma_beta <-  diag(5,pCov)
##     }
##     Sigma_beta = para$Sigma_beta

##     if (is.null(labelnp)) labelnp <- rep(0,length(Y))
    
    
    

##     KK <- para$Nclust
##     if (!is.null(KK)){
##       if (!is.null(para$a_D)) stop("Both KK and a_D are selected! must be one of them...")
##       ## the total number of patients observed by irev^th review
##       NtotID <- length(unique(ID))
##       ## the number of patients whose new scans are observed at the irev^th review
##       NnewID <- length(unique(ID[labelnp==1]))
##       ## parameters
##       para$a_D <- uniroot(f=K_NisK,interval=c(0.003,50),
##                           rD=para$r_D,N=(NtotID+NnewID)*2,K=KK)$root
      
##     }
##     a_D <- para$a_D
    
##     if (is.null(M)) M  <- round(1 + log(0.01)/log(5/(1+5)))

##     ## change the index of ID to numeric from 1 to # patients
##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }
    
##     p <- pCov
##     if( is.vector(X)) p <- 1 
##     X <- c(X) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }
##     mID <- ID-1

##     ## The patients with labelnp = 1 (no old scans) for all repeated measures
##     ## are treated as lack of new scans and used to estimate beta only.
##     ## skip the computation of H2
##     ## All patients have old scans
    
##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##     patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##     for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0
    
##     patwoNS <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
##     if (length(patwoNS)==0) patwoNS <- -999;
##     patwoNS <- patwoNS - 1
    
##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))

    
##     evalue_sigma_beta <- eigen(Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
##     if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
##     Inv_sigma_beta <- c( solve(Sigma_beta) )
    
##     if (model=="para-constantRE"| model==1)
##       {
        
        
##         if (is.null(para$mu_aG )){
##           para$mu_aG <- 0.5
##           cat("\n mu_aG needs to be specified!! set to 0.5")
##         }
##         mu_aG = para$mu_aG
        
##         if (is.null(para$mu_rG)){
##           para$mu_rG <- 0.5
##           cat("\n mu_rG needs to be specified!! set to 0.5")
##         }
##         mu_rG = para$mu_rG
        
##         if (is.null(para$sd_aG) ){
##           para$sd_aG <- 2
##           cat("\n mu_rG needs to be specified!! set to 2")
##         }
##         sd_aG = para$sd_aG
        
##         if (is.null(para$sd_rG)) {
##           para$sd_rG <- 2
##           cat("\n mu_rG needs to be specified!! set to 2")
##         }
##         sd_rG = para$sd_rG
        
##         ##parametric model: gi is constant over time 
##         labelnp <- rep(0,length(Y))    
##         re <- .Call("Beta1",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(mu_aG),
##                     as.numeric(sd_aG),
##                     as.numeric(mu_rG),
##                     as.numeric(sd_rG),
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(evalue_sigma_beta),  ## REAL
##                     as.numeric(Inv_sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     package = "lmeNBBayes"
##                     )
        
##         re[[3]] <- matrix(re[[3]],B,N,byrow=TRUE)
##         re[[4]] <- matrix(re[[4]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","gPre","beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      Sigma_beta = Sigma_beta,
##                      B=B,
##                      model="para-constantRE",
##                      burnin = burnin,
##                      mu_aG = mu_aG,
##                      sd_aG = sd_aG,
##                      mu_rG = mu_rG,
##                      sd_rG = sd_rG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", "beta")
##         return(re)
        
##       }else if (model=="para"){
##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ beta(aG,rG)
##         ## // beta ~ rnorm(mu_beta,sigma_beta)
##         ## // aG, rG ~ lognorm(mu_aG,sd_aG),lognorm(mu_rG,sd_rG)
        
        
##         if (is.null(para$mu_aG )){
##           para$mu_aG <- 0.5
##           cat("\n mu_aG needs to be specified!! set to 0.5")
##         }
##         mu_aG = para$mu_aG
        
##         if (is.null(para$mu_rG)){
##           para$mu_rG <- 0.5
##           cat("\n mu_rG needs to be specified!! set to 0.5")
##         }
##         mu_rG = para$mu_rG
        
##         if (is.null(para$sd_aG) ){
##           para$sd_aG <- 2
##           cat("\n mu_rG needs to be specified!! set to 2")
##         }
##         sd_aG = para$sd_aG
        
##         if (is.null(para$sd_rG) ) {
##           para$sd_rG <- 2
##           cat("\n mu_rG needs to be specified!! set to 2")
##         }
##         sd_rG = para$sd_rG
        
##         re <- .Call("Beta14",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),
##                     as.numeric(r_qG),
##                     as.numeric(mu_aG),
##                     as.numeric(sd_aG),
##                     as.numeric(mu_rG),
##                     as.numeric(sd_rG),
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(evalue_sigma_beta),  ## REAL
##                     as.numeric(Inv_sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 4 : 7 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[8]] <- matrix(re[[8]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","qG",
##                        "gPre","g2s","gNew","js",
##                        "beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      Sigma_beta = Sigma_beta,
##                      B=B,
##                      model="para",
##                      burnin = burnin,
##                      mu_aG = mu_aG,
##                      sd_aG = sd_aG,
##                      mu_rG = mu_rG,
##                      sd_rG = sd_rG,
##                      a_qG=a_qG,
##                      r_qG=r_qG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", "beta")
##         return(re)
        
        
##       }else if (model=="para-unif"){
##         if (is.null(para$max_aG) )  para$max_aG <- 30
##         max_aG <- para$max_aG
##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ beta(aG,rG)
##         ## // beta ~ rnorm(mu_beta,sigma_beta)
##         ## // aG, rG ~ unif(0,max_aG)

##         re <- .Call("Beta14Unif",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),
##                     as.numeric(r_qG),
##                     as.numeric(max_aG),
##                     ## as.numeric(mu_aG),
##                     ## as.numeric(sd_aG),
##                     ## as.numeric(mu_rG),
##                     ## as.numeric(sd_rG),
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(evalue_sigma_beta),  ## REAL
##                     as.numeric(Inv_sigma_beta),  ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     package = "lmeNBBayes"
##                     )
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 5 : 8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         re[[9]] <- matrix(re[[9]],B,p,byrow=TRUE)
##         names(re) <- c("aG","rG","qG","logL",
##                        "gPre","g2s","gNew","js",
##                        "beta","AR","prp")
##         para <- list(mu_beta = mu_beta,
##                      Sigma_beta = Sigma_beta,
##                      B=B,
##                      model="para",
##                      burnin = burnin,
##                      max_aG = max_aG,
##                      ## mu_aG = mu_aG,
##                      ## sd_aG = sd_aG,
##                      ## mu_rG = mu_rG,
##                      ## sd_rG = sd_rG,
##                      a_qG=a_qG,
##                      r_qG=r_qG
##                      )
##         re$para <- para
##         names(re$AR) <-c("aG", "rG", "beta")
##         return(re)
##       }else if (model== "nonpara" | model==8){

##         ## nonparametric model
##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ sum_{h=1}^M pi_h beta(aG_h,rG_h)
##         ## // beta ~ mvrnorm(mu_beta,sigma_beta)
##         ## // aG ~lognorm(mu_aG,sd_aG)
##         ## // rG ~ lognorm(mu_rG,sd_rG)
##         ##  nonparametric model,
        
        
##         if (is.null(para$mu_aG )){
##           para$mu_aG <- 0.5
##           cat("\n mu_aG needs to be specified!! set to 0.5")
##         }
##         mu_aG = para$mu_aG
        
##         if (is.null(para$mu_rG)){
##           para$mu_rG <- 0.5
##           cat("\n mu_rG needs to be specified!! set to 0.5")
##         }
##         mu_rG = para$mu_rG
        
##         if (is.null(para$sd_aG) ){
##           para$sd_aG <- 2
##           cat("\n mu_rG needs to be specified!! set to 2")
##         }
##         sd_aG = para$sd_aG
        
##         if (is.null(para$sd_rG) ) {
##           para$sd_rG <- 2
##           cat("\n mu_rG needs to be specified!! set to 2")
##         }
##         sd_rG = para$sd_rG
##         if (is.null(para$r_D) ){
##           cat("\n r_D needs to be specified!! set to 1")
##           para$r_D <- 1
##         }
##         r_D = para$r_D
##         re <- .Call("Beta24",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.integer(M),           ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),        ## REAL
##                     as.numeric(r_qG),        ## REAL
##                     as.numeric(mu_aG),       ## REAL
##                     as.numeric(sd_aG),       ## REAL
##                     as.numeric(mu_rG),       ## REAL
##                     as.numeric(sd_rG),       ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(evalue_sigma_beta),  ## REAL
##                     as.numeric(Inv_sigma_beta),  ## REAL
##                     as.numeric(a_D),         ## REAL
##                     as.numeric(r_D),         ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     package = "lmeNBBayes"
##                     )
##         for ( i in 3 : 8 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         for ( i in 9 : 12 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         re[[13]] <- matrix(re[[13]],B,p,byrow=TRUE)
##         names(re) <- c("qG","D",
##                        "gPre","g2s","gNew","js","h1s","h2s",
##                        "weightH1","vs","aGs","rGs",
##                        "beta","AR","prp")
##         re$para <- list(
##                         burnin = burnin,
##                         M=M,
##                         B=B,
##                         model=model,
##                         Npat=N,
##                         Ntot=length(Y),
##                         para
##                         )
##         names(re$AR) <- names(re$prp) <- c(paste("aG",1:M,sep=""),
##                                            paste("rG",1:M,sep=""),
##                                            paste("beta",sep="")
##                                            )
##         re$h1s <- re$h1s + 1
##         re$h2s <- re$h2s + 1
##         re$MeanWH <- colMeans(re$weightH1);
##         names(re$MeanWH) <- paste("cluster",1:M)
##         return(re)
        
##       }else if (model=="nonpara-unif"){
        
##         ## nonparametric model D: is fixed
##         ## // Y_ij | Gij = gij ~ NB(size=exp(X_{ij}^T beta),prob=gij)
##         ## // gij = g1 if j is in pre-scan 
##         ## //     = g_new if j is in old-scan 
##         ## // g_new = Ji * g1 + (1-Ji) * g2 
##         ## // Ji ~ ber(qG)
##         ## // qG ~ beta(a_qG,r_qG)
##         ## // g1, g2 ~ sum_{h=1}^M pi_h beta(aG_h,rG_h)
##         ## // beta ~ mvrnorm(mu_beta,sigma_beta)
##         ## // aG,rG ~ unif(0.0001,max_aG)
##         ##  nonparametric model,
##         if (is.null(para$max_aG) )  para$max_aG <- 30
##         max_aG <- para$max_aG

##         if (is.null(para$max_D))  para$max_D <- 5
##         maxD <- para$max_D 

##         if (is.null(initBeta)) initBeta <- rep(0, pCov)
##         if (is.null(M)) M  <- 10 ## Need to fix this 
        
##         re <- .Call("Beta24Unif",
##                     as.numeric(Y),           ## REAL
##                     as.numeric(X),           ## REAL
##                     as.integer(mID),         ## INTEGER
##                     as.integer(B),           ## INTEGER
##                     as.integer(maxni),       ## INTEGER
##                     as.integer(Npat),        ## INTEGER
##                     as.integer(M),           ## INTEGER
##                     as.numeric(labelnp),     ## REAL
##                     as.numeric(a_qG),        ## REAL
##                     as.numeric(r_qG),        ## REAL
##                     as.numeric(max_aG),
##                     ## as.numeric(mu_aG),    ## REAL
##                     ## as.numeric(sd_aG),    ## REAL
##                     ## as.numeric(mu_rG),    ## REAL
##                     ## as.numeric(sd_rG),    ## REAL
##                     as.numeric(mu_beta),     ## REAL
##                     as.numeric(evalue_sigma_beta),  ## REAL
##                     as.numeric(Inv_sigma_beta),  ## REAL
##                     as.numeric(maxD),         ## REAL
##                     ## as.numeric(r_D),         ## REAL
##                     as.integer(burnin),      ## INTEGER
##                     as.integer(printFreq),
##                     as.integer(patwoNS),
##                     as.numeric(initBeta),
##                     package = "lmeNBBayes"
##                     )
##         for ( i in 4 : 9 ) re[[i]] <- matrix(re[[i]],B,N,byrow=TRUE)
##         for ( i in 10 : 13 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         re[[14]] <- matrix(re[[14]],B,p,byrow=TRUE)
##         names(re) <- c("qG","D","logL",
##                        "gPre","g2s","gNew","js","h1s","h2s",
##                        "weightH1","vs","aGs","rGs",
##                        "beta","AR","prp")
        
##         para$mu_aG <- para$mu_rG <- para$sd_aG <- para$sd_rG <- para$a_D <- para$r_D <- NULL
        
##         re$para <- list(
##                         burnin = burnin,
##                         M=M,
##                         B=B,
##                         model=model,
##                         Npat=N,
##                         Ntot=length(Y),
##                         a_qG=a_qG,
##                         r_qG=r_qG,
##                         max_aG=max_aG,
##                         mu_beta=mu_beta,
##                         Sigma_beta=Sigma_beta,
##                         maxD=maxD
##                         )
##         names(re$AR) <- names(re$prp) <- c(paste("aG",1:M,sep=""),
##                                            paste("rG",1:M,sep=""),
##                                            paste("beta",sep=""),
##                                            "D"
##                                            )
##         re$h1s <- re$h1s + 1
##         re$h2s <- re$h2s + 1
##         re$MeanWH <- colMeans(re$weightH1);
##         names(re$MeanWH) <- paste("cluster",1:M)
##         return(re)

##       }
##   }



## DevRatio <- function(gPre,gNew) ## gPre, gNew must be a matrix of length(useSample) by N
##   {
##     ## The deviation of CEL counts of single patient i from overall cohort at week j could be measured by:
##     ## R_ij=E(Y_ij|G_ij)/E(Y_ij)= r_ij (1-g_ij)/gij * 1/(rij*(mu_{1/G}-1)) = (1-g_ij)/(gij*(mu_{1/G}-1))
##     ## The change in the deviation R_ij at two time point could be measured by
##     ## R_ij/R_ij' = (1-g_ij)/(gij*(mu_{1/G}-1))*(gij'*(mu_{1/G}-1))/(1-g_ij) = g_ij'*(1-g_ij)/(g_ij*(1-g_ij'))
##     ## Hence the deviation ratio between the pre scan period and new scan period is:
##     ## R_inew/R_ipre
##     ## How much patient i increases the deviation from the overall trend between the pre-scan period and new-scan period

##     ## the elemnt-wise operations
##     m1G_1 <- mean(c(1/gPre,1/gNew)) - 1
    
##     EYGpre <- 1/gPre - 1
##     EYGnew <- 1/gNew - 1


##     CI95_EYGpre <- apply(EYGpre,2,quantile,prob=c(0.5,0.025,0.975)) / m1G_1
##     CI95_EYGnew <- apply(EYGnew,2,quantile,prob=c(0.5,0.025,0.975)) / m1G_1

##     CI95_EYGpre[1,] <- colMeans(EYGpre)
##     CI95_EYGnew[1,] <- colMeans(EYGnew)
    
##     ##DevMat <- EYGnew/EYGpre ## B-Bburn by N
##     ##CI95 <- apply(DevMat,2,quantile,prob=c(0.5,0.025,0.975))

##     NewPre <-  CI95_EYGnew[1,]/CI95_EYGpre[1,]
    
##     ##CI95[1,] <- colMeans(DevMat)
    
##     ##rownames(CI95)[1] <- "mean"
##     rownames(CI95_EYGpre)[1] <- "mean"
##     rownames(CI95_EYGnew)[1] <- "mean"
##     if (!is.null(rownames(gPre)))
##       {
##         names(CI95) <- colnames(gPre)
##         colnames(CI95_EYGpre) <- colnames(gPre)
##         colnames(CI95_EYGnew) <- colnames(gPre)
##       }
##     re <- list(Ratio=NewPre,
##                Rpre=CI95_EYGpre,
##                Rnew=CI95_EYGnew)
    
##     return (re)
##   }


## test <- function(x,mu,Sigma)
##   {
##     InvSigma <- solve(Sigma)
##     evalue <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
##     .Call("tempFun",
##           as.numeric(x),
##           as.numeric(mu),
##           as.integer(length(mu)),
##           as.numeric(evalue),
##           as.numeric(c(InvSigma)),
##           package ="lmeNBBayes")
##   }
## infoPriorTrt <- function(info.mean,info.var,p0s,add=0)
##   {
##     ## develop Informative prior under the random sample assumption

##     ## p0s = c(Pr(plcb),Pr(trt))
##     ## Save parameters
##     input <- list(info.mean=info.mean,info.var=info.var,p0s=p0s,add=add)
##     ## The number of 3-month interval during the followup
##     Nint <- 3
##     ## info.mean is a vector of length 1 + Nint*2 containing
##     ## (Intercept) trt010:timeInt1 trt011:timeInt1 trt010:timeInt2 trt011:timeInt2 trt010:timeInt3 trt011:timeInt3
##     ##  alpha       beta_plcb,1     beta_trt,1     beta_plcb,2      beta_trt,2 ,....
##     ## which column of info.var corresponds to coefficients of placebo patients
##     which.plcb <- c(2,4,6)
##     which.trt <- which.plcb + 1
##     ## the first entry of info.mean is intercept coefficient (alpha)
##     m.alpha <- info.mean[1]
##     m.betas.plcb <- info.mean[which.plcb]
##     m.betas.trt <- info.mean[which.trt]
##     m.beta01 <- rep(NA,Nint+1)

##     ## The mean of the intercept term does not change
##     m.tild.beta01[1] <- info.mean[1]
##     s.tild.beta01 <- matrix(NA,Nint+1,Nint+1) ## + 1 for alpha (intercept)
##     ## The variance of the intercept term does not change
##     s.tild.beta01[1,1] <- info.var[1,1]
    
##     for (itime1 in 1 : Nint)
##       {
##         m.bet0 <- m.betas.plcb[itime1]
##         m.bet1 <- m.betas.trt[itime1]
##         ## Extract the index to pick
##         which.pt <- c(which.plcb[itime1],which.trt[itime1])
##         Sigma <- info.var[which.pt, which.pt]
##         mt <- moments.tildeBeta(mu=c(m.bet0,m.bet1),Sigma=Sigma,p0=p0s)
##         m.tild.beta01[itime1+1] <- mt[1] ## mean tilde.beta[itime1]
##         s.tild.beta01[itime1+1,itime1+1] <- mt[3] ## var(tilde.beta)
##         ## Compute the covariance between tilde.beta and alpha
##         mt <- moments.tildeBeta(mu=c(m.bet0,m.bet1,m.alpha),
##                                 Sigma=info.var[c(which.pt,1), c(which.pt,1)], 
##                                 p0=p0s)
##         s.tild.beta01[1,itime1+1] <- s.tild.beta01[itime1+1,1] <- mt - m.tild.beta01[1]*m.tild.beta01[itime1 + 1]
##         if (itime1 > 1){
##           for (itime2 in 1 : (itime1-1)){
##             m.bet0.2 <- m.betas.plcb[itime2]
##             m.bet1.2 <- m.betas.trt[itime2]
##             which.pt2 <- c(which.plcb[itime2],which.trt[itime2])
##             Sigma <- info.var[c(which.pt,which.pt2), c(which.pt,which.pt2)]
##             mt <- moments.tildeBeta(mu=c(m.bet0,m.bet1,m.bet0.2,m.bet1.2),Sigma=Sigma,p0=p0s) 
##             s.tild.beta01[itime1+1,itime2+1] <- s.tild.beta01[itime2+1,itime1+1] <- mt - m.tild.beta01[itime1]*m.tild.beta01[itime2]
##           }
##         }
##       }
##     print("info.var before any truncation")
##     print(s.tild.beta01)
##     re <- adjustPosDef(mat=s.tild.beta01,zero=add)
##     info.var <- re$adjust.mat;
##     colnames(info.var) <-  rownames(info.var) <- colnames(s.tild.beta01) <- rownames( s.tild.beta01) <- names(m.tild.beta01) <- c("alpha (intercept)",paste("tildeBeta",1:Nint,sep=""))
    
##     print(paste("eigenvalue of the informative prior variance"))
##     print(eigen(info.var))
##     print(paste("correlation matrix of the informative prior variance"))
##     print(cov2cor(info.var))
    
##     return(list(info.mean.trt=m.tild.beta01,info.var.trt=info.var,
##                 para=list(input,orig.var=s.tild.beta01),note=re$note))
##   }


## moments.tildeBeta <- function(mu,Sigma,p0s){
##   ## mu = c(mu_{beta_{plcb,t}} mu_{beta_{trt,t}})
##   ## p0s = c( Pr(plcb), Pr(trt) ) 2 by 1
##   ## Sigma = Var(c(beta_{plcb,t},beta_{trt,t})) 2 by 2
##   Ndim <- length(mu)
 
##   if (length(p0s)==1) p0s[2] <- 1 - p0s[1]
##   else if( sum(p0s)!=1 ) stop("Pr(placebo) + Pr(TRT) must be 1")
##   if (length(p0s)!=2) stop("The length of p0s must be 2, each containing the pr(placebo), pr(trt)")
  
##   if (dim(Sigma)[1]!=Ndim || dim(Sigma)[2]!=Ndim)
##     stop("The dimension of Sigma does not much with the length of mu")
  
##   evalue_sigma <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
##   if (min(evalue_sigma) <= 0) stop("Sigma_beta must be positive definite!")
##   Inv_sigma <- c( solve(Sigma) )
  
##   if (Ndim==2){
##     ## Compute E(beta.tilde) = int int log(exp(a0)*pi0 + exp(a1)*pi1)*f(a0,a1) da0 da1 where f(x,y) = joint pdf of (mu.beta_{0,t},mu.beta_{1,t})
##     re <- .Call("EYP1or2",
##                 as.numeric(mu), as.numeric(evalue_sigma),
##                 as.numeric(Inv_sigma),
##                 as.numeric(p0s),
##                 package ="lmeNBBayes")[[1]]
##     re[3] <- re[2] - re[1]^2 ## Variance
##     names(re) <-c("E( log[exp(X1)P1+exp(X2)P2] )","E( log[exp(X1)P1+exp(X2)P2]^2 )","Var(log[exp(X1)P1+exp(X2)P2])")
    
##   }else if (Ndim %in% 3:4){
##     ## If ndim==4 compute E{ log(exp(X1)Pr(A1)+exp(X2)Pr(A2))*log(exp(X1')Pr(A1)+exp(X2')Pr(A2)) }
##     ## mu[1] = placebo beta at time interval k
##     ## mu[2] = trt beta at time interval k
##     ## mu[3] = placebo beta at time interval k'
##     ## mu[4] = trt beta at time interval k'

##     ## If ndim==3 compute E{ X1*log(exp(X1')Pr(A1)+exp(X2')Pr(A2)) }
##     ## mu[1:2] = c(E(beta.t.0), E(beta.t.1))
##     ## mu[3] = E(alpha),
##     ## Sigma[1:2,1:2]= Var(c(beta.t.0,beta.t.1)), Sigma[3,3] = Var(alpha) 
  
##     re <- .Call("EYPEXP",
##                 as.numeric(mu), as.numeric(evalue_sigma),
##                 as.numeric(Inv_sigma),
##                 as.numeric(p0s),
##                 as.integer(Ndim),
##                 package ="lmeNBBayes")[[1]]
##    if (Ndim == 4) names(re) <- "E( log[exp(X1)P1+exp(X2)P2]*log[exp(X1')P1+exp(X2')P2] )"
##    else if (Ndim == 3) names(re) <- "E( alpha*log[exp(X1')P1+exp(X2')P2] )"
##   }
##   return(re)
## }



## momentsBeta <- function(aG1,rG1,beta)
##   {
##     ## Y ~ NB(size=rij,prob=G)
##     ## G ~ Beta(aG1,rG1)
##     ## E(1/G) = ...some calculations... = (aG+rG-1)/(aG-1)
##     ## E(Y)=exp(beta0)*(mu_{1/G}+1)=exp(1)*((aG1+rG1-1)/(aG1-1)+1)
##     ## Var(1/G) = E(1/G^2) - E(1/G)^2 = (aG1+rG1-1)/(aG1-1)*( (aG1+rG1-2)/(aG1-2) - (aG1+rG1-1)/(aG1-1) ) = 0.06559506
##     ## Var(Y) = rij*Var(1/G)*(rij+1) + rij*E(1/G)*(E(1/G)-1)
##     ## calculate E(Y), Var(Y), E(1/G), and E(1/G^2) at initial time
##     rij <- exp(beta) ## value of size parmeter at time zero
##     E1G <- (aG1+rG1-1)/(aG1-1)
##     V1G <- (aG1+rG1-1)/(aG1-1)*( (aG1+rG1-2)/(aG1-2) - (aG1+rG1-1)/(aG1-1) )
##     EY <- rij*(E1G-1)
##     VY <- rij*V1G*(rij+1)+rij*E1G*(E1G-1)
##     return (c(EY=EY,SDY=sqrt(VY),
##               E1G=E1G,V1G=V1G))
##   }


## momentsGL <- function(EG,VG,beta)
##   {
##     ## Y ~ NB(size=rij,prob=1/(G+1))
##     ## G ~ dist(mean=EG,var=VG)
##     rij <- exp(beta)
##     ## E(Y|G) = rij*G
##     ## Var(Y|G) = rij*G*(G+1)
##     ## E(Y)=E(E(Y|G))=rij*E(G)
##     EY <- rij*EG
##     ## Var(Y) = Var(E(Y|G)) + E(Var(Y|G))
##     ##        = Var(rij*G) + E(G*(G+1)*rij)
##     ##        = rij^2*Var(G) + rij*(E(G^2)+E(G))
##     ##        = rij^2*Var(G) + rij*(Var(G)+E(G)^2+E(G))
##     VY <- rij^2*VG + rij*(VG+EG^2+EG)
##     return (c(EY=EY,SDY=sqrt(VY),
##               EG=EG,VG=VG))
##   }

## paralnorm <- function(E,V=15^2){
##   ## E: expectation of log-normally distributed r.v.
##   ## V: variance of log-normally distributed r.v

##   sd <- sqrt(log(V/E+1))
##   mu <- log(E)-sd/2
##   return(list(mu=mu,sd=sd))
## }
## F_inv_gam <- function(p,sp1,sc1,sp2,sc2,pi)
##   {
##     G = function (t)
##       pi*pgamma(t,shape=sp1,scale=sc1)+(1-pi)*pgamma(t,shape=sp2,scale=sc2) - p
##                 return(uniroot(G,c(0,100))$root)
##   }

## getSample <- function(
##                       iseed = 911,
##                       rev = 4,
##                       dist = "b",
##                       mod = 0,
##                       probs = seq(0,0.99,0.01),
##                       ts = seq(0.001,0.99,0.001),
##                       full = FALSE
##                       )
##   {
##     ## mod = 0: generate sample
##     ## mod = 1: quantiles of the true populations
##     ## mod = 2: densities of the true populations
##     ## mod = 3: parameters of the simulation model
##     ## mod = 4: parameters for uninformative prior

##     ## dist = "b","b2","g" 

##     ## if full = TRUE then full dataset is returned and rev is ignored
##     ni <- 11
##     Npat <- 180; ## upto review 4, Npat=160 in total this number must be divisible by NpatEnterPerMonth
##     NpatEnterPerMonth <- 15
##     DSMBVisitInterval <- 4 ## months
    
##     varInfoP <- c(0.1,0.01)
##     ## === necessary for mod= 0, 3 and mod = 4
##     ## d contains a full dataset
##     days <- NULL
##     for (ipatGroup in 1 : (Npat/NpatEnterPerMonth))
##       {
##         ScandaysForSingleGroup <- ipatGroup:(ipatGroup+ni-1)
##         days <- c(days,rep(ScandaysForSingleGroup,NpatEnterPerMonth))
##       }

##     set.seed(iseed)
##     if (dist=="g") ## need to went through by all mod=0,...,4
##       {
##         ## r.e. mixture of two gamma distributions
##         pi <- 0.7
##         sp1 <- 0.3; sc1 <- 0.5
##         sp2 <- 4; sc2 <- 0.8
##         beta0 <- 0.5
##         beta1 <- 0
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)

##         gsBASE <- rep(NA,Npat)
##         gsBASE[hs==1] <- rgamma(Npat_dist1,shape=sp1,scale=sc1)
##         gsBASE[hs==2] <- rgamma(Npat_dist2,shape=sp2,scale=sc2)
##         gsBASE <- 1/(1+gsBASE)
        
##         if (mod==1)
##           {
##             lq <- length(probs);
##             quan <- rep(NA,lq);
##             for (i in 1 : lq) quan[i] <- F_inv_gam(p=probs[i],sp1=sp1,sc1=sc1,sp2=sp2,sc2=sc2,pi=pi)
##             quan = sort(1/(1+quan))
            
##             return (rbind(probs=probs,quantile=quan))

##           }else if (mod == 2){
            
##             ts.trans <-1/ts-1
##             return ((pi*dgamma(ts.trans,shape=sp1,scale=sc1)+(1-pi)*dgamma(ts.trans,shape=sp2,scale=sc2) )*(1/ts^2))
            
##           }else if (mod == 3){
            
##             c1 <- momentsGL(EG=sp1*sc1,VG=sp1*sc1^2,beta=beta0) ## expectation and variance of Y
##             c2 <- momentsGL(EG=sp2*sc2,VG=sp2*sc2^2,beta=beta0)
            
##             outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
##                                  Sigma_beta=diag(varInfoP),a_qG=100,r_qG=1),
##                                beta0=beta0,beta1=beta1,K=2,
##                                scales=c(sc1=sc1,sc2=sc2),
##                                shapes=c(sp1=sp1,sp2=sp2),
##                                c1=c1,c2=c2
##                                )
##           }
##       }else if (dist == "b"){
       
##         beta0 <- 1
##         beta1 <- -0.05
##         aG1 <- 3
##         rG1 <- 0.8
        
##         gsBASE <- rbeta(Npat,aG1,rG1)
##         hs <- rep(1,Npat)
##         if (mod==1){
##           return( 
##                  rbind(probs=probs,
##                        quantile=qbeta(probs,shape1=aG1,shape2=rG1)
##                        )
##                  )
##         }else if (mod==2){
          
##           return (dbeta(ts,shape1=aG1,shape2=rG1))
          
##         }else if (mod == 3){
          
##           c1 <- momentsBeta(aG1,rG1,beta0)
##           outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
##                                Sigma_beta=diag(varInfoP),
##                                max_aG=30
##                                ),##a_qG=10,r_qG=1),
##                              beta0=beta0,beta1=beta1,K=1,c1=c1, aGs=c(aG1=aG1),
##                              rGs=c(rG1=rG1))
##         }

##       }else if (dist == "b2")
##         {
##         ## mixture of two beta distributions
##         ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##         ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##         pi <- 0.3
        
##         beta0 <- 1.5
##         beta1 <- -0.05
        
##         aG1 <- 10
##         rG1 <- 10
##         aG2 <- 20
##         rG2 <- 1

##         ## generate the initial random effect values of everyone
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gsBASE <- rep(NA,Npat)
##         gsBASE[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
##         gsBASE[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
##         if (mod==1){
##           F_inv_beta2 <- function(p,aG1,rG1,aG2,rG2,pi)
##             {
##               G = function (t) pi*pbeta(t,shape1=aG1,shape2=rG1)+(1-pi)*pbeta(t,shape1=aG2,shape2=rG2) - p
##               return(uniroot(G,c(0,1000))$root)
##             }
##           lq <- length(probs);
##           quant <- rep(NA,lq);
##           for (i in 1 : lq) quant[i] <- F_inv_beta2(p=probs[i],aG1,rG1,aG2,rG2,pi=pi)
##           return (rbind(probs=probs,quantile=quant))
##         }
##         else if (mod==2){

##           return ((pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2)) )
##         }else if (mod == 3){

##           c1 <- momentsBeta(aG1,rG1,beta0)
##           c2 <- momentsBeta(aG2,rG2,beta0)
##           outputMod3 <-  list(infoPara = list(mu_beta=c(beta0,beta1),
##                                 Sigma_beta=diag(varInfoP),max_aG=30),
##                               beta0=beta0,beta1=beta1,K=2,c1=c1,c2=c2,
##                               aGs=c(aG1=aG1,aG2=aG2),
##                               rGs=c(rG1=rG1,rG2=rG2),
##                               pi=pi
##                               )
##         }
##       }
##     if ( dist == "YZ")
##       {
##         pi <- 0.85

##         alpha <- exp(-0.5)
##         ## logalpha <- -0.5
##         beta0 <- 0.905 ##0.405+0.5 beta - logalpha
##         beta1 <- 0
##         ## a bimodal distribution with 85 % of Gi from a gamma distribution with mean 0.647 and variance 2.374
##         scale <- 2.374/0.647*alpha 
##         shape <- 0.647^2/2.374
##         mu <- 3*alpha
##         sd <- sqrt(0.25)*alpha
       
##         ## generate the initial random effect values of everyone
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gsBASE <- rep(NA,Npat)
##         gsBASE[hs==1] <- rgamma(Npat_dist1,shape=shape,scale=scale);
##         gsBASE[hs==2] <- rnorm(Npat_dist2,mean=mu,sd=sd);
##         gsBASE[gsBASE < 0 ] <- 0
##         gsBASE <- 1/(1+gsBASE)
        
##         if (mod==1)
##           {
##             return (NULL)

##           }else if (mod == 2){
            
##             ts.trans <-1/ts-1
##             return ((pi*dgamma(ts.trans,shape=shape,scale=scale)+(1-pi)*dnorm(ts.trans,mean=mu,sd=sd) )*(1/ts^2))
            
##           }else if (mod == 3){
           
##             outputMod3 <- list(infoPara = list(mu_beta=c(beta0,beta1),
##                                  Sigma_beta=diag(varInfoP),max_aG=30),
##                                beta0=beta0,beta1=beta1,scale=scale,shape=shape,mu=mu,sd=sd,pi=pi,K=2
##                                )
##           }
##       }
      
##     ## generate samples for dist = "b","b3", "g", and "l" (except "b2")
##     Y <- NULL
##     timecov <- c(0,0,(1:(ni-2)))
##     for ( ipat in 1 : Npat)
##       {
##         ## the number of repeated measures are the same
##         ## we assume that the time effects occurs once after the treatments are in effect
##         got <- rnbinom(ni,size = exp(beta0+beta1*timecov), prob = gsBASE[ipat])         
##         Y <- c(Y,got)
##       }
##     ##cat("beta0",beta0,"\n")
##     d <- data.frame(Y=Y,
##                     X1=rep(1,Npat*ni),
##                     X2=rep(timecov,Npat),
##                     ID=rep(1:Npat,each=ni),
##                     gsBASE = rep(gsBASE,each=ni),
##                     scan = rep(1:ni,Npat),
##                     ## day contains the day when the scan was taken
##                     ## 10 patients enter a trial every month
##                     days = days,
##                     hs = rep(hs,each=ni)
##                     )

##     ## dSMB visit is assumed to be every 4 months
##     if (full) return (d) 
##     d <- subset(d,subset= days <= DSMBVisitInterval*rev)
##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ DSMBVisitInterval*(rev-1) < d$days ] <- 1
##     ## The first two scans (screening and base-line scans are treated as pre-scans)
##     d$labelnp[ d$X2==0 ] <- 0

##     if (dist=="b2")
##       {
##         temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
##                              betas=c(beta0,beta1),aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))
##         d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp) ))
##       }else if (dist=="b")
##         {
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
##                                betas=c(beta0,beta1),aGs=c(aG1),rGs=c(rG1),pis=c(1))
##           d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
##         }else if (dist=="g"){

##           temp <- index.gammas(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
##                               betas=c(beta0,beta1),
##                               shapes=c(sc1,sc2), ## shape
##                               scales=c(sp1,sp2), ## scale
##                               pis=c(pi,1-pi))
##           d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
##         }else if (dist == "YZ"){
          
##           temp <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=cbind(d$X1,d$X2),
##                            betas=c(beta0,beta1),
##                            shape=shape, ## shape
##                            scale=scale, ## scale
##                            mu=mu,
##                            sd=sd,
##                            pi=pi)
##           d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp)) )
            
##           }

##     if (mod==0) return(d)
    
##     Npat <- length(unique(d$ID))
##     maxDs <- seq(0.01,5,0.01); Ktilde <- 2
##     ekn <- E.KN(maxDs=maxDs,N=Npat)
##     ib_D <- maxDs[min(which(ekn >= Ktilde))]
##     a_D <- 0.0001
##     if (mod==3){

##       outputMod3$infoPara$ib_D <- ib_D
##       outputMod3$infoPara$a_D <- a_D
##       ##outputMod3$infoPara$maxD <- 5
##       return(outputMod3)

##     }else if (mod == 4){
      
##       mu_beta <- rep(0,irev)
##       Sigma_beta <- diag(rep(5,irev))
      
##       outputMod4 <- list(
##                          mu_beta=mu_beta,Sigma_beta=Sigma_beta,
##                          ##a_qG=a_qG,r_qG=r_qG,
##                          a_D = a_D, ib_D=ib_D,
##                          ##maxD=5,
##                          max_aG = 30
##                          )
##       return(list(UinfoPara=outputMod4))
##     }
##   }


## Meta <- function(full.mean, full.Sigma, freq=ceiling(B/5),B=10000,tooLarge=1.5,
##                  Jacknife=TRUE){
##   ## function to compute the element-wise confidence intervals of correlation matrix for (beta, logD)
##   ## via bootstrap.

##   ## If the element-wise 95%CI of some elements of the correlation matrix are too large (i.e., larger than tooLarge),
##   ## such estimates are set to 0.
  
##   ## full.mean: a N by p matrix
##   ## full.Sigma: List with N element, each element is p by p matrix
##   ## freq: The frequency of printing
##   Nstudy <- nrow(full.mean)
##   p <- ncol(full.mean)

##   ## Obtain the point estimates of mean vec. sd vec. and covariances using all the Nstudy

##   if ( Jacknife ){
##     multfull <- mvmeta(full.mean,S=full.Sigma)
##     est.mean <- multfull$coefficient
##     est.sd <- sqrt(diag(multfull$Psi))
##     est.cov <- multfull$Psi
##     est.cor <- cov2cor(est.cov)
##   }else{
##     multfull <- multMeta(hat.betas=t(full.mean),hat.sigmas2=full.Sigma,want.mu=TRUE)
##     est.mean <- multfull$info.mean
##     est.sd <- sqrt(diag(multfull$info.var))
##     est.cov <- multfull$info.var
##     est.cor <- cov2cor(est.cov)
##   }

##   ## Bootstrap starts here!
##   if (Jacknife){
##     ## Storages to correct the bootstrap estimates of mean vec. sd vec. and covariance matrix
##     boot.cor <- matrix(NA,nrow=Nstudy,ncol=p^2)
##     boot.mean <- boot.sd <- matrix(NA,nrow=Nstudy,ncol=p)
##     oneToNstudy <- 1 : Nstudy
##     for (i.remove in 1 : Nstudy)
##       {
##         if (i.remove %% freq == 0) cat("\n sim ",ib)
##         boot.Sigmas <- list()
##         ib.temp <- 1
##         for (ib in oneToNstudy[-i.remove])
##           {
##             ## A list with Nstudy elements containing bootstrap p by p sigmas
##             boot.Sigmas[[ib.temp]] <- full.Sigma[[ib]]
##             ib.temp <- ib.temp + 1
##           }
##         mult <- mvmeta( full.mean[-i.remove,], S=boot.Sigmas)
##         boot.cor[i.remove,] <- (  c(cov2cor(mult$Psi)) - est.cor )^2
##         boot.sd[i.remove,] <- ( sqrt(diag(mult$Psi)) - est.sd )^2
##         boot.mean[i.remove,] <- ( mult$coefficient - est.mean )^2
##       }
##     jSD.boot.cor <- matrix( sqrt((Nstudy-1)*colMeans(boot.cor)),
##                            nrow=p,ncol=p)
##     jSD.boot.sd <- sqrt( (Nstudy-1)*colMeans(boot.sd) )
##     jSD.boot.mean <- sqrt( (Nstudy-1)*colMeans(boot.mean) )

##     untrust <- jSD.boot.cor > 0.25
##   }else{
##     ## Storages to correct the bootstrap estimates of mean vec. sd vec. and covariance matrix
##     boot.cor <- matrix(NA,nrow=B,ncol=p^2)
##     boot.mean <- boot.sd <- matrix(NA,nrow=B,ncol=p)
##     for (ib in 1 : B)
##       {
##         if (ib %% freq == 0) cat("\n sim ",ib)
##         samp.ind <- sample(1:Nstudy,Nstudy,replace=TRUE)
##         boot.Sigmas <- list()
##         for (i in 1 : Nstudy)
##           {
##             pick <- samp.ind[i]
##             ## A list with Nstudy elements containing bootstrap p by p sigmas
##             boot.Sigmas[[i]] <- full.Sigma[[pick]]
##           }
##         ## A Nstudy by p matrix. Each row contains bootstrap mean coefficient vector of length p
##         boot.means <- full.mean[samp.ind,]

##         ## Evaluate the mean vec., sd vec. and correlation matrix based on the bootstrap samples
##         ## mult <- mvmeta(boot.means,S=boot.Sigmas)
##         ## boot.cor[ib,] <- c(cov2cor(mult$Psi))
##         ## boot.sd[ib,] <- sqrt(diag(mult$Psi))
##         ## boot.mean[ib,] <- mult$coefficient
##         mult <- multMeta(hat.betas=t(boot.means),hat.sigmas2=boot.Sigmas,want.mu=TRUE)
##         boot.cor[ib,] <- c(cov2cor(mult$info.var))
##         boot.sd[ib,] <- sqrt(diag(mult$info.var))
##         boot.mean[ib,] <- mult$info.mean
##       }
##     cat("\n")
##     probs <- c(0.025,0.975)
##     CI.mean <- apply(boot.mean,2,quantile,probs=probs)
    
##     CI.sd <- apply(boot.sd,2,quantile,probs=probs)
    
##     CI.cor <- apply(boot.cor,2,quantile,probs=probs)
##     CIl.cor <- matrix(c(CI.cor[1,]),  nrow=p, ncol=p)
##     CIu.cor <- matrix(c(CI.cor[2,]),  nrow=p, ncol=p)
##     colnames(CIl.cor) <- rownames(CIl.cor) <- colnames(CIu.cor) <- rownames(CIu.cor) <- colnames(full.mean)
##     untrust <- CIu.cor-CIl.cor > tooLarge
##   }

##   ## The elements of covariance matrix whose bootstrap CI of corresponding correlations are tooLarge
##   ## receive zero estimates


##   est.cov.mod0 <- est.cov
##   est.cov.mod0[untrust] <- 0
##   ## 
##   if (sum(eigen(est.cov.mod0)$values < 0) > 0){
##     note <- "The modified covariance matrix which replaces un-trustable estimates with zero is not positive definite"  
##   }else{
##     note <- "The modified covariance matrix which replaces un-trustable estimates with zero is positive definite"
##   }
##   ## If the modified covariance matrix could be non-positive definite..
##   est.cov.mod <- adjustPosDef(est.cov.mod0)$adjust
##   muM <- muMult(hat.betas=t(full.mean),hat.sigmas2=full.Sigma,info.var=est.cov.mod)
##   est.mean.mod <- muM$info.mean
##   est.sd.mod <- sqrt(diag(est.cov.mod))

##   if (Jacknife){
##         return(
##            list(mean =
##                 list(est=est.mean,jSD = jSD.boot.mean,
##                      est.mod=est.mean.mod),
##                 sd   =
##                 list(est=est.sd,jSD = jSD.boot.sd,
##                      est.mod=est.sd.mod),
##                 cor  =
##                 list(est=est.cor,jSD = jSD.boot.cor,
##                      ## The elements of covariance matrix which
##                      est.mod0 = cov2cor(est.cov.mod0),
##                      est.mod = cov2cor(est.cov.mod)),
##                 cov  =
##                 list(est=est.cov,
##                      est.mod0 = est.cov.mod0,
##                      est.mod = est.cov.mod),
##                 para = list(Nstudy=Nstudy,B=B,tooLarge=tooLarge,note=note))
##            )
##   }else{
##     return(
##            list(mean =
##                 list(est=est.mean,lower=CI.mean[1,],upper=CI.mean[2,],
##                      est.mod=est.mean.mod),
##                 sd   =
##                 list(est=est.sd,lower=CI.sd[1,],upper=CI.sd[2,],
##                      est.mod=est.sd.mod),
##                 cor  =
##                 list(est=est.cor,lower=CIl.cor,upper=CIu.cor,
##                      ## The elements of covariance matrix which
##                      est.mod0 = cov2cor(est.cov.mod0),
##                      est.mod = cov2cor(est.cov.mod)),
##                 cov  =
##                 list(est=est.cov,
##                      est.mod0 = est.cov.mod0,
##                      est.mod = est.cov.mod),
##                 para = list(Nstudy=Nstudy,B=B,tooLarge=tooLarge,note=note))
##            )
##   }
## }

## multMeta <- function(hat.betas,hat.sigmas2,want.mu=TRUE)
##   {
##     ## hat.betas: Ndim by Ndat matrix, each column contains estimated betas from a study, missing values are accepted
##     ## hat.sigmas2: list of size Ndat, each element contains a estimated covariance matrix of estimator beta.
##     ##              The dimensions must agree with the maximum dimension of hat.betas among studies
##     ##              missing values are accepted

##     ## Compute correlation of hat.betas from their covariance matrix
##     corPhi <- lapply(hat.sigmas2,cov2corNA)

##     Ndim <- nrow(hat.betas)
##     Ndat <- ncol(hat.betas)
##     if (is.null(rownames(hat.betas)))  rownames(hat.betas) <- paste("beta",1:nrow(hat.betas),sep="")
    
##     ro.inv.phi2s <- inv.phi2s <- inv.phi4s <- Q <- bar.betas <- matrix(NA,nrow=Ndim,ncol=Ndim)
##     colnames(inv.phi4s) <- rownames(inv.phi4s) <- colnames(inv.phi2s) <- rownames(inv.phi2s) <-
##       colnames(Q) <- rownames(Q) <- colnames(bar.betas) <- rownames(bar.betas) <- rownames(hat.betas)
##     for (idim1 in 1 : Ndim)
##       {
##         betas <- hat.betas[idim1,] ## A vector of length N-study

##         ## inv.phi2 =|sum_s 1/(SE(h.bet1[s])*SE(h.bet1[s])) sum_s 1/(SE(h.bet1[s])*SE(h.bet2[s])) sum_s 1/(SE(h.bet1[s])*SE(h.bet3[s])) |
##         ##           |sum_s 1/(SE(h.bet2[s])*SE(h.bet1[s])) sum_s 1/(SE(h.bet2[s])*SE(h.bet2[s])) sum_s 1/(SE(h.bet2[s])*SE(h.bet3[s])) |
##         ##           |sum_s 1/(SE(h.bet3[s])*SE(h.bet1[s])) sum_s 1/(SE(h.bet3[s])*SE(h.bet2[s])) sum_s 1/(SE(h.bet3[s])*SE(h.bet3[s])) |
##         ## Extract the SD(hat.beta[idim1]) for each study then obtain the vector of 1/SE(hat.beta[idim1,istudy]) istudy=1,..,Ndat length Ndat
##         inv.phi1s.1 <- 1/sqrt(unlist(lapply(hat.sigmas2,function(x) x[idim1,idim1])))
##         for (idim2 in 1 : Ndim)
##           {
##             ## Extract the var(hat.beta[idim2]) similarly
##             inv.phi1s.2 <- 1/sqrt(unlist(lapply(hat.sigmas2,function(x) x[idim2,idim2])))
##             ## sum_{istudy=1}^{Ndat} 1/[ SE(hat.beta[idim1,istudy])*SE(hat.beta[idim2,istudy])]
##             inv.phi2s[idim1,idim2] <- sum(inv.phi1s.1*inv.phi1s.2,na.rm=TRUE)
##             ## inv.phi2 is Symmetric
##             ## will be used for the computation of var(beta) (out side of this loop)
##             ros <- unlist(lapply(corPhi,function(x) x[idim1,idim2]))
##             ro.inv.phi2s[idim1,idim2] <- sum(ros*(inv.phi1s.1*inv.phi1s.2),na.rm=TRUE)
##             inv.phi4s[idim1,idim2] <- sum((inv.phi1s.1*inv.phi1s.2)^2,na.rm=TRUE)
##             ## Contains bar.beta1 in its diagonal. The off-diagonal entries contain the estimates of bar.beta2 
##             bar.betas[idim1,idim2] <- sum(  betas*(inv.phi1s.1*inv.phi1s.2),na.rm=TRUE)/inv.phi2s[idim1,idim2]
##             ## bar.betas is NOT symmetric
##             ## bar.beta[i,j] = (sum_s h.beta_i[s]/(SE(h.bet_i[s])*SE(h.bet_j[s])) /[sum_s 1/{SE(h.bet_i[s])*SE(h.bet_j[s])}]
##           }
##       }
##     ## Compute the Q-statistics
##     ## (Computation of Q does not require the correlations among variables)
##     for (idim1 in 1 : Ndim)
##       {
##         betas1 <- hat.betas[idim1,]
##         phi1s.1 <- sqrt(unlist(lapply(hat.sigmas2,function(x) x[idim1,idim1])))
##         for (idim2 in 1 : Ndim)
##           {
##             betas2 <- hat.betas[idim2,]
##             phi1s.2 <- sqrt(unlist(lapply(hat.sigmas2,function(x) x[idim2,idim2])))
##             bar.beta1 <- bar.betas[idim1,idim2] 
##             bar.beta2 <- bar.betas[idim2,idim1]

##             Q[idim1,idim2] <- sum( (betas1-bar.beta1)*(betas2-bar.beta2)/(phi1s.1*phi1s.2),
##                                   na.rm=TRUE)
##           }
##       }

##     num <- Q - listSum(corPhi) + ro.inv.phi2s/inv.phi2s
##     den <- inv.phi2s - inv.phi4s/inv.phi2s
##     original.info.var <- num/den
##     ## original.info.var could be NOT positive definite
##     aPD<- adjustPosDef(original.info.var)
##     note <- aPD$note
##     info.var <- aPD$adjust.mat
##     colnames(info.var) <- rownames(info.var) <- rownames(hat.betas)
    
##     if (want.mu)
##       {
##         muM <- muMult(hat.betas=hat.betas,hat.sigmas2=hat.sigmas2,info.var=info.var)
##         info.mean <- muM$info.mean
##         SEmu <- muM$SEmu
##         ## varHatMu0 <- list()
##         ## nums <- matrix(NA,ncol=Ndat,nrow=Ndim)
##         ## dim.full <- nrow(hat.sigmas2[[1]])
##         ## for (istudy in 1 : Ndat)
##         ##   {
##         ##     sig <- hat.sigmas2[[istudy]]
##         ##     ## reduced matrix that only contains the non missing entries
##         ##     which.notNA <- !is.na(sig[,1])
##         ##     dim.red <- sum(which.notNA)
##         ##     sig.red <- matrix(sig[!is.na(sig)],nrow=dim.red)
##         ##     info.var.red <- matrix(info.var[!is.na(sig)],nrow=dim.red)
##         ##     inv.var <- solve(info.var.red + sig.red)
##         ##     varHatMu0[[istudy]] <- cbind(rbind(inv.var,matrix(NA,nrow=dim.full-dim.red,ncol=dim.red)),
##         ##                                  matrix(NA,nrow=dim.full,ncol=dim.full-dim.red))
##         ##     nums[which.notNA,istudy] <-inv.var%*%hat.betas[which.notNA,istudy]
##         ##   }
##         ## varHatMu <- solve(listSum(varHatMu0))
##         ## info.mean <- c( varHatMu%*%rowSums(nums,na.rm=TRUE))
        
##         ## SEmu <- sqrt(diag(varHatMu))
##         ## names(SEmu) <-  paste("SE:",rownames(hat.betas))
##         ## names(info.mean) <- rownames(hat.betas)
##       }else{
##         SEmu <- info.mean <- NULL
##       }
    
##     return(list(info.mean=info.mean,info.var=info.var,SEmu=SEmu,note=note,info.var0=original.info.var))
##   }

## getSNorm <- function(iseed = "random",
##                      rev = 4,
##                      dist = "b",
##                      mod = 0,
##                      probs = seq(0,0.99,0.01),
##                      ts = seq(0.001,0.99,0.001),
##                      trtAss = FALSE,
##                      trueCPI = FALSE,
##                      full=FALSE
##                      )
##   {
##     ## mod = 0: generate sample
##     ## mod = 1: quantiles of the true populations at given probs
##     ## mod = 2: densities of the true populations
##     ## mod = 3: parameters of the simulation model
##     ## dist = "b","b2","YZ" 
##     ## trtASS == FALSE then piTRT = 0 else piTRT = 0.674
##     if (trtAss) piTRT <- 0.6736842 else piTRT <- 0
##     ## if trueCPI == TRUE then returns the most precise conditional probability computed based on the
##     ## treatment assignment 
##     Npat <- 180; ## upto review 4, Npat=160 in total this number must be divisible by NpatEnterPerMonth

##     if (iseed=="random") set.seed(sample(1e+6,1)) else  set.seed(iseed)

##     ## The prior for regression coefficients beta

##     ## (Intercept) trt010:timeInt1 trt011:timeInt1 trt010:timeInt2 trt011:timeInt2
##     ## Estimated based on the 5 IFN-beta trials

##     ## mult.INF$info.mean
##     mu.beta <- round( c(1.38958538,  0.02353943,     -1.07433384,    -0.20072843,  -1.37750990,  -0.65131815),4 )

##     ## mult.INF$info.var
##     Sigma.beta <- matrix(round(c(0.08206687,      0.00000000,      0.00000000,      0.00000000,       0.0000000, -0.04902744,
##                                   0.00000000,      0.03136526,     -0.04367466,      0.00000000,       0.0000000,  0.00000000,
##                                   0.00000000,     -0.04367466,      0.06307239,      0.00000000,       0.0000000,  0.00000000,
##                                   0.00000000,      0.00000000,      0.00000000,      0.03811573,       0.0000000,  0.00000000,
##                                   0.00000000,      0.00000000,      0.00000000,      0.00000000,       0.1427152,  0.00000000,
##                                  -0.04902744,      0.00000000,      0.00000000,      0.00000000,       0.0000000,  0.03012234
##                                  ),4),nrow=length(mu.beta),byrow=TRUE)

##     if (mod == 3){
##       ## return parameter information of selected model
##       outputMod3 <- list(mu.beta=mu.beta, Sigma.beta = Sigma.beta)
##       mu.alpha <- mu.beta[1]
##       sigma.alpha <- Sigma.beta[1,1]
##     }
    
##     if (dist == "b"){
##       aG1 <- 3  
##       rG1 <- 0.8 ##0.8
##       if (mod %in% c(0,4:6)){
##         gs <- rbeta(Npat,aG1,rG1)
##         hs <- rep(1,Npat)
##       }else if (mod==1){
##         ## mod = 1: quantiles of the true populations at given probs
##         return(cbind(probs=probs,
##                      quantile=qbeta(probs,shape1=aG1,shape2=rG1)))
##       }else if (mod==2){
##         ## mod = 2: densities of the true populations
##         return (cbind(ts=ts,
##                       dens=dbeta(ts,shape1=aG1,shape2=rG1)) )
##       }else if (mod == 3){
##         ## mod = 3: parameters of the simulation model
##         c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- c(outputMod3,list(K=1,c1=c1,aGs=c(aG1=aG1), rGs=c(rG1=rG1)))
##       }

##     }else if (dist == "b2"){
##       ## mixture of two beta distributions
##       ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##       ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##       pi <- 0.3
      
##       aG1 <- 10
##       rG1 <- 10
##       aG2 <- 20
##       rG2 <- 1
##       ## generate the initial random effect values of everyone
##       if (mod %in% c(0,4:6)){
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gs <- rep(NA,Npat)
##         gs[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
##         gs[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
##       }else if (mod==1){
##         return (cbind(probs=probs,
##                       quantile=F_inv_beta2(ps=probs,aG1,rG1,aG2,rG2,pi=pi)))
##       }else if (mod==2){
##         return (cbind(ts=ts,
##                       dens=(pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2))))
##       }else if (mod == 3){

##         c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         c2 <- momBeta(aG2,rG2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- c(outputMod3,list(
##                                         K=2,c1=c1,c2=c2,
##                                         aGs=c(aG1=aG1,aG2=aG2),
##                                         rGs=c(rG1=rG1,rG2=rG2),
##                                         pi=pi))
##       }
##     }else if ( dist == "YZ"){
      
##       pi <- 0.85

##       alpha <- exp(-0.5)
##       ## a bimodal distribution with 85 % of Gi from a gamma distribution with mean 0.647 and variance 2.374
##       scale <- 2.374/0.647*alpha 
##       shape <- 0.647^2/2.374
##       mu <- 3*alpha
##       sd <- sqrt(0.25)*alpha
      
##       ## generate the initial random effect values of everyone
##       if (mod %in% c(0,4:6) ){
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gs <- rep(NA,Npat)
##         gs[hs==1] <- rgamma(Npat_dist1,shape=shape,scale=scale);
##         gs[hs==2] <- rnorm(Npat_dist2,mean=mu,sd=sd);
##         gs[gs < 0 ] <- 0
##         gs <- 1/(1+gs)
        
##       }else if (mod==1){
##         return (NULL)
##       }else if (mod == 2){
##         ts.trans <-1/ts-1
##         return (cbind(ts=ts,
##                       (pi*dgamma(ts.trans,shape=shape,scale=scale)+(1-pi)*dnorm(ts.trans,mean=mu,sd=sd) )*(1/ts^2)))
        
##       }else if (mod == 3){
##         c1 <- momGL(EG=shape*scale,VG=shape*(scale^2),mu.alpha=mu.alpha,sigma.alpha=sigma.alpha) ## expectation and variance of Y
##         c2 <- momGL(EG=mu,VG=sd^2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- outputMod3 <- c(outputMod3,list(c1=c1,c2=c2,scale=scale,shape=shape,mu=mu,sd=sd,pi=pi,K=2))
##       }
##     }
##     ## return parameter information of selected model
##     if (mod == 3) return (outputMod3)
    
##     ## Generate regression coefficients of this dataset
##     betasFull <- matrix(rmvnorm(n=1,mean=mu.beta[-length(mu.beta)],
##                                 sigma=Sigma.beta[-length(mu.beta),-length(mu.beta)]),ncol=1)
    
##     ## if (mod == 0) print(paste(c("Intercept",
##     ##       "timeInt1.plcb","timeInt1.trt",
##     ##       "timeInt2.plcb","timeInt2.trt"),
##     ##       sprintf("%1.4f",betasFull),collapse=":",sep=""))
##     ## Generate the treatment assignments 0 = placebo 1= trt
##     trtAss.pat <- sample(0:1,Npat,c(1-piTRT,piTRT),replace=TRUE)

##     ##print(trtAss.pat)
##     ##print(piTRT)
##     ## Sequential samples
##     ni <- 10
##     NpatEnterPerMonth <- 15
##     DSMBVisitInterval <- 4 ## months
    
##     ## === necessary for mod= 0, 3 and mod = 4
##     ## d contains a full dataset
##     days <- NULL
##     for (ipatGroup in 1 : (Npat/NpatEnterPerMonth))
##       {
##         ScandaysForSingleGroup <- ipatGroup:(ipatGroup+ni-1)
##         days <- c(days,rep(ScandaysForSingleGroup,NpatEnterPerMonth))
##       }
##     Y <- XforFit <- XforTCPI <- NULL
##     ## timeInt for all patients (used in model fit)
##     Xagg <- cbind(rep(1,ni),
##                   c(rep(0,2),rep(1,4),rep(0,ni-6)),
##                   c(rep(0,6),rep(1,ni-6)))
    
##     Xplcb <- cbind(Xagg[,1], ## ni = 9
##                    Xagg[,2],
##                    rep(0,ni), 
##                    Xagg[,3],
##                    rep(0,ni))
##     Xtrt <- cbind(Xagg[,1], ## ni = 9
##                   rep(0,ni), 
##                   Xagg[,2],
##                   rep(0,ni),
##                   Xagg[,3])
##     for ( ipat in 1 : Npat)
##       {
##         ## placebo
##         if (trtAss.pat[ipat]==0){
##           X <- Xplcb 
##         }else if (trtAss.pat[ipat]==1) X <- Xtrt ## trt
##         ## the number of repeated measures are the same
##         ## we assume that the time effects occurs once after the treatments are in effect
##         got <- rnbinom(ni,size = exp(X%*%betasFull), prob = gs[ipat])         
##         Y <- c(Y,got)
##         ## XforFit and XforTCPI must be the same if piTRT = 0
##         XforFit <- rbind(XforFit,Xagg)
##         if (trueCPI) XforTCPI <- rbind(XforTCPI,X)
##       }
##     colnames(XforFit) <- c("Intercept","timeInt1","timeInt2")

##     if (trtAss){
##       trtlabel <- factor(rep(trtAss.pat,each=ni),labels=c("plcb","trt"))
##     }else{
##       trtlabel <- factor(rep(trtAss.pat,each=ni),labels=c("plcb"))
##     }
##     d <- data.frame(Y=Y,
##                     XforFit,
##                     ID=rep(1:Npat,each=ni),
##                     gs = rep(gs,each=ni),
##                     scan = rep(-1:(ni-2),Npat),
##                     ## day contains the day when the scan was taken
##                     ## 10 patients enter a trial every month
##                     days = days,
##                     hs = rep(hs,each=ni),
##                     trtAss = trtlabel)

##     if (trueCPI){
##       d$XforTCPI <- XforTCPI
##     }
##     if (full) return (d) 
##     ## DSMB visit is assumed to be every 4 months
##     d <- subset(d,subset= days <= DSMBVisitInterval*rev)
##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ DSMBVisitInterval*(rev-1) < d$days ] <- 1
##     ## The first two scans (screening and base-line scans are treated as pre-scans)
##     d$labelnp[ d$scan <= 0 ] <- 0
##     betsAgg <- c(betasFull[1],
##                  aggBeta(bets=betasFull[2:3],piTRT=piTRT),
##                  aggBeta(bets=betasFull[4:5],piTRT=piTRT))
      
##     if (mod == 0){
##       XX <- cbind(d$Intercept,d$timeInt1, d$timeInt2)
##       if (dist=="b")
##         {
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                                betas=betsAgg,aGs=aG1,rGs=rG1,pis=1)
##           if (trueCPI){
##             tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                                  betas=betasFull,aGs=aG1,rGs=rG1,pis=1)
##           }
##         }else if (dist=="b2"){
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                                 betas=betsAgg,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))  
##           if (trueCPI){
##              tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                                   betas=betasFull,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))
##           }
##         }else if (dist == "YZ"){
          
##           temp <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                            betas=betsAgg,
##                            shape=shape, ## shape
##                            scale=scale, ## scale
##                            mu=mu,sd=sd,pi=pi)
##           if (trueCPI){
##             tCPI <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                              betas=betasFull,
##                              shape=shape, ## shape
##                              scale=scale, ## scale
##                              mu=mu,sd=sd,pi=pi)
##           }
##         }
      
##       d$betas <- c(betsAgg,rep(NA,nrow(d)-length(betsAgg)))
##       d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp) ))
##       if (trueCPI){
##         d$probIndexTRUE <- c(tCPI,rep(NA,nrow(d)-length(tCPI) ))
##       }
##       return(d)
##     }
##   }

## ## ========================================================

  
## Sigma2Cholv <- function(Sigma,loc.0=NULL,loc.SAME=NULL)
## {
##   if ( is.vector(loc.0)) loc.0 <- matrix(loc.0,nrow=1,ncol=4)
##   if ( is.vector(loc.SAME)) loc.SAME <- matrix(loc.SAME,nrow=1,ncol=4)
##   if ( !is.null(loc.SAME)){ if(ncol(loc.SAME)!= 4) stop()}

##   ## transform Sigma2 to vector of length sum(1:q) - Nzero containing the vectorized entries of
##   ## lower-triangle (nonzero) entries of the Cholesky decomposed Sigma
##   q <- nrow(Sigma)

##   Chol <- t(chol(Sigma))
##   if (!is.null(loc.0)){
##     for (irow in 1 : nrow(loc.0))
##       {
##         ir <- loc.0[irow,1]
##         ic <- loc.0[irow,2]
##         ## Sigma is symetric
##         Chol[ir,ic] <- Chol[ic,ir] <-  NA
##     }
##   }
##   if (!is.null(loc.SAME)){
##     for (irow in 1 : nrow(loc.SAME))
##       {
##         ## The Chol[loc.SAME[i,3], loc.SAME[i,4]] <- Chol[loc.SAME[i,1], loc.SAME[i,2] ]
##         ir <- loc.SAME[irow,3]
##         ic <- loc.SAME[irow,4]
##         ## Sigma is symetric
##         Chol[ir,ic] <- Chol[ic,ir] <- - Inf
##       }
##   }

##   ## diagonal entries of the Cholesky decomposed matrix must be positive
##   ## Cholv2Sigma takes exp() for all the diagonal entries
##   ## So Sigma2Cholv which is an inverse function of Cholv2Sigma,
##   ## must take log of all the diagonal entries of Cholvec.NA
##   diag(Chol) <- log( diag(Chol))
##   Cholvec.NA <- c(Chol[lower.tri(Chol, diag = TRUE)])
##   loc.vec.SAME <- Cholvec.NA == -Inf & !is.na(Cholvec.NA)
##   loc.vec.0 <- is.na(Cholvec.NA)
##   return(list(Cholvec=Cholvec.NA[ !loc.vec.SAME & !loc.vec.0],
##               q=q,
##               loc.vec.SAME = loc.vec.SAME,
##               loc.vec.0=loc.vec.0
##               )
##          )
## }


## Cholv2Sigma <- function(Cholvec,
##                         loc.0=NULL,loc.SAME=NULL,
##                         loc.vec.0,loc.vec.SAME,q)
## {
##   if (is.vector(loc.0)) loc.0 <- matrix(loc.0,nrow=1,ncol=2)
##   if (is.vector(loc.SAME))loc.SAME <- matrix(loc.SAME,nrow=1,ncol=4)
##   ## Cholvec: vector of sum(1:q)-Nzero, containing only the nonzero entries of vectored off-diagonal entries of
##   ##          Cholesky decomposed matrix
##   ## loc.0:  Nzero by 2 matrix, containing the location corresponding to zero in Sigma
##   Cholvec.NA <- rep(NA,sum(1:q))
##   Cholvec.NA[ !loc.vec.0 & !loc.vec.SAME ] <- Cholvec

##   Chol <- matrix(0,nrow=q,ncol=q)
##   Chol[lower.tri(Chol, diag = TRUE)] <- Cholvec.NA

##   ## Diagonal entries must be positive definite
##   diag(Chol) <- exp(diag(Chol))

##   ## Set 0 to the selected entries of covariance matrix noted in loc.0
##   if (! is.null(loc.0))
##     {
##       for (irow in 1 : nrow(loc.0))
##         {
##           ir <- loc.0[irow,1]
##           ic <- loc.0[irow,2]
##           if (ic > 1) ic.vec <- 1:(ic-1) else ic.vec <- NULL
##           Chol[ir,ic] <- - sum(Chol[ir,ic.vec]*Chol[ic,ic.vec])/Chol[ic,ic]
##         }
##     }
  
##   ## Keep the correlation of pairs of each row of Chol[irow,] being the same
##   if (! is.null(loc.SAME))
##     {
##       for (irow in 1 : nrow(loc.SAME))
##         {
##           ir1 <- loc.SAME[irow,1]
##           ic1 <- loc.SAME[irow,2]
##           ir2 <- loc.SAME[irow,3]
##           ic2 <- loc.SAME[irow,4]

##           ## The values located at (ir1,ic1) are copied into (ir2,ic2) 
##           ic1.vec <- 1:ic1
          
##           if (ic2 > 1) ic2.vec <- 1:(ic2-1) else ic2.vec <- NULL
##           ## Correlation term
##           if ( ir2 %in% c(ir1,ic1) )
##             {
##               if (ir2 == ir1){
##                 not.ir2 <- ic1
##               }else if(ir2 == ic1){
##                 not.ir2 <- ir1
##               } 
##               cor.scale <- sqrt(   sum(Chol[ic2,1:ic2]^2)/sum(Chol[not.ir2,1:not.ir2]^2)  )
                
##             }else{
##               ## sum(Chol[ir2,1:ir2]^2)*
##               ir2.vec <- (1:ir2)[! (1:ir2 %in% ic2)] 
##               cor.scale <- sqrt(sum(Chol[ic2,1:ic2]^2)/
##                                 (sum(Chol[ir1,1:ir1]^2)*sum(Chol[ic1,1:ic1]^2)))
##             }
##           C  <- cor.scale*(sum(Chol[ir1,ic1.vec]*Chol[ic1,ic1.vec])
##                            - sum(Chol[ir2,ic2.vec]*Chol[ic2,ic2.vec])
##                            )/Chol[ic2,ic2]
##           signC <-sign(C)
##           Chol[ir2,ic2] <- signC*C*sqrt(sum(Chol[ir2,ir2.vec]^2))/sqrt(1-C^2)
##         }
##     }
##     ## Set 0 to the selected entries of covariance matrix
##   while (sum(is.na(Chol)) > 0 )
##     {
##       if (! is.null(loc.0))
##         {
##           for (irow in 1 : nrow(loc.0) )
##             {
##               ir <- loc.0[irow,1]
##               ic <- loc.0[irow,2]
##               if (ic > 1) ic.vec <- 1:(ic-1) else ic.vec <- NULL
##               Chol[ir,ic] <- - sum(Chol[ir,ic.vec]*Chol[ic,ic.vec])/Chol[ic,ic]
##             }
##         }
##     }
##   ## Sigma
##   return(Chol%*%t(Chol))
## }


## Cholv2Sigma <- function(Cholvec,loc.0=NULL,loc.NAs,q)
## {
##   if (is.vector(loc.0)) loc.0 <- matrix(loc.0,nrow=1)
##   ## Cholvec: vector of sum(1:q)-Nzero, containing only the nonzero entries of vectored off-diagonal entries of
##   ##          Cholesky decomposed matrix
##   ## loc.0:  Nzero by 2 matrix, containing the location corresponding to zero in Sigma
##   ## loc.NAs: a vector of Nzero, containing the location corresponding to NA in Cholvec.NA
##   ##          where Cholvec.NA is an augmented vector of Cholvec with length sum(1:q)
##   Cholvec.NA <- rep(NA,sum(1:q))
##   Cholvec.NA[! loc.NAs] <- Cholvec

##   Chol <- matrix(0,nrow=q,ncol=q)
##   Chol[lower.tri(Chol, diag = TRUE)] <- Cholvec.NA
##   diag(Chol) <- exp(diag(Chol))
##   if (! is.null(loc.0))
##     {
##       for (irow in 1 : nrow(loc.0))
##         {
##           ir <- loc.0[irow,1]
##           ic <- loc.0[irow,2]
##           ic.vec <- (1:ic)[! (1:ic %in% ic)]
##           Chol[ir,ic] <- - sum(Chol[ir,ic.vec]*Chol[ic,ic.vec])/Chol[ic,ic]
##         }
##     }
##   ## Sigma
##   return(Chol%*%t(Chol))
## }



## listSum <- function(listobj,burnin=0)
##   {
##     objex <- listobj[[1]]
##     output <- rep(0,length(c(objex)))
##     for (i in (burnin+1) : length(listobj))
##       {
##         lsO <- c(listobj[[i]])
##         lsONA.TF <- !is.na(lsO)
##         output[lsONA.TF] <- output[lsONA.TF] + lsO[lsONA.TF]
##       }
##     return (matrix(output,
##                    nrow=nrow(objex),
##                    ncol=ncol(objex)))
##   }


## meta.rmlk <- function(ini.Psi=NULL,
##                       loc.0=NULL,loc.SAME=NULL,
##                       hat.sigmas2,hat.beta){

##   ## location
##   if (! is.null(loc.SAME)){
##     colnames(loc.SAME) <- c("irow_keep","icol_keep","irow_change","icol_change")
##     if (sum(loc.SAME[,1] <= loc.SAME[,2]) > 0 | sum(loc.SAME[,3] <= loc.SAME[,4]) >0 )
##       stop("The pairs of loc.SAME must be from off-diagonal lower-triangle matrix.")
##   }
##   if (is.null(ini.Psi))
##     { 
##       ini.Psi0 <- multMeta(hat.betas=hat.beta,hat.sigmas2=hat.sigmas2,want.mu=FALSE)$info.var
##       if (!is.null(loc.0)){
##         for (irow in 1 : nrow(loc.0))
##           {
##             ir <- loc.0[irow,1]
##             ic <- loc.0[irow,2]
##             ini.Psi0[ir,ic] <- ini.Psi0[ic,ir] <- 0
##           }
##         ini.Psi <- adjustPosDef(ini.Psi0)$adjust.mat
##       }
      
##       if (is.null(loc.0) & is.null(loc.SAME)){
##         ini.Psi <- ini.Psi0
##       }
##     }
  
##   temp <- Sigma2Cholv(Sigma=ini.Psi,loc.0=loc.0,loc.SAME=loc.SAME)
##   opt <- optim(par=temp$Cholvec,
##                fn=profile.nllk,
##                loc.0=loc.0,loc.SAME = loc.SAME,
##                loc.vec.0=temp$loc.vec.0,
##                loc.vec.SAME=temp$loc.vec.SAME,
##                q=temp$q,
##                hat.beta=hat.beta,hat.sigmas2=hat.sigmas2
##                )
##   Psi <- Cholv2Sigma(opt$par,
##                      loc.0=loc.0,
##                      loc.SAME=loc.SAME,
##                      q=temp$q,
##                      loc.vec.0=temp$loc.vec.0,
##                      loc.vec.SAME=temp$loc.vec.SAME)
##   if (!is.null(rownames(hat.beta)))rownames(Psi) <- colnames(Psi) <- rownames(hat.beta)
##   nllk <- opt$value
  
##   theta <- muMult(hat.betas=hat.beta,
##                   hat.sigmas2=hat.sigmas2,
##                   info.var=Psi,
##                   Ndat = ncol(hat.beta), Ndim=nrow(hat.beta))$info.mean
##   return(list(theta=theta,Psi=Psi,nllk=nllk,loc.SAME=loc.SAME))
## }



## profile.nllk <- function(Cholvec,
##                          hat.beta,
##                          loc.0,loc.SAME,
##                          loc.vec.0,loc.vec.SAME,
##                          hat.sigmas2,
##                          q=nrow(hat.beta),Nstudy =length(hat.sigmas2)
##                          ){
##   ## hat.beta is a q by Nstudy matrix

##   Psi <- Cholv2Sigma(Cholvec=Cholvec,
##                      loc.0=loc.0,loc.SAME=loc.SAME,
##                      loc.vec.0=loc.vec.0,loc.vec.SAME=loc.vec.SAME,
##                      q=q)
##   cat("\n\n Psi");print(Psi)
##   ## Matrix with large condition number is rejected because such matrix is computationally singular 
##   if (rcond(Psi) < 1e-15) return(Inf)
##   ## Update the overall mean with respect to the Psi
##   hat.theta <- muMult(hat.betas=hat.beta,
##                       hat.sigmas2=hat.sigmas2,
##                       info.var=Psi,
##                       Ndat = Nstudy, Ndim=q)$info.mean
  
##   Sigma.Psi <- lapply(hat.sigmas2,function(x) x + Psi)
##   ## A vector of length Nstudy
##   log.det.Sigma.Psi <- unlist(lapply(Sigma.Psi,function(x) log(det(x))))
##   ## A list of Nstudy elements containing q by q matrix
##   inv.Sigma.Psi <- lapply(Sigma.Psi,solve)
##   ## A vector of length Nstudy
##   det.inv.Sigma.Psi <- unlist(lapply(inv.Sigma.Psi, det))
##   mahara.sum <- 0
##   for (istudy in 1 : Nstudy)
##     {
##       disc <- matrix(hat.beta[,istudy] - hat.theta,ncol=1)
##       mahara.sum <- mahara.sum + t(disc)%*%inv.Sigma.Psi[[istudy]]%*%disc
##     }
##   tem1 <-  (Nstudy*q-q)*log(pi)
##   tem2 <-  sum(log.det.Sigma.Psi)
##   tem3 <-  log(sum(det.inv.Sigma.Psi))
##   nllk <- c( (tem1 + tem2 + tem3 + mahara.sum)/2 )
##   cat("\n nllk:",nllk)
##   return(nllk )
## }




## muMult <- function(hat.betas,hat.sigmas2,info.var,Ndat = length(hat.sigmas2),Ndim=nrow(hat.sigmas2[[1]]))
##   {
##     ## NA is not accepted
##     varHatMu0 <- list()
##     nums <- matrix(NA,ncol=Ndat,nrow=Ndim)
    
##      for (istudy in 1 : Ndat)
##        {
##          sig <- hat.sigmas2[[istudy]]
##          ## reduced matrix that only contains the non missing entries       
##          varHatMu0[[istudy]] <- inv.var <- solve(info.var + sig)
         
##          nums[,istudy] <-inv.var%*%hat.betas[,istudy]
##        }
##     varHatMu <- solve(listSum(varHatMu0))
##     info.mean <- c( varHatMu%*%rowSums(nums,na.rm=TRUE))
     
##     SEmu <- sqrt(diag(varHatMu))
##     names(SEmu) <-  paste("SE:",rownames(hat.betas))
##     names(info.mean) <- rownames(hat.betas)
##     return(list(SEmu=SEmu,info.mean=info.mean))
##   }
## ==============

## gsLabelnp <- function(d,olmeNBB,useSample){
##   ## d: simulated dataset from getSample
##   ## olmeNBB: must be the output of DPfit with model="nonpara"
  
##   ## This function returns the estimated random effects of patients at each time point.
##   ## the output is a vector of length sum ni
##   ## for example if a patient ipat has 5 scans and the first two are treated as pre-scans and the rest as new-scans,
##   ## then the d$labelnp[d$ID== opat] = 0,0,1,1,1
##   ## the output is output[d$ID==ipat] = gPre,gPre,gNew,gNew,gNew
##   ##            where gPre and gNew are the estimated random effects of patients 
##   ## nonparametric changing random effects
##   Npat <- length(unique(d$ID))
##   nolnp0TF <- tapply(d$labelnp==0,d$ID,sum) == 0 ## TRUE if patients do not have prescans d$labelnp==0
##   nolnp1TF <- tapply(d$labelnp==1,d$ID,sum) == 0 ## TRUE if patients do not have newscans d$labelnp==1
##   patwG2 <- which(! (nolnp0TF|nolnp1TF)) ## patients with both new scans and pre scans
  
##   gPres <- colMeans(olmeNBB$gPre[useSample,]); gPres_Index <- (1:Npat) - 0.5  ## everyone has g1
##   gNews <- colMeans(olmeNBB$gNew[useSample, olmeNBB$gNew[1,]!= -1000,drop=FALSE]); gNews_Index <- patwG2
##   ## combine gPres_gNews in the right order; in the order how they appear in the dataset d
##   gPres_gNews <- c(gPres,gNews)[order(c(gPres_Index,gNews_Index))]
##   g_long <- repeatAsID(values=gPres_gNews,ID=newCat(d$ID,d$labelnp))
##   return (g_long)
## }

## aggBeta <- function(bets,piTRT)
##   {
##     if (is.vector(bets)) bets <- matrix(bets,nrow=1)
##     log(exp(bets[,1])*(1-piTRT)+exp(bets[,2])*piTRT)
##   }

## rTildeBeta <- function(N.MC=1e+6,piTRT=0.6736842,
##               info.mean,info.var,logDMVN=FALSE){
##   samp <- rmvnorm(n=N.MC, mean = info.mean, sigma = info.var)
##   which.plcb <- seq(2,length(info.mean)-2,2)
##   which.trt <- which.plcb + 1
  
##   alpha <- samp[,1]
##   betas.plcb <- samp[,which.plcb]
##   betas.trt <- samp[,which.trt]

##   betatilde <- matrix(NA,nrow=N.MC,ncol=length(which.plcb))
##   for (itime1 in 1 : ncol(betatilde))
##     {
##       betatilde[,itime1]<-aggBeta(bets=cbind(betas.plcb[,itime1],
##                                     betas.trt[,itime1]),
##                                   piTRT=piTRT) 
##     }
##   re <- cbind(alpha,betatilde)

##   colnames(re) <- c("alpha",paste("tildeBeta",1:length(which.plcb),sep=""))

##   if (logDMVN){
##     logDs <- samp[,length(info.mean)]
##     re <- cbind(re,logD=logDs)
##   }
##   return(re)
## }


## aggBeta <- function(bets,piTRT)
##   {
##     if (is.vector(bets)) bets <- matrix(bets,nrow=1)
##     log(exp(bets[,1])*(1-piTRT)+exp(bets[,2])*piTRT)
##   }




## getS <- function(iseed = "random",
##                  rev = 4,
##                  dist = "b",
##                  mod = 0,
##                  probs = seq(0,0.99,0.01),
##                  ts = seq(0.001,0.99,0.001),
##                  full = FALSE,
##                  trtAss = FALSE,
##                  trueCPI = FALSE,
##                  IFN=FALSE)
##   {
##     ## mod = 0: generate sample
##     ## mod = 1: quantiles of the true populations at given probs
##     ## mod = 2: densities of the true populations
##     ## mod = 3: parameters of the simulation model
##     ## mod = 4: parameters for uninformative prior
##     ## mod = 5: parameters for informative prior under the placebo assumption
##     ## mod = 6: parameters for informative prior under the random sample assumption
##     ## dist = "b","b2","YZ" 
##     ## trtASS == FALSE then piTRT = 0 else piTRT = 0.674
##     if (trtAss) piTRT <- 0.6736842 else piTRT <- 0
##     ## if full = TRUE then full dataset is returned and rev is ignored
##     ## if trueCPI == TRUE then returns the most precise conditional probability computed based on the
##     ## treatment assignment 
##     Npat <- 180; ## upto review 4, Npat=160 in total this number must be divisible by NpatEnterPerMonth

##     if (iseed=="random") set.seed(sample(1e+6,1)) else  set.seed(iseed)

##     ## The prior for regression coefficients beta
##     if (IFN){
##       ## prior only based on IFN datasets
##       mu.beta <- c( 1.3863398,  0.0283592,  -1.0780041,  -0.1943629, -1.4005187)
##       Sigma.beta <- matrix(round(c( 0.07080114,  -0.03561139,   0.05512568,     -0.03559668,    0.08999295,
##                                    -0.03561139,   0.03365417,  -0.04243730,      0.02763631,   -0.05092112,
##                                     0.05512568,  -0.04243730,   0.06099776,     -0.02992216,    0.08058643,
##                                    -0.03559668,   0.02763631,  -0.02992216,      0.04543372,   -0.03128294,
##                                     0.08999295,  -0.05092112,   0.08058643,     -0.03128294,    0.13064430),4),
##                                  byrow=TRUE,5,5)
##     }else{
##       mu.beta <- c(1.29761146,  0.02762343,  -0.71042449,    -0.16054429,   -0.92415434)
##       Sigma.beta <- matrix(round(c( 0.06839318,   -0.0230865954,    -0.037123180,     -0.02048146,   -0.0697404833,
##                                    -0.02308660,    0.0173408469,    -0.000580806,      0.01097764,    0.0005858456,
##                                    -0.03712318,   -0.0005808060,     0.284828736,      0.02000585,    0.4364152673,
##                                    -0.02048146,    0.0109776373,     0.020005847,      0.02616263,    0.0312899440,
##                                    -0.06974048,    0.0005858456,     0.436415267,      0.03128994,    0.6989128826),4),
##                            byrow=TRUE,5,5)
##     }
##     if (mod == 3){
##       ## return parameter information of selected model
##       outputMod3 <- list(mu.beta=mu.beta, Sigma.beta = Sigma.beta)
##       mu.alpha <- mu.beta[1]
##       sigma.alpha <- Sigma.beta[1,1]
##     }
    
##     if (dist == "b"){
##       aG1 <- 3  
##       rG1 <- 0.8 ##0.8
##       if (mod %in% c(0,4:6)){
##         gs <- rbeta(Npat,aG1,rG1)
##         hs <- rep(1,Npat)
##       }else if (mod==1){
##         ## mod = 1: quantiles of the true populations at given probs
##         return(cbind(probs=probs,
##                      quantile=qbeta(probs,shape1=aG1,shape2=rG1)))
##       }else if (mod==2){
##         ## mod = 2: densities of the true populations
##         return (cbind(ts=ts,
##                       dens=dbeta(ts,shape1=aG1,shape2=rG1)) )
##       }else if (mod == 3){
##         ## mod = 3: parameters of the simulation model
##         c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- c(outputMod3,list(K=1,c1=c1,aGs=c(aG1=aG1), rGs=c(rG1=rG1)))
##       }

##     }else if (dist == "b2"){
##       ## mixture of two beta distributions
##       ## cluster 1: E(Y)=exp(beta0)*mu_G=exp(2)*((aG1+rG1-1)/(aG1-1)+1)=3.098636 ## at initial time
##       ## cluster 2: E(Y)=exp(beta0)*mu_G=exp(2)*((aG2+rG2-1)/(aG2-1)+1)=7.037196 ## at initial time
##       pi <- 0.3
      
##       aG1 <- 10
##       rG1 <- 10
##       aG2 <- 20
##       rG2 <- 1
##       ## generate the initial random effect values of everyone
##       if (mod %in% c(0,4:6)){
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gs <- rep(NA,Npat)
##         gs[hs==1] <- rbeta(Npat_dist1,aG1,rG1);
##         gs[hs==2] <- rbeta(Npat_dist2,aG2,rG2);
##       }else if (mod==1){
##         return (cbind(probs=probs,
##                       quantile=F_inv_beta2(ps=probs,aG1,rG1,aG2,rG2,pi=pi)))
##       }else if (mod==2){
##         return (cbind(ts=ts,
##                       dens=(pi*dbeta(ts,shape1=aG1,shape2=rG1)+(1-pi)*dbeta(ts,shape1=aG2,shape2=rG2))))
##       }else if (mod == 3){

##         c1 <- momBeta(aG1,rG1,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         c2 <- momBeta(aG2,rG2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- c(outputMod3,list(
##                                         K=2,c1=c1,c2=c2,
##                                         aGs=c(aG1=aG1,aG2=aG2),
##                                         rGs=c(rG1=rG1,rG2=rG2),
##                                         pi=pi))
##       }
##     }else if ( dist == "YZ"){
      
##       pi <- 0.85

##       alpha <- exp(-0.5)
##       ## a bimodal distribution with 85 % of Gi from a gamma distribution with mean 0.647 and variance 2.374
##       scale <- 2.374/0.647*alpha 
##       shape <- 0.647^2/2.374
##       mu <- 3*alpha
##       sd <- sqrt(0.25)*alpha
      
##       ## generate the initial random effect values of everyone
##       if (mod %in% c(0,4:6) ){
##         hs <- sample(1:2,Npat,c(pi,1-pi),replace=TRUE)
##         Npat_dist1 <- sum(hs==1)
##         Npat_dist2 <- sum(hs==2)
        
##         gs <- rep(NA,Npat)
##         gs[hs==1] <- rgamma(Npat_dist1,shape=shape,scale=scale);
##         gs[hs==2] <- rnorm(Npat_dist2,mean=mu,sd=sd);
##         gs[gs < 0 ] <- 0
##         gs <- 1/(1+gs)
        
##       }else if (mod==1){
##         return (NULL)
##       }else if (mod == 2){
##         ts.trans <-1/ts-1
##         return (cbind(ts=ts,
##                       (pi*dgamma(ts.trans,shape=shape,scale=scale)+(1-pi)*dnorm(ts.trans,mean=mu,sd=sd) )*(1/ts^2)))
        
##       }else if (mod == 3){
##         c1 <- momGL(EG=shape*scale,VG=shape*(scale^2),mu.alpha=mu.alpha,sigma.alpha=sigma.alpha) ## expectation and variance of Y
##         c2 <- momGL(EG=mu,VG=sd^2,mu.alpha=mu.alpha,sigma.alpha=sigma.alpha)
##         outputMod3 <- outputMod3 <- c(outputMod3,list(c1=c1,c2=c2,scale=scale,shape=shape,mu=mu,sd=sd,pi=pi,K=2))
##       }
##     }
##     ## return parameter information of selected model
##     if (mod == 3) return (outputMod3)
    
##     ## Generate regression coefficients of this dataset
##     betasFull <- matrix(rmvnorm(n=1,mean=mu.beta,sigma=Sigma.beta),ncol=1)
    
##     ## if (mod == 0) print(paste(c("Intercept",
##     ##       "timeInt1.plcb","timeInt1.trt",
##     ##       "timeInt2.plcb","timeInt2.trt"),
##     ##       sprintf("%1.4f",betasFull),collapse=":",sep=""))
##     ## Generate the treatment assignments 0 = placebo 1= trt
##     trtAss.pat <- sample(0:1,Npat,c(1-piTRT,piTRT),replace=TRUE)

##     ##print(trtAss.pat)
##     ##print(piTRT)
##     ## Sequential samples
##     ni <- 11
##     NpatEnterPerMonth <- 15
##     DSMBVisitInterval <- 4 ## months
    
##     ## === necessary for mod= 0, 3 and mod = 4
##     ## d contains a full dataset
##     days <- NULL
##     for (ipatGroup in 1 : (Npat/NpatEnterPerMonth))
##       {
##         ScandaysForSingleGroup <- ipatGroup:(ipatGroup+ni-1)
##         days <- c(days,rep(ScandaysForSingleGroup,NpatEnterPerMonth))
##       }
##     Y <- XforFit <- XforTCPI <- NULL
##     ## timeInt for all patients (used in model fit)
##     Xagg <- cbind(rep(1,ni),
##                   c(rep(0,2),rep(1,4),rep(0,ni-6)),
##                   c(rep(0,6),rep(1,ni-6)))
    
##     Xplcb <- cbind(Xagg[,1], ## ni = 9
##                    Xagg[,2],
##                    rep(0,ni), 
##                    Xagg[,3],
##                    rep(0,ni))
##     Xtrt <- cbind(Xagg[,1], ## ni = 9
##                   rep(0,ni), 
##                   Xagg[,2],
##                   rep(0,ni),
##                   Xagg[,3])
##     for ( ipat in 1 : Npat)
##       {
##         ## placebo
##         if (trtAss.pat[ipat]==0){
##           X <- Xplcb 
##         }else if (trtAss.pat[ipat]==1) X <- Xtrt ## trt
##         ## the number of repeated measures are the same
##         ## we assume that the time effects occurs once after the treatments are in effect
##         got <- rnbinom(ni,size = exp(X%*%betasFull), prob = gs[ipat])         
##         Y <- c(Y,got)
##         ## XforFit and XforTCPI must be the same if piTRT = 0
##         XforFit <- rbind(XforFit,Xagg)
##         if (trueCPI) XforTCPI <- rbind(XforTCPI,X)
##       }
##     colnames(XforFit) <- c("Intercept","timeInt1","timeInt2")

##     d <- data.frame(Y=Y,
##                     XforFit,
##                     ID=rep(1:Npat,each=ni),
##                     gs = rep(gs,each=ni),
##                     scan = rep(-1:(ni-2),Npat),
##                     ## day contains the day when the scan was taken
##                     ## 10 patients enter a trial every month
##                     days = days,
##                     hs = rep(hs,each=ni),
##                     trtAss =rep(trtAss.pat,each=ni))

##     if (trueCPI){
##       d$XforTCPI <- XforTCPI
##     }
##     ## DSMB visit is assumed to be every 4 months
##     if (full) return (d) 
##     d <- subset(d,subset= days <= DSMBVisitInterval*rev)
##     d$labelnp <- rep(0,nrow(d))
##     d$labelnp[ DSMBVisitInterval*(rev-1) < d$days ] <- 1
##     ## The first two scans (screening and base-line scans are treated as pre-scans)
##     d$labelnp[ d$scan <= 0 ] <- 0
##     betsAgg <- c(betasFull[1],
##                  aggBeta(bets=betasFull[2:3],piTRT=piTRT),
##                  aggBeta(bets=betasFull[4:5],piTRT=piTRT))
      
##     if (mod == 0){
##       XX <- cbind(d$Intercept,d$timeInt1, d$timeInt2)
##       if (dist=="b")
##         {
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                                betas=betsAgg,aGs=aG1,rGs=rG1,pis=1)
##           if (trueCPI){
##             tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                                  betas=betasFull,aGs=aG1,rGs=rG1,pis=1)
##           }
##         }else if (dist=="b2"){
##           temp <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                                 betas=betsAgg,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))  
##           if (trueCPI){
##              tCPI <- index.b.each(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                                   betas=betasFull,aGs=c(aG1,aG2),rGs=c(rG1,rG2),pis=c(pi,1-pi))
##           }
##         }else if (dist == "YZ"){
          
##           temp <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=XX,
##                            betas=betsAgg,
##                            shape=shape, ## shape
##                            scale=scale, ## scale
##                            mu=mu,sd=sd,pi=pi)
##           if (trueCPI){
##             tCPI <- index.YZ(Y=d$Y,ID=d$ID,labelnp=d$labelnp,X=d$XforTCPI,
##                              betas=betasFull,
##                              shape=shape, ## shape
##                              scale=scale, ## scale
##                              mu=mu,sd=sd,pi=pi)
##           }
##         }
      
##       d$betas <- c(betsAgg,rep(NA,nrow(d)-length(betsAgg)))
##       d$probIndex <- c(temp,rep(NA,nrow(d)-length(temp) ))
##       if (trueCPI){
##         d$probIndexTRUE <- c(tCPI,rep(NA,nrow(d)-length(tCPI) ))
##       }
##       return(d)
##     }
##     ## ======= mod == 4, 5 or 6 only get to here ===========
##     Npat <- length(unique(d$ID))
##     maxDs <- seq(0.01,5,0.01); Ktilde <- 2
##     ekn <- E.KN(maxDs=maxDs,N=Npat)
##     ib_D <- maxDs[min(which(ekn >= Ktilde))]
##     a_D <- 0.0001

##     if (mod == 4){
##       ## Uninformative priors 
##       outputMod4 <- list(a_D = a_D, ib_D=ib_D, max_aG = 30)
##       return(outputMod4)
##     }else if (mod== 5){
##       ## ======= mod == 5 or 6 only get to here ===========
##       ## Informative priors under the placebo assumption
##       if (IFN){
##         ## Data is developed based only on IFN
##         info.var <- matrix( round(c( 0.07080114,     -0.03561139,     -0.03559668,
##                                     -0.03561139,      0.03365417,      0.02763631,
##                                     -0.03559668,      0.02763631,      0.04543372
##                                     ),4),nrow=3,ncol=3)
##         info.mean <- c(1.3863398,       0.0283592,      -0.1943629)
##       }else{
##         ## Data is developed based on ALL datasets
##          info.var <- matrix( round(c( 0.06839318,     -0.02308660,     -0.02048146,
##                                      -0.02308660,      0.01734085,      0.01097764,
##                                      -0.02048146,      0.01097764,      0.02616263 
##                                      ),4),nrow=3,ncol=3)
##          info.mean <- c(1.29761146,    0.02762343,    -0.16054429)
##        }
##       modName <- "I.plcb"
##     }else if (mod==6){
##       ## Informative priors under the random sample assumption
##       if (IFN){
##         ## Data is developed based only on IFN
##         info.var <- matrix(round(c(0.070822670, 0.001545587, 0.013108930,
##                                    0.001545587, 0.002311349, 0.004566084,
##                                    0.013108930, 0.004566084, 0.023280450
##                                    ),4),nrow=3,ncol=3,byrow=TRUE)
##         info.mean <- c(1.3864922, -0.5497916, -0.8054787)
##       }else{
##         ## Data is developed based on ALL datasets
##         info.var <- matrix(round(c(0.06840140, -0.03005920, -0.04469143,
##                                   -0.03005920,  0.07694164,  0.11925741,
##                                   -0.04469143,  0.11925741,  0.20175429
##                                    ),4),nrow=3,ncol=3,byrow=TRUE)
##         info.mean <- c(1.2976255, -0.3690132, -0.5293237)
##       }
##       modName <- "I.RS"
##     }

##     if (sum(d$timeInt2) == 0)
##       {
##         pick.ind <- 1:2
##         info.var <- info.var[pick.ind,pick.ind]
##         info.mean <- info.mean[pick.ind]
##       }
##     infoPara <- list(max_aG=30, a_D=a_D, ib_D = ib_D,
##                      mu_beta=info.mean,
##                      Sigma_beta=info.var,
##                      mod=mod,modName=modName)
##     return(infoPara)
    
##   }



## lmeNBBayes <- function(formula,          ##   A vector of length sum ni, containing responses
##                        data,
##                        ##   A sum ni by p matrix, containing covariate values. The frist column must be 1 (Intercept)
##                        ID,         ##   A Vector of length sum ni, indicating patients
##                        B = 105000, ##     A scalar, the number of Gibbs iteration 
##                        burnin = 5000,  
##                        printFreq = B,
##                        M = NULL,
##                        probIndex = FALSE,
##                        thin =1, ## optional
##                        labelnp=NULL, ## necessary if probIndex ==1
##                        epsilonM = 1e-4,## nonpara
##                        para = list(mu_beta = NULL,Sigma_beta = NULL,max_aG=30),
##                        DP=TRUE,
##                        Reduce=1,
##                        thinned.sample=FALSE,
##                        Ddist=c("MVN","unif"),
##                        proposalSD = NULL
##                        ## Does not matter if DP=FALSE
##                        )
##   {
##     Xmodmat <- model.matrix(object=formula,data=data)
##     covariatesNames<- colnames(Xmodmat)
##     Y <- model.response(model.frame(formula=formula,data=data))
##     ## If ID is a character vector of length sum ni,
##     ## it is modified to an integer vector, indicating the first appearing patient
##     ## as 1, the second one as 2, and so on..

##     ## This code generate samples from NBRE model with constant random effect ~ DP mixture of Beta
##     if (is.vector(Xmodmat)) Xmodmat <- matrix(Xmodmat,ncol=1)
##     NtotAll <- length(Y)
##     if (nrow(Xmodmat)!= NtotAll) stop ("nrow(Xmodmat) != length(Y)")
##     if (length(ID)!= NtotAll)  stop ("length(ID)!= length(Y)")
##     if (!is.null(labelnp) & length(labelnp)!= NtotAll)  stop ("labelnp!= length(Y)")
##     if (thinned.sample & (!is.numeric(thin)) & (!is.numeric(burnin)))
##       stop("If you only want thinned samples, you must give the thinning and burnin parameters")
##     if (thinned.sample & (!Reduce))
##       stop("Non-reduced MCMC must return ALL samples")
##     if (length(Ddist) > 1) Ddist <- Ddist[1]
    
##     dims <- dim(Xmodmat)
##     Ntot <- dims[1]
##     pCov <- dims[2]
    
##     if (is.null(proposalSD)) proposalSD <- list() 
    
    
##     if (is.null(proposalSD$min$D))    proposalSD$min$D <- 0.02

##     if (is.null(proposalSD$max$D))    proposalSD$max$D <- 2
    
##     if (is.null(proposalSD$min$aG))   proposalSD$min$aG <- 0.05
##     if (is.null(proposalSD$min$rG))   proposalSD$min$rG <- 0.05
##     if (is.null(proposalSD$max$aG))   proposalSD$max$aG <- 5
##     if (is.null(proposalSD$max$rG))   proposalSD$max$rG <- 5
    
##     if (is.null(para$max_aG)) para$max_aG <- 30
##     if (DP){
      
##       if (Ddist == "unif")
##         {
##           if (is.null(para$mu_beta)) para$mu_beta <- rep(0,pCov)
##           if (is.null(para$Sigma_beta)){
##             para$Sigma_beta <- diag(10,pCov)
##           }
          
##           if (is.null(proposalSD$min$beta)) proposalSD$min$beta <- rep(diag(para$Sigma_beta)/1e+5,pCov)
##           if (is.null(proposalSD$max$beta)) proposalSD$max$beta <- rep(diag(para$Sigma_beta)/1e+3,pCov)
          
##           if (is.null(para$a_D ) ) para$a_D <- 1e-4
##           if (is.null(para$ib_D) ) para$ib_D <- 0.5
          
##           if (length(para$mu_beta)!=ncol(Xmodmat))
##             stop("The dimension of the fixed effect hyperparameter is wrong!")
          
##           if (is.null(M)) M  <- round(1 + log(epsilonM)/log(para$ib_D/(1+para$ib_D)))
          
##         }else if (Ddist == "MVN"){
          
##           EX <- 0.5
##           SDX <- 0.5
##           logDpara <- lnpara(EX=EX,SDX=SDX)
          
##           if (is.null(para$mu_beta)){ ## mu_beta is pCov + 1 length. The last entry corresponding to prior mean of log(D)
##             para$mu_beta <- rep(0,pCov)
##             para$mu_beta[pCov+1] <- logDpara$meanlog
##           }
##           if (is.null(para$Sigma_beta)){
##             para$Sigma_beta <- diag(10,pCov+1)
##             para$Sigma_beta[pCov+1,pCov+1] <- (logDpara$sdlog)^2
##           }
##           atlnD <- length(para$mu_beta)
##           if (is.null(proposalSD$min$beta)) proposalSD$min$beta <- rep(diag(para$Sigma_beta[-atlnD])/1e+5,pCov)
##           if (is.null(proposalSD$max$beta)) proposalSD$max$beta <- rep(diag(para$Sigma_beta[-atlnD])/1e+3,pCov)
##           if (is.null(proposalSD$min$D)) proposalSD$min$D <- rep(diag(para$Sigma_beta[atlnD])/1e+3,pCov)
##           if (is.null(proposalSD$max$D)) proposalSD$max$D <- rep(diag(para$Sigma_beta[atlnD])/1e+3,pCov)
          
          
##           if (length(para$mu_beta)!=(ncol(Xmodmat)+1))
##             stop("The dimension of the fixed effect and log(D) hyperparameter is wrong!")
          
##           if (is.null(M)){
##             for (M in 1 : 1000)
##               {
##                 EpiM <- piM(M=M,
##                             mean.norm=para$mu_beta[pCov+1],
##                             sd.norm=sqrt(para$Sigma_beta[pCov+1,pCov+1]))
                
##                 if (EpiM < epsilonM ) break
##               }
##           }
##         }
##     }else{
##       ## DP = FALSE
##       if (is.null(para$mu_beta)) para$mu_beta <- rep(0,pCov)
##       if (is.null(para$Sigma_beta))para$Sigma_beta <- diag(10,pCov)
##       if (length(para$mu_beta)!=ncol(Xmodmat))
##         stop("The dimension of the fixed effect hyperparameter is wrong!")
      
##     }
    
##     evalue_sigma_beta <- eigen(para$Sigma_beta, symmetric = TRUE, only.values = TRUE)$values
##     if (min(evalue_sigma_beta) <= 0) stop("Sigma_beta must be positive definite!")
##     Inv_sigma_beta <- c( solve(para$Sigma_beta) )

##     X <- c(Xmodmat) ## {xij} = { x_{1,1},x_{2,1},..,x_{Ntot,1},x_{1,2},....,x_{Ntot,p} }

##     ## change the index of ID to numeric from 1 to # patients
##     temID <- ID  
##     N <- length(unique(temID))
##     uniID <- unique(temID)
##     ID <- rep(NA,length(temID))
##     for (i in 1 : length(uniID))
##       {
##         ID[temID == uniID[i]] <- i
##       }
    
##     mID <- ID-1
##     ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
    
##     maxni <- max(tapply(rep(1,length(ID)),ID,sum))
##     Npat <- length(unique(ID))
    
##     if (probIndex)
##       {
##         ## the labelnp of patients only with 1 (new scans) labels are replaced by all 0 (old scans)
##         patwonew <- which(as.numeric(tapply((labelnp==0),ID,sum)==0)==1)
##         for (i in 1 : length(patwonew)) labelnp[ID == patwonew[i]] <- 0
        
##         patwoNorO <-  which(as.numeric(tapply((labelnp==1),ID,sum)==0)==1)
##         if (length(patwoNorO)==0) patwoNorO <- -1000;
##       }else{
        
##         labelnp <- rep(0,length(Y))
##       }
    
##     Btilde <- B
##     if (B %% thin != 0 ) stop("B %% thin !=0")
##     if (burnin %% thin !=0) stop("burnin %% thin !=0")
##     min_proposalSD <- c(aG=proposalSD$min$aG,rG=proposalSD$min$rG,beta=proposalSD$min$beta);
##     max_proposalSD <- c(aG=proposalSD$max$aG,rG=proposalSD$max$rG,beta=proposalSD$max$beta);
##     if (DP)
##       {
##         if (Reduce)
##           {
##             if (Ddist=="unif")
##               {

##                 re <- .Call("ReduceGibbs",
##                             as.numeric(Y),           ## REAL
##                             as.numeric(X),           ## REAL
##                             as.integer(mID),         ## INTEGER
##                             as.integer(B),           ## INTEGER
##                             as.integer(maxni),       ## INTEGER
##                             as.integer(Npat),        ## INTEGER
##                             as.numeric(labelnp),     ## REAL
##                             as.numeric(para$max_aG),
##                             as.numeric(para$mu_beta),     ## REAL
##                             as.numeric(evalue_sigma_beta),  ## REAL
##                             as.numeric(Inv_sigma_beta),  ## REAL
##                             as.numeric(para$a_D),
##                             as.numeric(para$ib_D),
##                             as.integer(M),
##                             as.integer(burnin),      ## INTEGER
##                             as.integer(printFreq),
##                             as.integer(probIndex),
##                             as.integer(thin),
##                             as.integer(thinned.sample),
##                             package = "lmeNBBayes"
##                             )
##               }else if (Ddist=="MVN")
##                 {
##                   ##cat("\n para: max_aG",para$max_aG)
##                   ##cat("\n theta:",para$mu_beta)
##                   cat("\n M:",M)
##                   min_proposalSD <- c( min_proposalSD,D=proposalSD$min$D)
##                   max_proposalSD <- c( max_proposalSD,D=proposalSD$max$D)
##                   re <- .Call("ReduceDmvn",
##                               as.numeric(Y),           ## REAL
##                               as.numeric(X),           ## REAL
##                               as.integer(mID),         ## INTEGER
##                               as.integer(B),           ## INTEGER
##                               as.integer(maxni),       ## INTEGER
##                               as.integer(Npat),        ## INTEGER
##                               as.numeric(labelnp),     ## REAL
##                               as.numeric(para$max_aG),
##                               as.numeric(para$mu_beta),     ## REAL
##                               as.numeric(evalue_sigma_beta),  ## REAL
##                               as.numeric(Inv_sigma_beta),  ## REAL
##                               as.integer(M),
##                               as.integer(burnin),      ## INTEGER
##                               as.integer(printFreq),
##                               as.integer(probIndex),
##                               as.integer(thin),
##                               as.integer(thinned.sample),
##                               as.double(min_proposalSD),
##                               as.double(max_proposalSD),
##                               package = "lmeNBBayes"
##                               )
##                 }
            
##             if (thinned.sample){
##               B <- (B - burnin)/thin
##             }
##             for ( i in 13 : 14 ) re[[i]] <- matrix(re[[i]],nrow=B,ncol=Npat,byrow=TRUE)
##             names(re) <- c("aGs","rGs","vs","weightH1",
##                            "condProb","h1s","g1s",
##                            "beta",
##                            "D","logL",
##                            "AR","prp","aGs_pat","rGs_pat")
##           }else{ ## Not Reduce
##             if (Ddist == "unif")
##               {
##                 re <- .Call("gibbs",
##                             as.numeric(Y),           ## REAL
##                             as.numeric(X),           ## REAL
##                             as.integer(mID),         ## INTEGER
##                             as.integer(B),           ## INTEGER
##                             as.integer(maxni),       ## INTEGER
##                             as.integer(Npat),        ## INTEGER
##                             as.numeric(labelnp),     ## REAL
##                             as.numeric(para$max_aG),
##                             as.numeric(para$mu_beta),     ## REAL
##                             as.numeric(evalue_sigma_beta),  ## REAL
##                             as.numeric(Inv_sigma_beta),  ## REAL
##                             as.numeric(para$a_D),
##                             as.numeric(para$ib_D),
##                             as.integer(M),
##                             as.integer(burnin),      ## INTEGER
##                             as.integer(printFreq),
##                             as.integer(probIndex),
##                             as.integer(thin),
##                             package = "lmeNBBayes"
##                             )
##                 names(re) <- c("aGs","rGs","vs","weightH1",
##                                "condProb","h1s","g1s",
##                                "beta",
##                                "D","logL",
##                                "AR","prp")
##                 re$para$a_D <- para$a_D
##                 re$para$ib_D <- para$ib_D
##               }
##           }
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         for ( i in 1 : 4 ) re[[i]] <- matrix(re[[i]],nrow=B,ncol=M,byrow=TRUE)
##         re[[5]]  <- matrix(re[[5]],nrow=(Btilde-burnin)/thin,ncol=Npat,byrow=TRUE)
##         for ( i in 6 : 7 ) re[[i]] <- matrix(re[[i]],nrow=B,ncol=Npat,byrow=TRUE)
##         re[[8]] <- matrix(re[[8]],nrow=B,byrow=TRUE)[,1:pCov] ## ncol= pCov or pCov + 1
##         if (probIndex)
##           {
##             ## patients with no new scans
##             if (sum(patwoNorO < 0) > 1) patwoNorO <- NULL
##             re$condProb[,patwoNorO] <- NA
##           }else{
##             re$condProb <- NULL
##           }
##         re$para$max_aG <- para$max_aG
##         re$para$M <- M
##         if (Ddist=="MVN")names(re$prp) <- names(re$AR) <-c("aG", "rG",paste("beta",1:pCov,sep=""),"D")

##       }else{ ## NOT DP

##         if (Reduce)
##           {

##             re <- .Call("Beta1reduce",
##                         as.numeric(Y),           ## REAL
##                         as.numeric(X),           ## REAL
##                         as.integer(mID),         ## INTEGER
##                         as.integer(B),           ## INTEGER
##                         as.integer(maxni),       ## INTEGER
##                         as.integer(Npat),        ## INTEGER
##                         as.numeric(labelnp),     ## REAL
##                         as.numeric(para$max_aG),
##                         as.numeric(para$mu_beta),     ## REAL
##                         as.numeric(evalue_sigma_beta),  ## REAL
##                         as.numeric(Inv_sigma_beta),  ## REAL
##                         as.integer(burnin),      ## INTEGER
##                         as.integer(printFreq),
##                         as.integer(probIndex),
##                         as.integer(thin),
##                         as.integer(thinned.sample),
##                         as.double(min_proposalSD),
##                         as.double(max_proposalSD),
##                         package = "lmeNBBayes"
##                         )
##             if (thinned.sample){
##               B <- (B - burnin)/thin
##             }
##           }else{
##             re <- .Call("Beta1",
##                         as.numeric(Y),           ## REAL
##                         as.numeric(X),           ## REAL
##                         as.integer(mID),         ## INTEGER
##                         as.integer(B),           ## INTEGER
##                         as.integer(maxni),       ## INTEGER
##                         as.integer(Npat),        ## INTEGER
##                         as.numeric(labelnp),     ## REAL
##                         as.numeric(para$max_aG),
##                         as.numeric(para$mu_beta),     ## REAL
##                         as.numeric(evalue_sigma_beta),  ## REAL
##                         as.numeric(Inv_sigma_beta),  ## REAL
##                         as.integer(burnin),      ## INTEGER
##                         as.integer(printFreq),
##                         as.integer(probIndex),
##                         as.integer(thin),
##                         package = "lmeNBBayes"
##                         )
##           }
##         ## http://stackoverflow.com/questions/8720550/how-to-return-array-of-structs-from-call-to-c-shared-library-in-r
##         ## for ( i in 1:4 ) re[[i]] <- matrix(re[[i]],B,M,byrow=TRUE)
##         re[[3]] <- matrix(re[[3]],nrow=B,ncol= Npat, byrow=TRUE) ## 
##         re[[4]] <- matrix(re[[4]],nrow=B,byrow=TRUE)[,1:pCov] ## ncol=pCov
##         re[[7]]  <- matrix(re[[7]],nrow=(Btilde-burnin)/thin,ncol=Npat,byrow=TRUE)
##         names(re) <- c("aG","rG",
##                        "g1s",
##                        "beta",
##                        "AR","prp",
##                        "condProb",
##                        "logL"
##                        )
##         if (probIndex)
##           {
##             ## patients with no new scans
##             if (sum(patwoNorO < 0) > 1) patwoNorO <- NULL
##             re$condProb[,patwoNorO] <- NA
##           }
##         re$para$max_aG <- para$max_aG
##         if (Reduce) names(re$prp) <- names(re$AR) <-c("aG", "rG",paste("beta",1:pCov,sep=""))
##         else{ names(re$prp) <- names(re$AR) <-c("aG", "rG","beta")}
##       }
##     if (thinned.sample){
##       thin <- 1
##       burnin <- 0
##     }

##     if (probIndex)
##       {
##         re$para$labelnp <- labelnp
##         re$condProbSummary <- condProbCI(ID,re$condProb)
##         rownames(re$condProbSummary) <- uniID
##       }
##     re$para$CEL <- Y
##     re$para$ID <- temID ## return original IDs
##     re$para$X <- Xmodmat
##     re$para$Sigma_beta <- para$Sigma_beta
##     re$para$mu_beta <- para$mu_beta
##     names(re$para$mu_beta[1:pCov]) <- rownames(re$para$Sigma_beta[1:pCov,])<-
##       colnames(re$para$Sigma_beta[,1:pCov] ) <- colnames(re$beta) <- covariatesNames
##     re$para$B <- B
##     re$para$burnin <- burnin
##     re$para$thin <- thin
##     re$para$probIndex <- probIndex
##     re$para$Reduce <- Reduce
##     re$para$burnin <- burnin
##     re$para$DP <- DP
##     re$para$thinned.sample <- thinned.sample
##     re$para$formula <- formula
##     class(re) <- "LinearMixedEffectNBBayes"
##     return (re)
##   }










## inpheatmap <- function(mat)
##   {
##     B = nrow(mat)
##     mat <- .Call("map_c",
##                  as.integer(c(mat)),
##                  as.integer(ncol(mat)),
##                  as.integer(B), package = "lmeNBBayes")/B
##     diag(mat) = 1;
##     return(mat);
##   }



## hm <- function(matBN,
##                aveCEL,
##                main="",
##                minbin=0.5)
##   {
##     ## This function is designed to summarize a dataset such that
##     ## the random effects of the first 100 pat are from dist1 and
##     ## the random effects of second 100 pat are from dist2.

##     ## Input is a B by N matrix, each row contains the cluster labels
##     ## Based on this input, first, compute the similarity matrix:
##     inphm <- inpheatmap(matBN)
##     m.hmps <- 1 - inphm
##     ord <- 1:ncol(matBN)

##     ## reorder the patients within dist1/dist2 based on the dissimilarity measure 
##     ## ord <- c(order.dendrogram(as.dendrogram(hclust(dist(m.hmps)))))

##     spacedID <- rep(NA,ncol(matBN))
##     for (i in 1 : ncol(matBN))
##       {
##         if (aveCEL[i]  < 0.5 )
##           spacedID[i] <- ""
##         if (aveCEL[i] >= 0.5 & aveCEL[i]  < 1 )
##           spacedID[i] <- "-"
##         else if (aveCEL[i] >= 1 & aveCEL[i]  < 2 )
##           spacedID[i] <- "--"
##         else if (aveCEL[i] >= 2 & aveCEL[i]  < 3 )
##           spacedID[i] <- "---"
##         else if (aveCEL[i] >= 3 & aveCEL[i]  < 4 )
##           spacedID[i] <- "----"
##         else if (aveCEL[i] >= 4 & aveCEL[i]  < 5 )
##           spacedID[i] <- "-----"
##         else if (aveCEL[i] >= 5 & aveCEL[i]  < 6 )
##           spacedID[i] <- "------"
##         else if (aveCEL[i] >= 6 & aveCEL[i]  < 7 )
##           spacedID[i] <- "-------"
##         else if (aveCEL[i] >= 7 & aveCEL[i]  < 8 )
##           spacedID[i] <- "--------"
##         else if (aveCEL[i] >= 8 & aveCEL[i]  < 9 )
##           spacedID[i] <- "---------"
##         else if (aveCEL[i] >= 9 & aveCEL[i]  < 10)
##           spacedID[i] <- "----------"
##         else if (aveCEL[i] >= 10)
##           spacedID[i] <- "-----------"
##       }
##     rownames(inphm) <- colnames(inphm) <- spacedID
    
##     ## Finally, plot the levelplot
##     library(lattice)
##     xmins <- c(0.01,0.5,0.01,0.5)
##     xmaxs <- c(0.49,1,0.49,1)
##     ymins <- c(0.5,0.5,0.01,0.01)
##     ymaxs <- c(1,1,0.51,0.51)
##     levelplot(inphm[ord,ord],
##               labels = list(cex=0.01),
##               at = do.breaks(c(minbin, 1.01), 20),
##               scales = list(x = list(rot = 90)),
##               aspect = "iso",
##               colorkey=list(text=3,labels=list(cex=2)),
##               xlab= "",ylab="",
##               region = TRUE, 
##               col.regions=gray.colors(100,start = 1, end = 0),
##               main=paste(main,sep="")
##               ##,par.settings=list(fontsize=list(text=15))
##               )

##     ## The posterior sample of label vector H1s contains the cluster labels.
##     ## The distinct advantage of our DP mixture model is its ability to classify the patients into 
##     ## the un-prespecified number of clusters.
##     ## 
##     ## Given B posterior samples of the label vectors of old scans H_1, for all combinations of the pairs of patients,
##     ## the average times that paired patients belong to the same cluster are recorded to create single Npat by Npat
##     ## similarity matrix. 
##     ## The level map is created with this similarity matrix where 
##     ## both row and column are reordered based on the result of hierarchical clustering on the similarity matrix (1-similarity matrix) 
##     ##
##   }

## proG1LeG2 <- function(gsBbyN,g2sBbyN,beta=TRUE)
##   {
##     if (is.vector(gsBbyN))
##       {
##         gsBbyN <- matrix(gsBbyN,ncol=1)
##         g2sBbyN <- matrix(g2sBbyN,ncol=1)
##       }
##     if(beta)
##       {
##         ## compare gPre >= gNew which is equivalent to 
##         ## (1-gPre)/gPre <= (1-gNew)/gNew
##         temp <- gsBbyN 
##         gsBbyN <-g2sBbyN
##         g2sBbyN <- temp
##       }
##     ##temp <- gsBbyN-g2sBbyN <= 0
##     ## pG1G2 <- colMeans(temp)
    
##     ##gsBbyN <- apply(gsBbyN,2,sort,decreasing=FALSE)
##     ##g2sBbyN <- apply(g2sBbyN,2,sort,decreasing=FALSE)
##     N <- ncol(gsBbyN)
##     pG1G2 <- colMeans(gsBbyN < g2sBbyN)
##     ## for (ipat in 1 : N)
##     ##   {
##     ##     g1 <- gsBbyN[,ipat]
##     ##     g2 <- g2sBbyN[,ipat]
##     ##     pG1G2[ipat] <-  .Call("pG1LeG2_c",
##     ##                           as.numeric(g1),
##     ##                           as.numeric(g2),
##     ##                           package = "lmeNBBayes"
##     ##                           )
##     ##   }
##     return(pG1G2)
##   }

## tempCmat <- function(A,B){
##   hi <- .Call("tempFun2",
##               as.numeric(c(A)), as.integer(nrow(A)),as.integer(ncol(A)),
##               as.numeric(c(B)), as.integer(ncol(B)),
##               package="lmeNBBayes")[[1]]
##   return(matrix(hi,nrow=Arow,ncol=Bcol))
## }
## tempC <- function(n,mu,Sigma){
##   temp <- eigen(Sigma)
##   evalue <- temp$value
##   evec <- c(temp$vector)
##   S <- matrix(NA,n,nrow(Sigma))
##   for (i in 1 : n)
##     S[i,] <- .Call("tempFun",as.numeric(mu),as.double(evalue),as.double(evec),package="lmeNBBayes")[[1]]
##   return (S)
## }
