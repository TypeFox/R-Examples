##' @export
summary.bptwin <- function(object,level=0.05,transform=FALSE,...) {
  logit <- function(p) log(p/(1-p))
  tigol <- function(z) 1/(1+exp(-z))
  dlogit <- function(p) 1/(p*(1-p))
  dtigol <- function(z) tigol(z)^2*exp(-z)
  trnam <- " "
  vcoef1 <- paste("var(",c("A","C","D"),")",sep="")
  if (object$transform$invname!="") {
      vcoef1 <- paste(object$transform$invname,"(",vcoef1,")",sep="")
  }
  vcoef2 <- paste("atanh(",
                  c(paste("rho)","MZ",sep=trnam),
                    paste("rho)","DZ",sep=trnam)),sep="")
  idx1 <- na.omit(match(vcoef1,names(coef(object))))
  idx2 <- na.omit(match(vcoef2,names(coef(object))))
  CIs <- c()
  alpha <- level/2
  CIlab <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  V <- c()
  if (length(idx2)>0) {
    idx <- idx2
    V <- vcov(object)[idx,idx]
    arho <- coef(object)[idx2[1:2]]
    mz <- multinomlogit(coef(object)[idx2[1]]); names(mz) <- c("U","E")
    dz <- multinomlogit(coef(object)[idx2[2]]); names(dz) <- c("U","E")
    cc <- tanh(arho)
    names(cc) <- c("Tetrachoric correlation MZ","Tetrachoric correlation DZ")
    corMZ <- cc[1]; corDZ <- cc[2]
    D <- diag(object$tr$dtr(arho))
    h <- function(x) 2*(x[1]-x[2])
    dh <- function(x) c(2,-2)
    i1 <- 1:2
    corr <- NULL
  }
  if (length(idx1)>0) {
      idx <- idx1
      V <- vcov(object)[idx,idx]
      ACD <- match(names(coef(object))[idx1],vcoef1)
      nn <- c(c("A","C","D")[ACD],"E")
      dzsc <- c(1/2,1,1/4)[ACD]
      pp <- coef(object)[idx1]
      cc <- multinomlogit(pp,object$transform$tr,object$transform$dtr); names(cc) <- nn
      D <- attributes(cc)$gradient
      vcovACDE <- (D%*%V%*%t(D))
      if (transform) {
          cc2 <- logit(cc)
          D2 <- diag(dlogit(cc))
          DD <- D2%*%D
          Vc2 <- DD%*%V%*%t(DD)
          CIs <- tigol(cc2%x%cbind(1,1)+diag(Vc2)^0.5%x%cbind(-1,1)*qnorm(1-alpha))
      } else {
          CIs <- cbind(cc-qnorm(1-alpha)*sqrt(diag(vcovACDE)),cc+qnorm(1-alpha)*sqrt(diag(vcovACDE)))
      }
      K <- length(ACD)
      Ki <- seq_len(K)
      corMZ <- sum(cc[Ki]); corDZ <- sum(cc[Ki]*dzsc)
      i1 <- na.omit(match(c("D","A"),nn))
      h <- function(x) sum(x)
      dh <- function(x) rep(1,length(i1))
      ##    dh <- function(x) { res <- rep(0,length(x)); res[i1] <- 1; res }
      ##    h <- function(x) 2*(sum(x[i1])-sum(x[i1]*dzsc))
      ##    dh <- function(x) 2*(1-dzsc)
      
  }
  Vc <- D%*%V%*%t(D)
  datanh <- function(r) 1/(1-r^2)
  if (length(idx1)>0) {
    pp <- coef(object)[idx]
    b <- cbind(rep(1,K))
    corMZ.sd <- (t(b)%*%Vc[Ki,Ki]%*%b)[1]^0.5
    corDZ.sd <- (t(dzsc)%*%Vc[Ki,Ki]%*%dzsc)[1]^0.5    
    corr <- rbind(c(corMZ,corMZ.sd),c(corDZ,corDZ.sd))
    zrho <- atanh(corr[,1])
    zrho.var <- datanh(corr[,1])^2*corr[,2]^2
    corr <- cbind(corr, tanh(zrho%x%cbind(1,1)+zrho.var^0.5%x%cbind(-1,1)*qnorm(1-alpha)))
    rownames(corr) <- c("MZ Tetrachoric Cor","DZ Tetrachoric Cor")    
  } else {   
    zrho <- atanh(cc)
    zrho.var <- datanh(cc)^2*diag(Vc)
    CIs <- tanh(zrho%x%cbind(1,1)+zrho.var^0.5%x%cbind(-1,1)*qnorm(1-alpha))
  }
  newcoef <- rbind(cbind(cc,diag(Vc)^0.5,CIs),corr);
  ##  CIs <- rbind(CIs,c(NA,NA),c(NA,NA))
  ##  newcoef <- cbind(newcoef,CIs)
  colnames(newcoef) <- c("Estimate","Std.Err",CIlab)
  H <- h(cc[i1])
  hstd <- (t(dh(cc[i1]))%*%Vc[i1,i1]%*%dh(cc[i1]))^.5
  if (!transform) {
      ci <- c(H-hstd*qnorm(1-alpha),H+hstd*qnorm(1-alpha))
  } else {
      logith <- function(x) logit(h(x))  
      dlogith <- function(x) dlogit(h(x))*dh(x)  
      Dlh <- dlogith(cc[i1])
      sdlh <- (t(Dlh)%*%Vc[i1,i1]%*%(Dlh))[1]^0.5      
      suppressWarnings(ci <- tigol(logith(cc[i1]) + qnorm(1-alpha)*c(-1,1)*sdlh))
  }

  rhoOS <- NULL
  if (object$OS) {
    rEst <- object$coef[nrow(object$coef),1]
    rSE <- object$coef[nrow(object$coef),2]
    rhoOS <- tanh(rbind(c(rEst,rEst+qnorm(1-alpha)*c(-1,1)*rSE)))
    colnames(rhoOS) <- c("Estimate",CIlab)
    if (length(idx1)==0) {
      rownames(rhoOS) <- "Tetrachoric correlation OS"
    } else {
      rownames(rhoOS) <- "Kinship OS"    
    }
  }
      
  concordance <-  conditional <- marg <- c()
  probs <- function(p,idx=1) {
    if (idx==0) {
      m <- p[1]
      ##else m <- p[length(object$midx0)+1]
      S <- object$SigmaFun(p)
      conc1 <- pmvn(upper=c(m,m),sigma=S[[1]])
      conc2 <- pmvn(upper=c(m,m),sigma=S[[2]])
      marg <- pnorm(m,sd=S[[1]][1,1]^0.5)      
      return(logit((conc1-conc2)/(marg*(1-marg))))
    }
    S <- (object$SigmaFun(p))[[idx]]
    m <- 0
    if((object$npar$intercept==1 & idx==1) | object$eqmean) m <- p[1]
    else m <- p[length(object$midx0)+1]
    mu.cond <- function(x) m+S[1,2]/S[2,2]*(x-m)
    var.cond <- S[1,1]-S[1,2]^2/S[2,2]    
    conc <- pmvn(upper=c(m,m),sigma=S)
    disconc <- pmvn(lower=c(-Inf,m),upper=c(m,Inf),sigma=S)
    marg <- pnorm(m,sd=S[1,1]^0.5)
    cond <- conc/marg
    discond <- disconc/(1-marg)
    logOR <- log(cond)-log(1-cond)-log(discond)+log(1-discond)
    lambdaR <- cond/marg
    c(logit(c(conc,cond,marg)),lambdaR,logOR)
  }

  mycoef <- coef(object)
  ## formals(probs) <- alist(p=,idx=0)
  ## hp <- probs(mycoef)
  ## Dhp <- numDeriv::grad(probs,mycoef)
  ## shp <- diag(t(Dhp)%*%vcov(object)%*%(Dhp))^0.5
  
  formals(probs) <- alist(p=,idx=1)
  probMZ <- probs(mycoef)
  
  Dp0 <- numDeriv::jacobian(probs,mycoef)
  formals(probs) <- alist(p=,idx=2)
  probDZ <- probs(mycoef)
  Dp1 <- numDeriv::jacobian(probs,mycoef)
  sprobMZ <- diag((Dp0)%*%vcov(object)%*%t(Dp0))^0.5
  sprobDZ <- diag((Dp1)%*%vcov(object)%*%t(Dp1))^0.5
  probMZ <- cbind(probMZ,probMZ-qnorm(1-alpha)*sprobMZ,probMZ+qnorm(1-alpha)*sprobMZ)
  probMZ[1:3,] <- tigol(probMZ[1:3,])
  probDZ <- cbind(probDZ,probDZ-qnorm(1-alpha)*sprobDZ,probDZ+qnorm(1-alpha)*sprobDZ)
  probDZ[1:3,] <- tigol(probDZ[1:3,])
  rownames(probMZ) <- rownames(probDZ) <- c("Concordance","Casewise Concordance","Marginal","Rel.Recur.Risk","log(OR)")
  colnames(probMZ) <- colnames(probDZ) <- c("Estimate",CIlab)
 
  ## mu <- coef(object)[c(object$bidx0[1],object$bidx1[1])]
  ## Sigma <- list(object$Sigma0,object$Sigma1)
  ## for (i in 1:2) {
  ##   conc <- function()
  ##   mu.cond <- function(x) mu+Sigma[[i]][1,2]/Sigma[[i]][2,2]*(x-mu[i])
  ##   var.cond <- Sigma[[i]][1,1]-Sigma[[i]][1,2]^2/Sigma[[i]][2,2]    
  ##   cc0 <- pmvn(upper=c(mu[i],mu[i]),sigma=Sigma[[i]])
  ##   px <- pnorm(mu[i],sd=Sigma[[i]][2,2]^0.5)
  ##   concordance <- c(concordance,cc0)
  ##   marg <- c(marg,px)
  ##   conditional <- c(conditional,cc0/px)
  ## }
  ## names(concordance) <- names(conditional) <- c("MZ","DZ")
  
  hval <- rbind(c(H,hstd,ci)); colnames(hval) <- c("Estimate","Std.Err",CIlab);
  if (hval[1]>1) hval[1,] <- c(1,NaN,NaN,NaN)

##  hval <- rbind(hval, tigol(c(hp,NA,hp-qnorm(1-alpha)*shp,hp+qnorm(1-alpha)*shp)))
  rownames(hval) <- c("Broad-sense heritability")##,"Risk-scale Heritability")

  Nstr <- object$N
  nN <- ncol(object$N)
  ngroups <- ifelse(object$OS,3,2)
  postn  <- "MZ/DZ"; if (object$OS) postn <- paste(postn,"OS",sep="/")
  npos <- seq(ngroups)
  Nstr <- cbind(paste(Nstr[npos],collapse="/"),
                paste(Nstr[npos+ngroups],collapse="/"),
                paste(Nstr[npos+2*ngroups],collapse="/"))
  rownames(Nstr) <- ""
  colnames(Nstr) <-
    unlist(lapply(strsplit(colnames(object$N)[(1:3)*ngroups],".",fixed=TRUE),
                                  function(x) paste(x[1], postn)))

  all  <-  rbind(hval[,c(1,3,4),drop=FALSE],newcoef[,c(1,3,4),drop=FALSE])
  allprob <- rbind(probMZ,probDZ); rownames(allprob) <-
      c(paste("MZ",rownames(probMZ)),paste("DZ",rownames(probDZ)))
  all <- rbind(all,allprob)  
  
  cc <- object$coef; cc[,2] <- diag(vcov(object))^.5;
  cc[,3] <- cc[,1]/cc[,2]; cc[,4] <- 2*(pnorm(abs(cc[,1]/cc[,2]),lower.tail=FALSE))
  res <- list(heritability=hval,
              par=cc,
              probMZ=probMZ, probDZ=probDZ, Nstr=Nstr,
              rhoOS=rhoOS,
              coef=newcoef, all=all,
              vcov=vcov(object),
              AIC=AIC(object),
              time=attributes(object)$time,
              logLik=logLik(object)) ##, concordance=concordance, conditional=conditional)

  class(res) <- "summary.bptwin"
  res
}

##' @export
print.summary.bptwin <- function(x,digits = max(3, getOption("digits") - 2),...) {
  cat("\n")
  printCoefmat(x$par,digits=digits,...)
  cat("\n")
  ##  x$Nstr <- x$Nstr[,which((colnames(x$Nstr)!="Complete MZ/DZ")),drop=FALSE]
  NN <- x$Nstr[,c(1,3),drop=FALSE]; colnames(NN)[1] <- gsub("Complete ","",colnames(NN)[1])
  print(NN,quote=FALSE)
  cat("\n")
  cc <- rbind(x$coef[,-2,drop=FALSE],x$rhoOS)
  print(RoundMat(cc,digits=digits),quote=FALSE)
  cat("\nMZ:\n");
  print(RoundMat(x$probMZ,digits=digits),quote=FALSE)
  cat("DZ:\n")
  print(RoundMat(x$probDZ,digits=digits),quote=FALSE)
##  cat("\nConcordance (MZ; DZ):\t\t", x$concordance,"\n")
##  cat("Case-wise concordance (MZ; DZ):\t", x$conditional,"\n\n")
  cat("\n")
  print(RoundMat(x$heritability[,-2,drop=FALSE],digits=digits),quote=FALSE)
  cat("\n")
  if (!is.null(x$time)) {
      cat("\n")
      cat("Event of interest before time ", x$time, "\n", sep="")
  }

}
