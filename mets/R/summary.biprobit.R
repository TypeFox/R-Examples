##' @export
summary.biprobit <- function(object,level=0.05,transform,contrast,mean.contrast=NULL,mean.contrast2=NULL,cor.contrast=NULL,marg.idx=1,...) {
  alpha <- level/2
  varcomp <- object$coef[length(coef(object)),1:2]
  varcomp <- rbind(object$model$tr(c(varcomp[1],varcomp[1]%x%cbind(1,1) + qnorm(1-alpha)*varcomp[2]%x%cbind(-1,1))))
  colnames(varcomp)[2:3] <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  rownames(varcomp) <- ifelse(is.null(object$model$varcompname),"Variance component",object$model$varcompname)
  if (!missing(contrast)) {
      contrast <- rbind(contrast)
      mean.contrast <- contrast[,seq(object$model$blen),drop=FALSE]
      cor.contrast <- contrast[,seq(object$model$zlen)+object$model$blen,drop=FALSE]
      if (!object$model$eqmarg) {          
          mean.contrast2 <- mean.contrast[,seq(ncol(mean.contrast)/2)+ncol(mean.contrast)/2,drop=FALSE]
          mean.contrast <- mean.contrast[,seq(ncol(mean.contrast)/2),drop=FALSE]
      }
  }
  
  h <- function(p) log(p/(1-p)) ## logit
  ih <- function(z) 1/(1+exp(-z)) ## expit
  ##dlogit <- function(p) 1/(p*(1-p))
  if (!missing(transform)) {
      h <- asin; ih <- sin      
      if (is.null(transform)) {
          h <- ih <- identity
      }
      if (is.list(transform)) {
          h <- transform[[1]]; ih <- transform[[2]]
      }     
  }
  convval <- function(val) {
      i1 <- which(val==1)
      i2 <- which(val==-1)
      val <- as.character(val)
      val[seq_len(length(val)-1)+1] <-
          paste("+ ", val[seq_len(length(val)-1)+1],sep="")
      val[i2] <- "- ";
      val[i1] <- "+ "
      val[intersect(1,i1)] <- ""
      return(val)
  }
  parfun <- function(p,ref=FALSE,mean.contrast,mean.contrast2,cor.contrast) {
      nn <- paste("[",gsub("r:","",rownames(object$coef),fixed=TRUE),"]",sep="")
      m <- rep(p[1],2)
      r <- p[object$model$blen+1]
      corref <- mref1 <- mref2 <- NULL
      if (ref) {
          corref <- nn[object$model$blen+1]
          mref1 <- mref2 <- nn[1]
      }
      if (!is.null(mean.contrast)) {
          m[1] <- sum(p[seq_along(mean.contrast)]*mean.contrast)
          if (ref) {
              idx1 <- which(mean.contrast!=0)
              mref <- nn[idx1]              
              mref1 <- mref2 <- paste(convval(mean.contrast[idx1]),mref,sep="")
          }
          if (is.null(mean.contrast2) && !object$model$eqmarg) {
              idx1 <- which(mean.contrast!=0)
              idx2 <- idx1+1
              if (length(object$npar$pred)>0) idx2 <- idx2+object$npar$pred/2
              mean.contrast2 <- rep(0,object$model$blen)
              mean.contrast2[idx2] <- mean.contrast[idx1]              
          }
      }
      if (!is.null(mean.contrast2)) {
          m[2] <- sum(p[seq_len(object$model$blen)]*mean.contrast2[seq_len(object$model$blen)])
          if (ref) {
              idx1 <- which(mean.contrast2!=0)
              mref <- nn[idx1]
              mref2 <- paste(convval(mean.contrast2[idx1]),mref,sep="")
          }          
      } else {
          if (object$model$eqmarg) { m <- rep(m[1],2) }          
      }
      if (!object$model$eqmarg & is.null(mean.contrast) & is.null(mean.contrast2)) {
          idx <- 2
          if (length(object$npar$pred)>0) idx <- object$npar$pred/2+1
          m[2] <- p[idx]
          mref2 <- nn[idx]
      }
      if (!is.null(cor.contrast)) {
          p.idx <- seq_len(object$model$zlen)+object$model$blen
          if (length(cor.contrast)==length(p)) p.idx <- seq(length(p))          
          r <- sum(p[p.idx]*cor.contrast)
          if (ref) {
              idx1 <- which(cor.contrast!=0)
              corref <- nn[p.idx[idx1]]
              corref <- paste(convval(cor.contrast[idx1]),corref,sep="")
          }
      }
      return(list(m=m,r=r,mref1=mref1,mref2=mref2,corref=corref))
  }
  probs <- function(p,...) {
      pp <- parfun(p,...)
      m <- pp[[1]]
      r <- pp[[2]]
      S <- object$SigmaFun(r,cor=FALSE)
      ##mu.cond <- function(x) m[1]+S[1,2]/S[2,2]*(x-m[2])
      ##var.cond <- S[1,1]-S[1,2]^2/S[2,2]
      p11 <- pmvn(lower=c(0,0),mu=m,sigma=S) 
      p01 <- pmvn(lower=c(-Inf,0),upper=c(0,Inf),mu=m,sigma=S)
      p10 <- pmvn(lower=c(0,-Inf),upper=c(Inf,0),mu=m,sigma=S)
      p00 <- 1-p11-p10-p01
      marg1 <- p11+p10
      marg2 <- p11+p01
      cond1 <- p11/marg2
      lambda <- cond1/marg1
      discond1 <- p10/(1-marg2)
      logOR <- log(cond1)-log(1-cond1)-log(discond1)+log(1-discond1)
      ##rho <- S[1,2]/S[1,1]
      if (object$model$eqmarg) {
          return(c(h(c(p11,cond1,marg1)),lambda,logOR,r))
      }
      return(c(h(c(p11,p10,p01,p00,marg1,marg2)),logOR,r))
  }  
  alpha <- level/2
  CIlab <- paste(c(alpha*100,100*(1-alpha)),"%",sep="")
  mycoef <- coef(object)

  cor.contrast <- rbind(cor.contrast)
  mean.contrast <- rbind(mean.contrast)
  mean.contrast2 <- rbind(mean.contrast2)
  KK <- lapply(list(cor.contrast,mean.contrast,mean.contrast2),nrow)
  if (all(is.null(unlist(KK)))) K <- 1 else  K <- max(unlist(KK))
  res <- pa <- c()
  for (i in seq(K)) {
      prob <- probs(mycoef,cor.contrast=cor.contrast[i,],mean.contrast=mean.contrast[i,],mean.contrast2=mean.contrast2[i,])  
      Dprob <- numDeriv::jacobian(probs,mycoef,cor.contrast=cor.contrast[i,],mean.contrast=mean.contrast[i,],mean.contrast2=mean.contrast2[i,])
      sprob <- diag((Dprob)%*%vcov(object)%*%t(Dprob))^0.5
      pp <- cbind(prob,prob-qnorm(1-alpha)*sprob,prob+qnorm(1-alpha)*sprob)
      pp[nrow(pp),] <- object$model$tr(pp[nrow(pp),])
      pp[nrow(pp)-1,] <- exp(pp[nrow(pp)-1,])      
      if (!object$model$eqmarg) {
          pp[1:6,] <- ih(pp[1:6,])      
          nn <- c("P(Y1=1,Y2=1)","P(Y1=1,Y2=0)","P(Y1=0,Y2=1)","P(Y1=0,Y2=0)","P(Y1=1)","P(Y2=1)","OR","Tetrachoric correlation")
      } else {
          pp[1:3,] <- ih(pp[1:3,])      
          nn <- c("Concordance","Casewise Concordance","Marginal","Rel.Recur.Risk","OR","Tetrachoric correlation")
      }
      if (K>1) nn <- paste("c",i,":",nn,sep="")
      if (nrow(pp)-length(nn)>0) nn <- c(nn,rep("",nrow(pp)-length(nn)))
      rownames(pp) <- nn
      colnames(pp) <- c("Estimate",CIlab)
      
      P <- nrow(pp)
      pa <- c(pa, list(parfun(object$coef[,1],ref=TRUE,cor.contrast=cor.contrast[i,],mean.contrast[i,],mean.contrast2[i,])))
      res <- rbind(res,pp)
  }      
  
  contrast <- any(c(!is.null(cor.contrast),!is.null(mean.contrast),!is.null(mean.contrast2)))
  res <- list(all=res,varcomp=varcomp,prob=res,coef=object$coef,score=colSums(object$score),logLik=object$logLik,msg=object$msg,N=object$N,ncontrasts=K,nstat=P,
              par=pa,model=object$model,contrast=contrast, time=attributes(object)$time)
  class(res) <- "summary.biprobit"
  res
}
