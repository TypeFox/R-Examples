####
## Collection of functions extracting information from an object of class PK
##
####

estimator <- function(obj, se=FALSE){

  if(is.null(attr(obj,"class")) | attr(obj,"class")!="PK"){
    stop('Object supplied is not from class PK')
  }
  
  if(se){
    return(data.frame(est=obj$est,se=obj$CIs[obj$CIs[,'method']==obj$CIs[1,'method'],2]))
  }
  return(obj$est)

}

ci <- function(obj, method=NULL){

  if(is.null(attr(obj,"class")) | attr(obj,"class")!="PK"){
    stop('Object supplied is not from class PK')
  }
  if(is.null(method)){
    return(obj$CIs[,c('lower','upper')])
  }else{
    if(any(obj$CIs[,'method']==method[1])){
      return(obj$CIs[obj$CIs[,'method']==method[1],c('lower','upper')])
    }else{
      stop('Method(s) requested not part of the original estimation.')
    }
  }
}

test <- function(obj, theta=0, method=c("t","fieller", "z","resample"), nsample=1000) {

  method <- match.arg(method)

  if(is.null(attr(obj,"class")) | attr(obj,"class")!="PK"){
    stop('Object supplied is not from class PK')
  }
  resample <-FALSE
  if(method=="resample") resample <- TRUE
  if(!any(as.character(obj$CIs[,'method'])=="t") & method=="t"){ #
    warning("Method 't' requested, but not part of original call. Method 'z' used instead.")
    method <- "z"
  }
  if(!any(as.character(obj$CIs[,'method'])=="fieller") & method=="fieller"){
    warning("Method 'fieller' requested, but not part of original call. Method 'z' used instead.")
    method <- 'fieller'
  }
  if(!any(as.character(obj$CIs[,'method'])=="z") & method=="z"){ #
    stop("Method 'z' not part of original call. ")
  }
  if(method=="resample") method <- obj$CIs[1,'method']

  res <- obj$CIs[which(obj$CIs[,'method']==method),]  

  if(length(theta)>1 & length(theta)!=nrow(res)){
    stop('dimension of theta does not match the number of parameters to be tested')
  }

  stat <- (res[,1]-theta)/res[,2]
  if(method=="t" & !resample){
    df <- res[,'df']
    if(is.infinite(res[1,'lower'])){ ## greater than hypothesis test
      pval <- pt(stat,df=df)
      alternative <- 'greater'
    }else{
      if(is.infinite(res[1,'upper'])){ ## less than hypothesis test
        pval <- 1-pt(stat,df=df)
        alternative <- 'less'
      }else{
        pval <- 2*(1-pt(abs(stat),df=df))
        alternative <- 'two-sided'
      }
    }
  }else{
    df <- rep(NA,nrow(res))
    if(!resample){
      if(any(is.infinite(res[,'lower']))){ ## greater than hypothesis test
        pval <- pnorm(stat)
        alternative <- 'greater'
      }else{
        if(any(is.infinite(res[,'upper']))){ ## less than hypothesis test
          pval <- 1-pnorm(stat)
          alternative <- 'less'
        }else{
          pval <- 2*(1-pnorm(abs(stat)))
          alternative <- 'two-sided'
        }
      }
    }else{  # resampling based method
      if(any(is.infinite(res[,'lower']))){ ## greater than hypothesis test
        pval <- .resampling.test(obj=obj, theta=theta, alternative="greater", nsample=nsample)
        alternative <- 'greater'
      }else{
        if(any(is.infinite(res[,'upper']))){ ## less than hypothesis test
          pval <- .resampling.test(obj=obj, theta=theta, alternative="less", nsample=nsample) 
          alternative <- 'less'
        }else{
          pval <- .resampling.test(obj=obj, theta=theta, alternative="two.sided", nsample=nsample)
          alternative <- 'two-sided'
        }
      }
      method <- "resample"
    }
  }
  out <- NULL

  out$stat <- matrix(stat,ncol=1)
  out$p.value <- matrix(pval,ncol=1)
  out$theta <- matrix(theta,ncol=1,nrow=nrow(out$stat))
  rownames(out$stat) <- rownames(obj$est)
  rownames(out$p.value) <- rownames(obj$est)
  rownames(out$theta) <- rownames(obj$est)
  out$alternative <- alternative
  out$df <- df
  out$conf.level <- obj$conf.level
  out$design <- obj$design
  out$method <- method
  class(out) <- "PKtest"
  return(out)
}
