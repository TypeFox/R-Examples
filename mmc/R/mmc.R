#' @export
#' @importFrom stats na.omit glm quantile formula sd
#' @importFrom survival coxph
#' @importFrom MASS mvrnorm
mmc <- function(model,type,mdata,rdata,rep,evar,rvar,bootstrap,boot=NULL){
  
  mdata <- na.omit(mdata)
  rdata <- na.omit(rdata)
  relib <- rdata[rvar] 
  mm <- formula(model)
  
  if(type=='cox' | type=='Cox' | type=='COX'){
    fitted <- coxph(formula=mm, mdata) 
  }
  
  if(type=='linear' | type=='Linear' | type=='LINEAR'){
    fitted <- glm(formula=mm, family='gaussian', mdata) 
  }
  
  if(type=='logistic' | type=='Logistic' | type=='LOGISTIC'){
    fitted <- glm(formula=mm, family='binomial', mdata) 
  }
  
  n <- dim(relib)[1]
  r <- rep
  k <- length(evar)
  
  rm <- matrix(rep(NA,n*k),ncol=k)
  colnum <- seq(from=1,to=r*k,by=r)
    
  for(i in 1:k){
    rm[,i] <-  rowMeans(relib[,colnum[i]:(colnum[i]+(r-1))])
  }
  
  rmr <- matrix(rep(NA,n*k*r),ncol=r*k)
  colnum2 <- rep(1:k, each=r)
  
  for(m in 1:(r*k)){
    rmr[,m] <- rm[,colnum2[m]]
  }
  
  difmat <- relib - rmr
  difmatsum <- t(difmat)%*%as.matrix(difmat)/(n*(r-1))
  colnum <- seq(from=1,to=r*k,by=r)
  within <- matrix(rep(NA,(k*k)),ncol=k)
  
  for(j in 1:k){
    for(i in 1:k){
      within[i,j] <- sum(diag(difmatsum[colnum[i]:(colnum[i]+(r-1)),colnum[j]:(colnum[j]+(r-1))]))
    }
  }
  
  mdataz <- data.frame(mdata[ ,which(colnames(mdata)%in%evar)]) # View(mdataz) # cov(mdataz)
  km <- ncol(mdataz)
  nm <- nrow(mdataz)
  difmatmdata <- matrix(rep(NA,nm*km),ncol=km)
  
  for(i in 1:km){
    difmatmdata[,i] <-  mdataz[,i] - mean(mdataz[,i])
  }
  
  total <- t(difmatmdata)%*%as.matrix(difmatmdata)/nm # View(total) # View(difmatmdata)
  
  between <- total - within
  
  coefficients <- summary(fitted)$coefficients
  if(k==1){
    est.coef <- coefficients[ which(rownames(coefficients)%in%evar), ][1]
  }
  
  if(k>1){
    est.coef <- coefficients[ which(rownames(coefficients)%in%evar), ][,1]
  }
  
  rinv <- diag(1,km) + within%*%solve(between)
  corrected <- t(as.matrix(est.coef))%*%rinv
  colnames(corrected) <- c(evar)
  rownames(corrected) <- c('Corrected coef')
  
  uncorrectedSummary <- summary(fitted)$coefficients
  coefsSummary <- c(t(coefficients[,1]),corrected)
  names(coefsSummary) <- c(rownames(coefficients),evar)
  
  output=list()
  
  output$uncorrected <- summary(fitted)$coefficients
  output$total <- total
  output$within <- within
  output$between <- between
  output$corrected <- t(corrected)
  
  if(bootstrap=='TRUE'){
    output$corrected <- NULL
    cat(noquote(paste("Now gerating ",boot," samples to obtain 95% confidence intervals for uncorrected and corrected estimates.",sep=""))) 
    cat("\n")
    
    ncolBoot <- dim(summary(fitted)$coefficients)[1] + k
    bootRep <- matrix(rep(NA,ncolBoot*boot),ncol=ncolBoot)
    
    for(tt in 1:boot){
      # tt <- 1
      x <- mdata[sample(1:nrow(mdata), nrow(mdata), replace=TRUE),] 
      
      if(type=='cox' | type=='Cox' | type=='COX'){
        fitted <- coxph(formula=mm, x) 
      }
      
      if(type=='linear' | type=='Linear' | type=='LINEAR'){
        fitted <- glm(formula=mm, family='gaussian', x) 
      }
      
      if(type=='logistic' | type=='Logistic' | type=='LOGISTIC'){
        fitted <- glm(formula=mm, family='binomial', x) 
      }
      
      mdataz <- data.frame(x[ ,which(colnames(x)%in%evar)])
      km <- ncol(mdataz)
      nm <- nrow(mdataz)
      difmatmdata <- matrix(rep(NA,nm*km),ncol=km)
      
      for(i in 1:km){
        difmatmdata[,i] <-  mdataz[,i] - mean(mdataz[,i])
      }
      
      
      total <- t(difmatmdata)%*%as.matrix(difmatmdata)/nm
      between <- total - within
      coefficients <- summary(fitted)$coefficients
      
      if(k==1){
        est.coef <- coefficients[ which(rownames(coefficients)%in%evar), ][1]
      }
      
      if(k>1){
        est.coef <- coefficients[ which(rownames(coefficients)%in%evar), ][,1]
      }
      
      rinv <- diag(1,km) + within%*%solve(between)
      corrected <- t(as.matrix(est.coef))%*%rinv
      
      bootRep[tt,] <-  c(coefficients[,1],corrected)
      
      if(tt/100==round(tt/100,0)){
        percentComplete <- (tt/boot)*100
        cat(noquote(paste("      Bootstrap sampling is ",percentComplete,"% complete.",sep=""))) 
        cat("\n")
      }
    }
    
    colnames(bootRep) <- c(rownames(coefficients),evar)
    bootCI <- matrix(rep(NA,ncolBoot*2),ncol=2)
    
    for(yy in 1:ncolBoot){
      bootCI[yy,] <- quantile(bootRep[,yy],c(.025,.975))
    }
    
    bootSE <- matrix(rep(NA,ncolBoot),ncol=1)
    for(ss in 1:ncolBoot){
      bootSE[ss,] <- sd(bootRep[,ss])
    }
    
    bootstrapped <- cbind(round(coefsSummary,6),round(bootSE,6),round(bootCI,6))
    
    vtype <- ifelse(rownames(summary(fitted)$coefficients)%in%evar, '(Uncorrected)','')
    vtype <- (c(vtype,rep('(Corrected)',length(evar))))
    rrnames <- paste(rownames(bootstrapped),vtype)
    rownames(bootstrapped) <- c(rrnames)
    colnames(bootstrapped) <- c('coef','SE','lower.limit','upper.limit')
    output$corrected <- bootstrapped
  }
  output
}