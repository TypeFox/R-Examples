

######################################################################################
### Purpose: the wrapper of S3 method "summary" for object "singletable"
### Input:   S3 object "singletable" 
### Output:  return value of summary.single (refer to help file)
### Note:    the implement is function "summary_single_sampling"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
######################################################################################
summary.singletable <- function(object,...) {
  if (!inherits(object, "singletable"))
    stop("Use only with 'singletable' objects.\n")
	result <- summary_single_sampling(object)
  invisible(result)
}


######################################################################################
### Purpose: the implement of S3 method "summary" for object "multipletables"
### Input:   S3 object "multipletable" 
### Output:  return value of summary.multiple (refer to help file)
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
######################################################################################
summary.multipletables <- function(object,...) {
  if (!inherits(object, "multipletables"))
    stop("Use only with 'multipletables' objects.\n")

  model <- object$model
  measure <- object$measure
  measurename <- object$measurename
  cov.matrix <- object$cov.matrix
  overall <- object$overall
  alpha <- object$alpha
  specific <- study_specifc(object)
  studynames <- object$studynames

  siglevel <- paste(as.character((1-alpha)*100),"%",sep="")

  if (model=="Sarmanov") cat("Model: Sarmanov Beta-Binomial Model",fill=TRUE)
  if (model=="Independent") cat("Model: Independent Beta-Binomial Model",fill=TRUE)
   
  cat("Overall",measurename,fill=TRUE)
  cat("Estimate:", round(overall$overall,3),fill=TRUE)
  cat(siglevel," CI:[",round(overall$CI[1],3),",", round(overall$CI[2],3),"]", sep="",fill=TRUE)
  cat("",fill=TRUE)
  cat("Maximum likelihood estimates of hyperparameters:",fill=TRUE)

  if (model=="Sarmanov"|model=="Independent") {
    chi2 <- object$chi2
    p.value <- object$pvalue
    a1 <- object$MLE[1]; b1 <- object$MLE[2]
    a2 <- object$MLE[3]; b2 <- object$MLE[4]
    rho <- object$MLE[5]

    cat("a1 =",round(a1,3),", ","b1 =",round(b1,3),", ","a2 =",round(a2,3),
        ", ","b2 =",round(b2,3),", ","rho =",round(rho,3),sep="",fill=TRUE)
    cat("Likelihood ratio test for within-group correlation (H0: rho=0):",sep="",fill=TRUE)
    cat("chi2: ",round(chi2,3),"; ","p-value: ",round(p.value,3),sep="",fill=TRUE)

    out <- list(model=model,measure=measure,cov.matrix=cov.matrix,hessian=object$hessian,
                overall=overall,studynames=studynames,chi2=chi2,pvalue=p.value,
                alpha=alpha,MLE=object$MLE,studyspecific=specific)
  }
  cat("Study-Specifc",measurename,":",sep="",fill=TRUE)
  print(round(specific[[1]],3))
  cat("\n")  
  invisible(out)
}

######################################################################################
### Purpose: the implement of S3 method "summary" for object "singletable"
### Input:   S3 object "singletable" 
### Output:  return value of summary.single (refer to help file)
### Note:    the wrapper is function "summary.single"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
######################################################################################
summary_single_sampling <- function(object) {
  alpha <- object$alpha
  measurename <- object$measurename
  model <- object$model
  measure <- object$measure
  siglevel <- paste(as.character((1-alpha)*100),"%",sep="")
  measure<-object$measure
  studynames<-object$studynames
  
  cat("Measure:",measurename,fill=TRUE)
  if (model=="Sarmanov") {
    cat("Model: Sarmanov Beta-Binomial Model",fill=TRUE)
    cat("Prior:",model,"\n")
  }
  if (model=="Independent")  cat("Model: Independent Beta-Binomial Model",fill=TRUE)
  cat("\n",fill=TRUE)
     
  result <- list()

  result[[1]] <- c(mean(object$sample[[1]]),median(object$sample[[1]]),
                   quantile(object$sample[[1]],probs = c(alpha/2,1-alpha/2),na.rm=TRUE),
                   minlength.CI(object$sample[[1]],alpha))

  names(result[[1]]) <- c("Posterior Mean","Posterior Median",paste(siglevel,"ET CI left"),
                          paste(siglevel,"ET CI right"),paste(siglevel,"HDR CI left"),
                          paste(siglevel,"HDR CI right"))

  ## print the object
  cat("Mean:",round(result[[1]][1],3),fill=TRUE)
  cat("Median:",round(result[[1]][2],3),fill=TRUE)
  cat(siglevel," ET CI: [",round(result[[1]][3],3),"," ,round(result[[1]][4],3),"]",sep="", fill=TRUE)
  cat(siglevel," HDR CI: [",round(result[[1]][5],3),"," ,round(result[[1]][6],3),"]",sep="", fill=TRUE)
  cat("\n")
  invisible(result)
}

######################################################################################
### Purpose: compute the study specific measure
### Input:   S3 object "multipletables" 
### Output:  study-specific measure(mean, lower bound, upper bound)
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
######################################################################################
study_specifc <- function(object) {
  alpha <- object$alpha
  k <- length(object$sample)
  measure <- object$measure
  model <- object$model
  studynames <- object$studynames
  reports <- matrix(NA,ncol=3,nrow=(k+1))
  out <- list()
   
  # compute study-specific mean and ETCI
    reports[1:k,1] <- sapply(object$sample,mean)   
    reports[1:k,2] <- sapply(object$sample, quantile,probs=alpha/2,na.rm=TRUE)
    reports[1:k,3] <- sapply(object$sample,quantile,probs=1-alpha/2,na.rm=TRUE)
    reports[(k+1),1:3] <- c(object$overall$overall,object$overall$CI[1],object$overall$CI[2])
    rownames(reports) <- c(studynames,"Overall")
    colnames(reports) <- c("Mean","Lower Bound","Upper Bound")  
    out[[1]] <- reports
  
  return(out)
} 
