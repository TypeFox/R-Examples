#$Id: print.R,v 1.4 2004/02/09 08:08:21 peters Exp $

print.classbagg <- function(x, digits=4, ...)
{
    cat("\n")
    B <- length(x$mtrees)
    if (B > 1)
      method <- paste("Bagging classification trees with", B, 
                      "bootstrap replications")
    else 
      method <- "Classification tree"
    cat(method, "\n")
    if (!is.null(x$call)) {
      cat("\nCall: ")
      print(x$call)
      cat("\n")
    }
    if (x$OOB) {
      cat("Out-of-bag estimate of misclassification error: ",
           round(x$err, digits), "\n")
    }
    cat("\n")
}

print.regbagg <- function(x, digits=4, ...)
{
    cat("\n")
    B <- length(x$mtrees)
    if (B > 1)
      method <- paste("Bagging regression trees with", B, 
                    "bootstrap replications")
    else
      method <- "Regression tree"
    cat(method, "\n")
    if (!is.null(x$call)) {
      cat("\nCall: ")
      print(x$call)
      cat("\n")
    }
    if (x$OOB)
      cat("Out-of-bag estimate of root mean squared error: ",
           round(x$err, digits), "\n")
    cat("\n")

}

print.survbagg <- function(x, digits=4, ...)
{
    cat("\n")
    B <- length(x$mtrees)
    if (B > 1)
      method <- paste("Bagging survival trees with", B, 
                      "bootstrap replications")
    else
      method <- "Survival tree"
    cat(method, "\n")
    if (!is.null(x$call)) {
      cat("\nCall: ")
      print(x$call)
      cat("\n")
    }
    if (x$OOB)
      cat("Out-of-bag estimate of Brier's score: ",
           round(x$err, digits), "\n")
    cat("\n")

}

summary.classbagg <- function(object, ...)
{
     print(object, ...)
     class(object) <- "summary.bagging"
     object
}

summary.regbagg <- function(object, ...)
{
     print(object, ...)
     class(object) <- "summary.bagging"
     object
}

summary.survbagg <- function(object, ...)
{
     print(object, ...)
     class(object) <- "summary.bagging"
     object
}

print.summary.bagging <- function(x, digits = max(3, getOption("digits")-3),
                                 ...)
{
     cat("Trees: \n")
     print(x$mtrees)
     invisible(x$mtrees)
}

print.cvclass <- function(x, digits=4, ...)
{
  cat("\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("\t", paste(x$k, "-fold", sep=""), 
      "cross-validation estimator of misclassification error \n") 
  cat("\n")
  cat("Misclassification error: ", round(x$error, digits), "\n")
  cat("\n")
}

print.bootestclass <- function(x, digits=4, ...) {
  if(all(names(x)[names(x)!="call"] %in% c("boot", "632plus"))) {
    XX <- x
    for(i in c("boot", "632plus")) {
      x <- XX[[i]]
      x$call <- XX[["call"]]
      cat("\n")
      if (!is.null(x$call)) {
        cat("Call:\n")
        print(x$call)
        cat("\n")
      }
      if (x$bc632plus) {
        cat("\t", ".632+ Bootstrap estimator of misclassification error \n")
      } else { 
        cat("\t", "Bootstrap estimator of misclassification error \n")
      }
      cat("\t with" , x$nboot, "bootstrap replications\n")
      cat("\n")
      cat("Misclassification error: ", round(x$error, digits), "\n")
      if (!x$bc632plus) cat("Standard deviation:", round(x$sd, digits), "\n")
      cat("\n")
    }
  } else {
# if(!all(names(x) %in% c("boot", "632plus"))){
    cat("\n")
    if (!is.null(x$call)) {
      cat("Call:\n")
      print(x$call)
      cat("\n")
    }
    if (x$bc632plus) 
      cat("\t", ".632+ Bootstrap estimator of misclassification error \n")
    else 
      cat("\t", "Bootstrap estimator of misclassification error \n")
    cat("\t with" , x$nboot, "bootstrap replications\n")
    cat("\n")
    cat("Misclassification error: ", round(x$error, digits), "\n")
    if (!x$bc632plus)
      cat("Standard deviation:", round(x$sd, digits), "\n")
    cat("\n")   
  }
}
  


print.cvreg <- function(x, digits=4, ...)
{
  cat("\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call) 
    cat("\n")
  }
  cat("\t", paste(x$k, "-fold", sep=""),
      "cross-validation estimator of root mean squared error\n")
  cat("\n")
  cat("Root mean squared error: ", round(x$error, digits), "\n")
  cat("\n")
}

print.bootestreg <- function(x, digits=4, ...)
{
  cat("\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("\t", "Bootstrap estimator of root mean squared error \n")
  cat("\t with" , x$nboot, "bootstrap replications\n")
  cat("\n")
  cat("Root mean squared error: ", round(x$error, digits), "\n")  
  cat("\n")
}


print.cvsurv <- function(x, digits=4, ...)
{
  cat("\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call) 
    cat("\n")
  }
  cat("\t", paste(x$k, "-fold", sep=""),
      "cross-validation estimator of Brier's score\n")
  cat("\n")
  cat("Brier's score: ", round(x$error, digits), "\n")
  cat("\n")
}

print.bootestsurv <- function(x, digits=4, ...)
{
  cat("\n")
  if (!is.null(x$call)) {
    cat("Call:\n")
    print(x$call)
    cat("\n")
  }
  cat("\t", "Bootstrap estimator of Brier's score\n")
  cat("\t with" , x$nboot, "bootstrap replications\n")
  cat("\n")
  cat("Brier's score: ", round(x$error, digits), "\n")  
  cat("\n")
}
