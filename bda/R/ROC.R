## roc(y,x,z)

roc <- function(y, x, z=NULL, event="less"){
  ny <- length(y);
  nx <- length(x);
  stopifnot(nx == ny);
  if(!is.null(z)){
    stopifnot(length(z)==nx);
  }
  
  x.na <- is.na(x);
  y.na <- is.na(y);
  if(any(x.na|y.na)){
    y <- y[!(x.na|y.na)];
    x <- x[!x.na|y.na];
  }
  if(!is.null(z)){
    z <- z[!x.na|y.na];
  }
  nx <- length(x);

  stopifnot(is.logical(y));
  stopifnot(is.numeric(x));
  
  if(!is.null(z)){
    z.na <- is.na(z)
    if(any(z.na)){
      x <- x[!z.na];
      y <- y[!z.na];
      z <- z[!z.na];
    }
    stopifnot(length(x)>1)
    stopifnot(is.logical(z));
  }else{
    z <- rep(FALSE, nx);
  }

  if(tolower(substr(event,1,1))=='l'){
    xe <- TRUE;
  }else{
    xe <- FALSE;
  }
  
  cutoff <- seq(min(x), max(x), length=102);
  Sensitivity <- NULL
  Specificity <- NULL
  for(i in 1:100){
    yx <- x < cutoff[i+1];
    if(!xe) yx <- !yx;
    yxz <- yx | z;
    TP <- sum(y & yxz);
    TN <- sum(!y & !yxz);
    FN <- sum(y & !yxz);
    FP <- sum(!y & yxz);
    Sensitivity <- c(Sensitivity, TP/(TP+FN));
    Specificity <- c(Specificity, TN/(FP+TN));
  }
  tmp <- diff(Specificity)
  auc <- (Sensitivity + Specificity -1)*c(tmp[1], tmp)
  auc <- sum(auc) + 0.5

  return(structure(list(TPR=Sensitivity, FPR= 1-Specificity,
                        Cutoff=cutoff[2:101], AUC = auc,
                        call = match.call()),
                   class = "broc"))
}

print.broc <- function (x, digits = NULL, ...) 
{
    cat("\nCall:\n\t", deparse(x$call),  
        "\n\t'AUC' = ", formatC(x$AUC, digits = digits),
        "\n\n", sep = "")
    print(as.data.frame(x[c("TPR", "FPR","Cutoff")]),
          digits = digits, ...)
    invisible(x)
}

plot.broc  <-function (x, ...) 
{
    mtitle <- paste("ROC Curve (AUC = ", round(x$AUC,3),")", sep='');
    plot.default(x$TPR~x$FPR, type='l',
                 xlim=c(0,1), ylim=c(0,1),
                 xlab = "FPR or (1-specificity)",
                 ylab = "TPR or sensitivity",
                 main = mtitle,...);
    abline(a=0, b=1, col='red', lty=2);
    invisible(NULL)
}


