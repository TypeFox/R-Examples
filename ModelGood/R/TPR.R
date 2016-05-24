# {{{ UseMethod

TPR <- function(x,...){
    UseMethod("TPR",x)
}

TPR.numeric <- function(x,...){
    tmp <- binom.test(x[4],(x[3]+x[4]),...)
    ci.tpr <- c(tmp$estimate,tmp$conf.int)
    names(ci.tpr) <- c("TPR (Sensitivity)","Lower","Upper")
    ci.tpr
}

TPR.table <- function(x,...){
  tmp <- binom.test(x[4],(x[3]+x[4]),...)
  ci.tpr <- c(tmp$estimate,tmp$conf.int)
  names(ci.tpr) <- c("TPR (Sensitivity)","Lower","Upper")
  ci.tpr
}

TNR <- function(x,...){
  UseMethod("TNR",x)
}

TNR.numeric <- function(x,...){
    tmp <- binom.test(x[1],(x[1]+x[2]),...)
    ci.tnr <- c(tmp$estimate,tmp$conf.int)
    names(ci.tnr) <- c("TPR (Specificity)","Lower","Upper")
    ci.tnr
}
TNR.table <- function(x,...){
    tmp <- binom.test(x[1],(x[1]+x[2]),...)
    ci.tnr <- c(tmp$estimate,tmp$conf.int)
    names(ci.tnr) <- c("TPR (Specificity)","Lower","Upper")
    ci.tnr
}

# }}}
