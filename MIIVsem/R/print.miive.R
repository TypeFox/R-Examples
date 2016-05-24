#' @method print miive 
#' @export
print.miive <- function(x, digits = 2,...){
  #options(scipen=10, digits=digits)
  cat("\n")
  cat("MIIVsem results \n")
  cat("\n")
  
  vcov     <- x$vcov
  ctrlopts <- x$ctrlopts
  dat   <- x$dat
  model <- x$model  #paste("t", i, )
  modeqns <- x$modeqns
  
  if (is.null(ctrlopts$bootstrap.se)) { `$`(dat , "P(|Z|)") <- round(`$`(dat , "P(|Z|)"),digits) }
  if (!is.null(ctrlopts$bootstrap.se)) { `$`(dat , "P(t)") <- round(`$`(dat , "P(t)"),digits) }
  `$`(dat , "P(Chi)") <- round(`$`(dat , "P(Chi)"),digits)
  dat$DV[duplicated(dat$DV)] <- NA
  cf <- format(dat, digits = digits) ## use format to set other options like digits, justify , ...
  
  if (!is.null(ctrlopts$varcov)){
    if (!is.null(vcov$var)){
      var <- vcov$var
      `$`(var , "P(|Z|)") <- round(`$`(var , "P(|Z|)"),digits)
      var <- format(var, digits = digits)
    }
    if (!is.null(vcov$cov)){
      cov <- vcov$cov
      `$`(cov, "P(|Z|)") <- round(`$`(cov , "P(|Z|)"),digits)
      cov <- format(cov, digits = digits)
    }
  }
  
  suppressWarnings(cf[is.na(dat)] <- "")
  
  if (ctrlopts$covariance == TRUE) { 
    cf <- cf[,c("DV", "EV", "Estimate", "StdErr", "Z", "P(|Z|)")]
  }
  
  print(cf, 
      quote = FALSE, right = FALSE, row.names = FALSE,print.gap = 2, na.print = "")
  
  z_ind <- 1
  for (i in 1:length(model)) {
    if (model[[i]]$NOTE != "") {
      if (z_ind == 1){
          cat("\n")
          cat("Notes:")
          cat("\n")
      }
      cat(model[[i]]$NOTE, "\n")
      z_ind = z_ind + 1
    }
  }
  
  if (!is.null(ctrlopts$varcov)){
    if (!is.null(vcov$cov)){
      cat("\n")
      cat("Conditional covariances \n")
      cat("\n")
      print(cov, quote = FALSE, right = FALSE, row.names = FALSE, print.gap=1)
    }
    if (!is.null(vcov$var)){
      cat("\n")
      cat("Conditional variances \n")
      cat("\n")
      print(var, quote = FALSE, right = FALSE, row.names = FALSE, print.gap=1)
    }
  }
  
  if (!is.null(x$restests)){
    cat("\n")
    cat("Tests of the Linear Restrictions \n")
    cat("\n")
    cat("Restrictions: \n")
    #dft <- data.frame("Restrictions" = cbind(unlist(x$restests[3][[1]])))
    dft <- table(cbind(unlist(x$restests[3][[1]])))
    write.table(dft,row.names=F, col.names=F, quote=F)
    #print.table(dft,quote = FALSE, row.names = FALSE, right=FALSE)
    cat("\n")
    cat(unlist(x$restests[2][[1]][4]))
    cat("\n")
    dfw <- data.frame(as.numeric(x$restests[2][[1]][3]),
                      as.numeric(x$restests[2][[1]][1]),
                      as.numeric(x$restests[2][[1]][2]))
    colnames(dfw) <- c("df", "Chi-sq", "P(Chi)")
    #`$`(dfw, "Chi-sq") <- round(`$`(dfw, "Chi-sq"),digits)
    #`$`(dfw, "P(Chi)") <- round(`$`(dfw, "P(Chi)"),digits)
    dfw <- format(dfw, digits = digits)
    print(dfw,quote = FALSE, row.names = FALSE, right=FALSE, print.gap=1)
    #cat("\n")
    #cat(unlist(x$restests[1][[1]][4]))
    #cat("\n")
    #dfl <- data.frame(as.numeric(x$restests[1][[1]][3]),
    #                  as.numeric(x$restests[1][[1]][1]),
    #                  as.numeric(x$restests[1][[1]][2]))
    #colnames(dfl) <- c("df", "Chi-sq", "P(Chi)")
    ##`$`(dfw, "Chi-sq") <- round(`$`(dfw, "Chi-sq"),digits)
    ##`$`(dfw, "P(Chi)") <- round(`$`(dfw, "P(Chi)"),digits)
    #dfl <- format(dfl, digits = digits)
    #print(dfl,quote = FALSE, row.names = FALSE, right=FALSE, print.gap=1)
  }
  if (ctrlopts$print.miivs == TRUE){
    cat("\n")
    cat("Model Equation Information \n")
    cat("\n")
    print(modeqns, quote = FALSE, right = FALSE, row.names = FALSE, print.gap=1)
  }
}