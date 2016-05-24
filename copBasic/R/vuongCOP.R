"vuongCOP" <- function(u, v=NULL, cop1=NULL, cop2=NULL, para1=NULL, para2=NULL,
                                  alpha=0.05, method=c("D12", "AIC", "BIC"), ...) {
   method <- match.arg(method)
   if(is.null(v)) {
      if(length(names(u)) != 2) {
         warning("a data.frame having only two columns is required")
         return(NULL)
      } else {
         v <- u[,2]
         u <- u[,1]
      }
   }
   if(length(u) != length(v)) {
      warning("argument(s) or implied arguments u and v are unequal in length, returning NULL")
      return(NULL)
   }

   n <- length(u)
   the.qt <- qt(1 - alpha/2, df=(n-2))

   "Di" <- function(x,y) {
       f1 <- densityCOP(x,y, cop=cop1, para=para1, ...)
       f2 <- densityCOP(x,y, cop=cop2, para=para2, ...)
       tmp <- log(f2/f1)
       return(tmp)
   }

   hatD12 <- sum(sapply(1:n, function(i) {  Di(u[i],v[i])             })) / n
   varD12 <- sum(sapply(1:n, function(i) { (Di(u[i],v[i]) - hatD12)^2 })) / (n-1)
   sigmaD12 <- sqrt(varD12)

   parameters <- c(alpha, as.integer(n), the.qt, sigmaD12)
   names(parameters) <- c("alpha", "sample.size", "qt(1-alpha/2, df)", "sigmaD12")

      lo  <-  hi  <- the.qt*sigmaD12/sqrt(n)
   D12lo  <- hatD12 - lo; D12hi <- hatD12 + hi
   D12    <- c(VuongD.lower=D12lo, VuongD=hatD12, VuongD.upper=D12hi,
               conf.interval=(100*(1-alpha)))
   names(D12) <- c("VuongD.lower", "VuongD", "VuongD.upper", "conf.interval")
   D12pvalue <- pt(hatD12*sqrt(n)/sigmaD12, df=(n-2))
   if(hatD12 == 0) {
     D12pvalue <- 1 # median
   } else if(hatD12 < 0) {
     D12pvalue <-    D12pvalue  * 2
   } else {
     D12pvalue <- (1-D12pvalue) * 2
   }

   dim1 <- dim(para1)
   dim1 <- ifelse(is.null(dim1), length(para1), sum(dim1))
   dim2 <- dim(para2)
   dim2 <- ifelse(is.null(dim2), length(para2), sum(dim2))
   # AIC and BIC will match if the dimension of the parameter vector/matrix
   # for the two copulas equal each other and thus the difference is zero.
   AICm <- hatD12 -              (dim2 - dim1) / n
   AIClo <- AICm - lo; AIChi <- AICm + hi
   AIC <- c(AIC.lower=AIClo, AIC=AICm, AIC.upper=AIChi,
                     conf.interval=(100*(1-alpha)))
   names(AIC) <- c("AIC.lower", "AIC", "AIC.upper", "conf.interval")
   AICpvalue <- pt(AICm*sqrt(n)/sigmaD12, df=(n-2))
   if(AICm == 0) {
     AICpvalue <- 1 # median
   } else if(hatD12 < 0) {
     AICpvalue <-    AICpvalue  * 2
   } else {
     AICpvalue <- (1-AICpvalue) * 2
   }

   BICm <- hatD12 - (1/2)*log(n)*(dim2 - dim1) / n
   BIClo <- BICm - lo; BIChi <- BICm + hi
   BIC <- c(BIC.lower=BIClo, BIC=BICm, BIC.upper=BIChi,
                     conf.interval=(100*(1-alpha)))
   names(BIC) <- c("BIC.lower", "BIC", "BIC.upper", "conf.interval")
   BICpvalue <- pt(BICm*sqrt(n)/sigmaD12, df=(n-2))
   if(BICm == 0) {
     BICpvalue <- 1 # median
   } else if(hatD12 < 0) {
     BICpvalue <-    BICpvalue  * 2
   } else {
     BICpvalue <- (1-BICpvalue) * 2
   }

   if(method == "D12") {
      CKlo <- D12lo; CKhi <- D12hi
      txt <- "The D12 method used for relative fit judgement (not fit adequacy)"
   } else if(method == "AIC") {
      CKlo <- AIClo; CKhi <- AIChi
      txt <- "The Akaike (AIC) method used for relative fit judgement (not fit adequacy)"
   } else {
      CKlo <- BIClo; CKhi <- BIChi
      txt <- "The Schwarz (BIC) method used for relative fit judgement (not fit adequacy)"
   }
   if(CKlo < 0 & CKhi > 0) { # The interval contains zero
      message <- "Copulas 1 and 2 are not significantly different at 100(1-alpha)"
      result  <- 0
   } else if(CKhi < 0) {
      message <- "Copula 1 has better fit than Copula 2 at 100x(1-alpha) level"
      result  <- 1
   } else if(CKlo > 0) {
      message <- "Copula 2 has better fit than Copula 1 at 100x(1-alpha) level"
      result  <- 2
   } else {
      message <- "Error: Evidently Copula 1 is the same as Copula 2"
      result  <- NA
   }

   pvalues <- c(D12pvalue, AICpvalue, BICpvalue)
   names(pvalues) <- c("VuongD p-value (two side)",
                          "AIC p-value (two side)",
                          "BIC p-value (two side)")

   zz <- list(title="Vuong's Procedure for Parametric Copula Comparison",
              method=txt, message=message, result=result, p.value=pvalues,
              D12=D12, AIC=AIC, BIC=BIC, parameters=parameters)
   return(zz)
}
