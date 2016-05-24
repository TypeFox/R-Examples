"rhoCOP" <-
function(cop=NULL, para=NULL, method=c("default", "joe21", "joe12"), as.sample=FALSE,
                              brute=FALSE, delta=0.002, ...) {

   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Spearman's Rho desired but para is NULL, ",
                 "returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, ",
                "returning NULL")
        return(NULL)
      }
      return(cor(para[,1], para[,2], method="spearman"))
   }


   method <- match.arg(method)
   if(brute) {
      # Following logic would implement the concordance function via tauCOP
      # if the previous if(brute) was not there, but here we only resort to if
      # brute is request.
      Q <- tauCOP(cop=cop, cop2=P, para=para, brute=brute, delta=delta, ...);
      rho <- 3*Q;
      return(rho)
   }
   myint <- NULL
   if(method == "default") {
      try(myint <- integrate(function(u) {
            sapply(u,function(u) { integrate(function(v) {
            COP( u, v, cop=cop, para=para, ...) - u*v}, 0, 1)$value })}, 0, 1) )
      rho1 <- ifelse(is.null(myint), NA, 12*myint$value)
      # This version seems the most stable given M and W copulas
      return(rho1)
   } else if(method == "joe21") {
      try(myint <- integrate(function(u) {
            sapply(u,function(u) { integrate(function(v) {
            u * derCOP( u, v, cop=cop, para=para, ...)},  0, 1)$value })}, 0, 1))
      ifelse(is.null(myint), return(NA), return(3 - 12*myint$value))
   } else if(method == "joe12") {
      try(myint <- integrate(function(u) {
            sapply(u,function(u) { integrate(function(v) {
            v * derCOP2( u, v, cop=cop, para=para, ...)}, 0, 1)$value })}, 0, 1))
      ifelse(is.null(myint), return(NA), return(3 - 12*myint$value))
   } else {
      stop("Should not be here in logic, EVER")
   }
}

