"rhobevCOP" <-
function(cop=NULL, para=NULL, as.sample=FALSE,
                              brute=FALSE, delta=0.002, ...) {

   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Bivariate Rho desired but para is NULL, ",
                 "returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, ",
                "returning NULL")
        return(NULL)
      }

      return(mean(-log(para[,1]) * -log(para[,2])) - 1 )
   }

   if(brute) {
      ws <- seq(.Machine$double.eps, 1-.Machine$double.eps, delta)
      sum <- sum(sapply(ws, function(w) {
              (-log(cop(exp(-w),exp(w-1), para=para, ...)))^(-2) }))
      rhoE <- sum*delta - 1
      return(rhoE)
   }

   myint <- NULL
   try(myint <- integrate(function(w) {
                  (-log(cop(exp(-w),exp(w-1), para=para, ...)))^(-2) }, 0, 1) )
   rhoE <- ifelse(is.null(myint), NA, myint$value - 1)
   return(rhoE)
}

