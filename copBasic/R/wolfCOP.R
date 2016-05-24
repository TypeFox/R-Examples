"wolfCOP" <-
function(cop=NULL, para=NULL, as.sample=FALSE, brute=FALSE, delta=0.002, ...) {

    if(as.sample) {
      if(is.null(para)) {
         warning("Sample Schweizer and Wolff's Sigma desired but para is NULL, ",
                 "returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, ",
                "returning NULL")
        return(NULL)
      }
      if(as.sample == -1) message("Sample Schweizer and Wolff's Sigma",
                                  "---CPU intensive!")
      # http://www.cs.cmu.edu/~bapoczos/articles/poczos11nipscopula.pdf
      #                                               (August 11, 2015)
      n <- length(para[,1]); nn <- n^2
      R <- rank(para[,1]); S <- rank(para[,2])
      samSIG <- sum(sapply(1:n, function(i) {
             sum(sapply(1:n, function(j) {
                   abs((sum(as.numeric(R <= i & S <= j))/n) - (i*j/nn))
             } ))
           } ))
      samSIG <- (12/(nn - 1)) * samSIG
      return(samSIG)
   }

   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }

   if(brute) {
      us <- vs <- seq(.Machine$double.eps, 1-.Machine$double.eps, delta)
      wolf <- sum(sapply(us, function(u) {
                 sum(sapply(vs, function(v) {
                    abs(cop(u,v, para=para, ...) - u*v) })) }))
      return(12*wolf*delta^2)
   }

   myint <- NULL
   try(myint <- integrate(function(u) {
            sapply(u,function(u) {
                      integrate(function(v) {
                          abs(COP( u, v, cop=cop, para=para, ...) - u*v)},
                      0, 1)$value
            })}, 0, 1) )
   wolf <- ifelse(is.null(myint), NA, 12*myint$value)
   return(wolf)
}

