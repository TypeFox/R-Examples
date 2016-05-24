"concordCOP" <- function(cop=NULL, para=NULL, cop2=NULL, para2=NULL, ...) {
   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }
   if(is.null(cop2)) {
      warning("must provide the second copula [and optional parameters]")
      return()
   }
   tauCOP(cop=cop, para=para, cop2=cop2, para2=para2, ...)
}

"tauCOP" <-
function(cop=NULL,  para=NULL,
         cop2=NULL, para2=NULL, as.sample=FALSE,  brute=FALSE, delta=0.002, ...) {
   if(as.sample) {
      if(is.null(para)) {
         warning("Sample Kendall's Tau desired but para is NULL, ",
                 "returning NULL")
         return(NULL)
      }
      if(length(names(para)) != 2) {
        warning("para argument must be data.frame having only two columns, ",
                "returning NULL")
        return(NULL)
      }
      return(cor(para[,1], para[,2], method="kendall"))
   }

   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }

   if(is.null(cop2)) {
     cop2  <- cop  # This is the expected operation
     para2 <- para # as Kendall's Tau gets returned
   }

   us <- vs <- seq(.Machine$double.eps, 1-.Machine$double.eps, delta)
   if(brute) {
     Q <- sum(sapply(us, function(u) {
             sum(sapply(vs, function(v) {
                return( derCOP(u,v, cop=cop,  para=para,  ...) *
                        derCOP2(u,v, cop=cop2, para=para2, ...))
             }))
          }))
     Q <- 4*(0.5 - Q*delta^2) - 1 # SEE P 164 of NELSON 2006!
     if(Q >  1) Q <-  1 # assume rounding errors just breaking through
     if(Q < -1) Q <- -1 # again for rounding errors
     return(Q)
   }

   myint <- NULL
   try(myint <- integrate(function(u) {
               sapply(u,function(u) {
                 integrate(function(v) {
                  derCOP( u,v, cop=cop,  para=para,  ...) *
                 derCOP2( u,v, cop=cop2, para=para2, ...)}, 0,1)$value
             })}, 0,1))
    if(is.null(myint)) {
        warning("error on dual partial derivative integration ",
                "(copula singularity?), swapping copulas (Nelsen corollary ",
                "5.1.2)")
        try(myint <- integrate(function(u) {
               sapply(u,function(u) {
                 integrate(function(v) {
                  derCOP( u,v, cop=cop2, para=para2, ...)*
                 derCOP2( u,v, cop=cop,  para=para,  ...)}, 0,1)$value
             })}, 0,1))
        if(is.null(myint)) {
           warning("another error on dual partial derivative integration ",
                   "(again copula singularity?), kicking over to integration ",
                   "of the Kendall Function")
           tau<-3-4*integrate(function(f) kfuncCOP(f,cop=cop,para=para),0,1)$value
           return(tau)
        }
    }
    Q <- 4*(0.5 - myint$value) - 1
    if(Q >  1) Q <-  1 # assume rounding errors just breaking through
    if(Q < -1) Q <- -1 # again for rounding errors
    return(Q)
}

