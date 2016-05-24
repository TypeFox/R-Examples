"copBasic.fitpara.beta" <-
function(uv=NULL, popstat=NULL, statf=NULL, cop=NULL,
         paradim=1, interval=NULL, par.init=NULL, ...) {
   if(is.null(cop)) {
      warning("must have copula argument specified, returning NULL")
      return(NULL)
   }
   if(paradim == 1) {
      if(is.null(interval)) {
         warning("must have parameter interval specified, returning NULL")
         return(NULL)
      }
   } else {
      if(is.null(par.init)) {
         warning("must have initial guesses specified, returning NULL")
         return(NULL)
      }
   }
   if(is.null(statf)) {
      warning("must have measure of association function specified, ",
              "returning NULL")
      return(NULL)
   }

   #target.stat <- NULL
   ifelse(is.null(uv), target.stat <- popstat,
                       target.stat <- statf(para=uv, as.sample=TRUE))

   if(paradim == 1) {
      "objfuni" <- function(p, stat=NULL, cop=NULL, ...) {
         stat - statf(cop=cop, para=p, ...)
      }
   } else {
      "objfmulti" <- function(p, cop=cop, ...) statf(p, ...)
   }

   rt <- NULL
   if(paradim == 1) {
      try(rt <- uniroot(objfuni, interval, stat=target.stat, cop=cop))
   } else {
      try(rt <- optim(par.init, objfmulti, stat=target.stat, cop=cop))
   }
   if(is.null(rt)) {
      warning("could not find parameter(s) numerically, returning NA")
      return(NA)
   }

   ifelse(paradim > 1, para <- rt$par, para <- rt$root)
   return(para)
}
