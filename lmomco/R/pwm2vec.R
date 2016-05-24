"pwm2vec" <-
function(pwm,...) {
   if(length(pwm$betas) == 0) {
      warning("argument is not a pwm object (no betas field)")
      return(NULL)
   }
   return(pwm$betas)
}
