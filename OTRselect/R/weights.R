weights <- function(y, delta) {

  func1 <- function(x, time, surv, m){
             if( x <= time[1L] ) {
               Gt <- 1.0
             } else {
               if( x >= time[m] ) {
                 Gt <- surv[m]
               } else {
                 Gt <- surv[{which(x <= time)[1L] - 1L}]
               }
             }
             return(Gt)
           }

  fit <- summary( survival::survfit(survival::Surv(y,1.0-delta)~1))
  if( length(fit$time) < 1L ) {
    stop("Unable to obtain survival fit", call. = FALSE)
  }

  m <- length(fit$time)

  stime <- sapply(X = y, 
                  FUN = func1,  
                  time = fit$time,  
                  surv = fit$surv,  
                  m = m)

  res <- numeric(length(delta))
  res[stime > 0.0] <- sqrt( {delta == 1L} / stime )[stime > 0.0]

  return( res )

}
