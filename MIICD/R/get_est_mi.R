
get_est_mi<-function( x , status , trans , imp_sets , data = data , cens.code = cens.code , model = c('Cox','FG') , r2 ){

  t1 <- imp_sets[ , x ]
  fitCI <- survfit( Surv( time = t1 , event = r2 , type = "mstate"  ) ~ 1 )
  w <- which( fitCI$states == trans )
  pr <- fitCI$prev[ , w ]
  sd <- fitCI$std.err[ , w ]
  t1 <- fitCI$time
  CI <- list(time = t1 , est = pr , sd = sd )
  return( CI )
}


