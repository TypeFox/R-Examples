#imp_sets <- sets
#x<-1
get_est_mi_surv<-function( x , imp_sets , data = data  ){

  t1 <- imp_sets[ , x ]
  r3 <- data$right != Inf
  surv <-  Surv( time = t1 , event = r3 , type = "right"  )
  surv2 <- Surv( time = t1 , event = r3 , type = "mstate"  )
  surv2[,2] <- surv[,2]
  fitCI <- survfit( surv2 ~ 1 , conf.type = 'none') 
  sd <- fitCI$std.err
  pr <- fitCI$prev
  t0 <- fitCI$time
  CI <- list(time = t0 , est = pr , sd = sd )
  return( CI )
}
