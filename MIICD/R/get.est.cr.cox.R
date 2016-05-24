
get.est.cr.cox <-
function( x , data , status , trans , cens.code , formula ){
  keep <- as.vector(sapply(strsplit( as.character(formula)[2] , '\\+' ) , function(x) gsub(" " , "", x) ))
  x[x==0]<-1e-16
  df1 <- data.frame( data , 'ns' = x )
  df2 <- mstate::crprep( data = df1 , Tstop = x , status = status , trans = trans , cens = cens.code , keep = keep )
  s1 <- with( df2 , survival::Surv( time = Tstart , time2 = Tstop , type = 'counting' , event = status==trans ) )
  formula1 <- formula( paste( 's1' , deparse( formula ) ) )
  cph1 <- survival::coxph( formula1 , data = df2 , weights = df2$weight.cens , ties=c("efron") )
vc1 <- vcov( cph1 )
coef1 <- cph1$coefficients
return( list( sigma = vc1 , beta = coef1) )
}


