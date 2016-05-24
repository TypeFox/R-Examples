
post_point_est_CI<-function(  beta , sd , times , conf.int = F , alpha = 0.05 ){

m <- ncol(beta)
#mean of point estimates
mpe <- apply( beta , 1 , mean )
  
#mean of variances
if(conf.int) {
  #
#WIV <- apply( sigma^2 , 1 , mean )
WIV <- sd^2
#Between imputation variance
#Matrix of point estimates at single times
M <- matrix( rep( mpe , m ) , ncol = m )

#between-imputation variance
BIV <- apply(beta,1,var)
#With inflation factor
IBIV <- ( 1 + ( 1/m ) ) * BIV
#log-log transformation
z <- qnorm( 1 - alpha / 2 )
sigma <- ( WIV + IBIV )**.5
#Confidance interval

Lt1 <- mpe + z * sigma
Lt0 <- mpe^( exp( (  + z*sigma / (mpe * abs(log( mpe ) ) ) ) ) )
Lt<-apply(data.frame(Lt1,Lt0),1,min)
Ut0 <- mpe^( exp( (  - z*sigma / (mpe * abs(log( mpe ) ) ) ) ) )
Ut1 <- mpe + z * sigma
Ut<-apply(data.frame(Ut1,Ut0),1,min)
  
  
  #Get results in a data frame
CI<-data.frame( time = times , est = mpe ,  sd = sigma , lci = Lt , uci = Ut )
l1<-c(0,0,0,0,0)
#limit of times
#dl<-tail(CI,1)
#CI<-rbind(l1,CI,dl)
CI<-rbind(l1,CI)  
return(CI)
}else{
CI<-data.frame( time = times , est = mpe )
l1<-c(0,0)
#limit of times
#dl<-tail(CI,1)
#CI<-rbind(l1,CI,dl)
CI<-rbind(l1,CI)  
return(CI)
}
}





  