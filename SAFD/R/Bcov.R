Bcov <-
function(XX,YY,theta=1/3){
  #calculates the empirical covariance of two lists of polygonal fuzzy numbers with same levels
  #if necessary just use translator first to assure same alpha levels
  #theta ... is weight in the def of the bertoluzza metric
  kx<-length(XX)
  ky<-length(YY)
  if(kx!=ky){
   print("lists have to have same length (i.e. same sample sizes)")
   }
   #------------ calculate integrals by hand as sums -----------------
   int_product<-function(x,y){
    #x,y vector (first column of fuzzy set)
    #calculate integral of the product of x and y (equidistant alpha levels assumed)
    #product of x and y is piecewise quadratic function - integrate via simpson rule
    if(length(x)!=length(y)){return(print("input vectors must have same length"))}
    if(length(x)==length(y)){
     k<-length(x)-1
     delta<-1/k
     pr<-x*y
     middle<-(x[1:k]+x[2:(k+1)])*(y[1:k]+y[2:(k+1)])
     values<-pr[1:k]+pr[2:(k+1)]+middle
     integral<-sum(values)*delta/6
     invisible(integral)
     }
    }
 #----------------------------------------------------------
 #if all ok continue:
 if(kx==ky){
  ZZ<-vector("list",length=(2*kx))
  ZZ[1:kx]<-XX[1:kx]
  ZZ[(kx+1):(2*kx)]<-YY[1:kx]

  temp_sum<-Msum(ZZ)
  if(nrow(temp_sum)>1){
   k<-length(XX)
   EX<-Mmean(XX)
   EY<-Mmean(YY)
   nl<-nrow(XX[[1]])/2

   midEX<-0.5*(EX$x[1:nl]+EX$x[(2*nl):(nl+1)])
   sprEX<-0.5*(EX$x[1:nl]-EX$x[(2*nl):(nl+1)])
   midEY<-0.5*(EY$x[1:nl]+EY$x[(2*nl):(nl+1)])
   sprEY<-0.5*(EY$x[1:nl]-EY$x[(2*nl):(nl+1)])
   contr_Emids<-int_product(midEX,midEY)
   contr_Espreads<-int_product(sprEX,sprEY)

   contr_mids<-rep(0,k)
   contr_spreads<-rep(0,k)

 for(i in 1:k){
   z<-XX[[i]]$x
   w<-YY[[i]]$x
   midX<-0.5*(z[1:nl]+z[(2*nl):(nl+1)])
   sprX<-0.5*(z[1:nl]-z[(2*nl):(nl+1)])
   midY<-0.5*(w[1:nl]+w[(2*nl):(nl+1)])
   sprY<-0.5*(w[1:nl]-w[(2*nl):(nl+1)])
   
   contr_mids[i]<-int_product(midX,midY)/k
   contr_spreads[i]<-int_product(sprX,sprY)/k
   }
   
 cova<-(sum(contr_mids)-contr_Emids + theta*(sum(contr_spreads)-contr_Espreads))
 invisible(cova)
  }
 }
}
