bertoluzza <-
function(X,Y,theta=1/3,pic=0){
  #function calculates bertoluzza distance with weigth theta of X and Y 
  #(both polygonal fuzzy numbers, calculated via translator 
  #->same discretization, same matrix dimensions)
  #function returns NA in case of any mistake, otherwise the Bertoluzza distance
  #if necessary just use translator(X,nl), translator(Y,nl) first
  temp_sum<-Msum(list(X,Y))
  if(is.null(temp_sum)==0){
   nl<-nrow(X)/2
   z<-X$x-Y$x
   integrand1<-z[1:nl]+z[(2*nl):(nl+1)]
   integrand2<-z[1:nl]-z[(2*nl):(nl+1)]
   
   #------------ calculate integrals by hand as sums -----------------
   l2dist2<-function(x){
    #x is vector (first column of fuzzy set)
    k<-length(x)-1
    delta<-1/k
    y<-x[1:k]+x[2:(k+1)]
    values<-x[1:k]^2+x[2:(k+1)]^2+y^2
    integral<-sum(values)*delta/6
    invisible(integral)
    }
  
   d1<-(l2dist2(integrand1)+theta*l2dist2(integrand2))*0.25
   distance<-(d1)^0.5
   dis<-round(distance,2)
   th<-round(theta,2)
   
   if(pic==1){
   limx<-c(min(c(X$x,Y$x))-0.25,max(c(X$x,Y$x))+0.25)
    plot(X,type="l", xlim=limx,lwd=2,xlab=NA, ylab=expression(alpha),cex.main=1)
      titletxt <- substitute(D[th] * "=" * dis ,
      list(th = as.character(th),dis=dis))
      title(main=titletxt,cex.main=1)
     lines(Y,type="l",lwd=2)
    }
   invisible(distance)
   }
 }
