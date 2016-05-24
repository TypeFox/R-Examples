btest.mean <-
function(XX,V,theta=1/3,B=100,pic=1){
 #XX...sample (list as always)
 #V...expectation we want to test for
 #theta ... is weight in the def of the bertoluzza metric
 #B...number of bootstrap replicates
 #alpha...confidence level
 k<-length(XX)
 #bind XX and V in list to simplify check for compatibility
 YY<-vector("list",length=(k+1))
 YY[1:k]<-XX[1:k]
 YY[[k+1]]<-V
 
 #it seems easier to simply use the implicit error checking in Msum  
 #to check if everything is in YY  
 temp_sum<-Msum(YY)
 if(nrow(temp_sum)>1){
  nl<-nrow(V)/2
  nobs<-k
  sample_mean <-Mmean(XX,0)
  sample_variance <- Bvar(XX, theta)*(nobs/(nobs-1))
  test_statistic <-bertoluzza(sample_mean,V,theta)^2/sample_variance
  #print(test_statistic)
  if(pic==1){
  lower<-min(sample_mean$x[1],V$x[1])
  upper<-max(sample_mean$x[2*nl],V$x[2*nl])
  limx<-c(min(lower)-0.25,max(upper)+0.25)
  plot(sample_mean,type="l", xlim=limx,lwd=3,xlab=NA, ylab=expression(alpha),cex.main=1, col="black",
           main=paste("Sample mean (in black) and V (in red)",sep=""))
  lines(V,type="l", lwd=3,col="red")  
  }
     
  boot_sample<-replicate(B,sample(XX, nobs,replace=TRUE))
  boot_sample_mean <- apply(boot_sample,2,Mmean)   #list with the means of the sample bootstrap
  boot_sample_variance <- apply(boot_sample,2,Bvar,theta)*(nobs/(nobs-1)) #vector with the sample variances of the sample bootstrap
 
  boot_test_statistic<-rep(0,B)
   for (i in 1:B){
    boot_test_statistic[i] <-bertoluzza(boot_sample_mean[[i]],sample_mean,theta)^2/boot_sample_variance[i]
    #print(boot_test_statistic[i])
   }
  if(pic==1){
    dev.new()
   limx<-c(min(c(boot_test_statistic,test_statistic)),max(c(boot_test_statistic,test_statistic)))
   plot(ecdf(boot_test_statistic),xlab=NA,ylab=NA,xlim=limx, do.points = FALSE, main=paste("Ecdf of T*"),cex.main=1,lwd=1.5)
    #cex.axis=1.3,cex.lab=1.3)
   abline(a = NULL, b = NULL, v = test_statistic,col="red")
   TS<-test_statistic
   mtext(paste("T=",round(TS,2),sep=""), at = TS,  side = 1, line = 2, col = "red", bg="white",cex=1.3)
   }
  pvalue<-mean(test_statistic<boot_test_statistic)
  invisible(pvalue)
 }
}
