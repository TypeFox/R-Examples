btest2.mean <-
function(XX,YY,theta=1/3,B=100,pic=1){
 #XX, YY ... independent samples
 #theta ... is weight in the def of the bertoluzza metric
 #B...number of bootstrap replicates
 kx<-length(XX)
 ky<-length(YY)
 #bind XX and YY in list to simplify check for compatibility
 ZZ<-vector("list",length=(kx+ky))
 ZZ[1:kx]<-XX[1:kx]
 ZZ[(kx+1):(kx+ky)]<-YY[1:ky]
 
 temp_sum<-Msum(ZZ)
 if(nrow(temp_sum)>1){
  n1obs<-kx
  n2obs<-ky
 
  nl<-nrow(XX[[1]])/2
  #compute the test statistic
  sample_mean_XX <-Mmean(XX,0)
  sample_mean_YY <-Mmean(YY,0)
  sample_variance <- Bvar(XX, theta)/(n1obs-1) + Bvar(YY, theta)/(n2obs-1)
  test_statistic <-bertoluzza(sample_mean_XX,sample_mean_YY,theta)^2/sample_variance
  #print(test_statistic)
  if(pic==1){
   lower<-min(sample_mean_XX$x[1],sample_mean_YY$x[1])
   upper<-max(sample_mean_XX$x[2*nl],sample_mean_YY$x[2*nl])
   limx<-c(min(lower)-0.25,max(upper)+0.25)
   plot(sample_mean_XX,type="l", xlim=limx,lwd=2,xlab=NA, ylab=expression(alpha),cex.main=1, col="black",

          main=paste("Sample mean 1st sample (in black) and 2nd sample (in red)",sep=""))
   lines(sample_mean_YY,type="l", lwd=2,col="red")
   }
  XXstar<-vector("list",length=n1obs)
   for (i in 1:n1obs){
     XXstar[[i]]<-Msum(list(XX[[i]],sample_mean_YY))
   }
  YYstar<-vector("list",length=n2obs)
   for (i in 1:n2obs){
     YYstar[[i]]<-Msum(list(YY[[i]],sample_mean_XX))
   }

  #####bootstrap technique to test the mean equality

  #sample bootstrap of XX
  boot_sample_XX<-replicate(B,sample(XXstar, n1obs,replace=TRUE))
  #list with the means of the sample bootstrap of XX
  boot_sample_mean_XX <- apply(boot_sample_XX,2,Mmean)
  #vector with the sample variances of the sample bootstrap of XX
  boot_sample_variance_XX <- apply(boot_sample_XX,2,Bvar,theta)

  #sample bootstrap of YY
  boot_sample_YY<-replicate(B,sample(YYstar, n2obs,replace=TRUE))
  #list with the means of the sample bootstrap of YY
  boot_sample_mean_YY <- apply(boot_sample_YY,2,Mmean)
  #vector with the sample variances of the sample bootstrap of YY
  boot_sample_variance_YY <- apply(boot_sample_YY,2,Bvar,theta)

  boot_test_statistic<-rep(0,B)
   for (i in 1:B){
	   boot_test_statistic[i] <-bertoluzza(boot_sample_mean_XX[[i]],boot_sample_mean_YY[[i]],theta)^2/(boot_sample_variance_XX[i]/(n1obs-1)+boot_sample_variance_YY[i]/(n2obs-1))
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
