btestk.mean <-
function(XXX, sel, theta=1/3,B=100,pic=1){
 #XXX ...list of independent samples
 #sel ...selection of the variables to be considered
 #theta ... is weight in the def of the bertoluzza metric
 #B...number of bootstrap replicates
 K<-length(XXX)
 ks<-length(sel)

 #checking
  if(ks>K){
   print("you can not select more variables than the ones contained in the sample XXX")
   }
  if(ks<=1){
   print("you have to select at least two variables (in XXX)")
   }
 if(ks<=K&ks>1){
  #checking samples in XXX are compatible
  YYY<-vector("list",length=ks)
  nobs<-rep(0,ks)
  sel<-sort(sel)
  for(i in 1:ks){
    YYY[[i]]<-XXX[[sel[i]]]
    nobs[i]<-length(YYY[[i]])
   }
  ZZ<-vector("list",length=sum(nobs))
  ZZ[1:nobs[1]]<-YYY[[1]][1:nobs[1]]
  selsum<-cumsum(nobs)
  for(i in 1:(ks-1)){
    ZZ[(selsum[i]+1):(selsum[i+1])] <-YYY[[i+1]][1:nobs[i+1]]
  }
 #checking done  
 temp_sum<-Msum(ZZ)
 if(nrow(temp_sum)>1){
   nl<-nrow(temp_sum)/2
  #compute the test statistic
   sample_mean<-vector("list",length=ks)
   sample_sum<-vector("list",length=ks)
   total_mean <-vector("list",length=ks)
   sample_variance<-rep(0,ks)

   for (i in 1:ks){
 	  sample_mean[[i]]<-Mmean(YYY[[i]],0)
	  sample_sum[[i]]<-sc_mult(sample_mean[[i]],nobs[i])
	  sample_variance[i]<-Bvar(YYY[[i]],theta)
   }
  total_mean <- sc_mult(Msum(sample_sum),1/sum(nobs))
 
  #optional plotting of the group means
  #in case of at most 10 groups a legend is plotted
  if(pic==1){
    lower<-sample_mean[[1]]$x[1]
    upper<-sample_mean[[1]]$x[2*nl]
    for (i in 2:ks){
       lower<-min(lower,sample_mean[[i]]$x[1])
       upper<-max(upper,sample_mean[[i]]$x[2*nl])
    }
   legend_name<-paste(rep("group ",ks),sel,sep="")
   limx<-c(lower,upper)+c(0,(upper-lower)/4)
   color<- colorRampPalette( c("green","blue","red"))(ks)
   plot(total_mean,type="l", xlim=limx,lwd=3,xlab=NA, ylab=expression(alpha), col="black",
          main=paste("Total mean (black) and group means","\n", 
          "(group mean colour ranging from green to blue to red)",sep=""),
          cex.main=1)
   for (i in 1:ks){
   lines(sample_mean[[i]],type="l", lwd=2,col=color[i])
   }
   if(ks<=10){
    legend(upper, 1, legend_name, col = color, text.col = "black", lty = rep(1,ks),cex=0.8)
   }
  }

  total_variance <- sum(sample_variance)

  temp<-rep(0,ks)
  for (i in 1:ks){
   temp[i]<-bertoluzza(sample_mean[[i]],total_mean,theta)^2
  }
  test_statistic <-sum(nobs*temp)/total_variance
  #print(test_statistic)
 #sample under H0
 samplestar<-list()
  for (i in 1:ks){
  samplestar[[i]]<-list()
  #calculate Mmeans of all groups different to the i-th and their sum
  relevant<-setdiff(seq(1,ks,by=1),i)
  Mean_list<-vector("list",length=length(relevant))
  for (m in 1:length(relevant)){
    Mean_list[[m]]<-Mmean(YYY[[relevant[m]]])
  }
  suplement<-Msum(Mean_list)
  #add suplement to all observations of variable i
  for (j in 1:nobs[i]){
   samplestar[[i]][[j]] <- Msum(list(YYY[[i]][[j]],suplement))
  }
  }

  ########bootstraping

  boot_sample<-list()
  boot_sample_mean<-list()
  boot_sample_sum<-list()
  boot_sample_variance<-list()
  boot_total_mean<-list()
  boot_total_variance<-rep(0,B)
  boot_test_statistic<-rep(0,B)

	for (b in 1:B){
	 print(b)
	 boot_sample[[b]]<-list()
	 boot_sample_mean[[b]]<-list()
	 boot_sample_sum[[b]]<-list()
	 boot_sample_variance[[b]]<-rep(0,ks)

			for (i in 1:ks){
			boot_sample[[b]][[i]]<-vector("list",length=nobs[i])
			boot_sample[[b]][[i]] <- sample(samplestar[[i]], nobs[i],replace=TRUE) #sample bootstrap of Xi
			boot_sample_mean[[b]][[i]] <- Mmean(boot_sample[[b]][[i]])	
			boot_sample_sum[[b]][[i]] <- sc_mult(boot_sample_mean[[b]][[i]],nobs[i])
			boot_sample_variance[[b]][[i]] <- Bvar(boot_sample[[b]][[i]],theta)	
			}
  
   boot_total_mean[[b]] <- sc_mult(Msum(boot_sample_sum[[b]]),1/sum(nobs))
	 boot_total_variance[[b]] <- sum(boot_sample_variance[[b]])
	 temp<-rep(0,ks)
	  for (i in 1:ks){
	   temp[i]<-bertoluzza(boot_sample_mean[[b]][[i]],boot_total_mean[[b]],theta)^2
    }
	 boot_test_statistic[b] <- sum(nobs*temp)/boot_total_variance[b]
  }
  if(pic==1){
   dev.new()
   limx<-c(min(c(boot_test_statistic,test_statistic)),max(c(boot_test_statistic,test_statistic)))
      plot(ecdf(boot_test_statistic),xlab=NA,ylab=NA,xlim=limx, do.points = FALSE, main=paste("Ecdf of T*"),cex.main=1,lwd=1.5)
    #cex.axis=1.3,cex.lab=1.3)
   abline(a = NULL, b = NULL, v = test_statistic,,col="red")
   TS<-test_statistic
   mtext(paste("T=",round(TS,2),sep=""), at = TS,  side = 1, line = 2, col = "red", bg="white",cex=1.3)
  }
  pvalue<-mean(test_statistic<boot_test_statistic)
  invisible(pvalue)
 }
 }
}
