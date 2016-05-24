getR <-
function(checked,resmass,nknots=13,spar=0.1,plotit=TRUE){

    ############################################################################
    # (1) issue warnings #######################################################
	if(length(resmass[,1])<10){stop("stop: not enough data points in resmass\n");}
    if(any(checked[,3]<min(resmass[,1])) || any(checked[,3]>max(resmass[,1]))){stop("stop: some mean_mass out of range of resmass\n");}
    if(nknots<3 || nknots>length(resmass[,1])){stop("stop: invalid nknots\n")}
    if(spar<=0 || spar>=1){stop("stop: invalid spar; spar=(0,1]")}
    if(plotit!="TRUE"&plotit!="FALSE"){stop("stop: plotit invalid. TRUE, FALSE.\n")}
    if(any(checked[,1]==TRUE)){stop("stop: in checked[,1] ... invalid chemical formula!")}
    options(digits=10);
    ############################################################################
    # (2) fit gam model ########################################################
    model<-smooth.spline(resmass[,1],resmass[,2],cv=TRUE,all.knots=FALSE,nknots=nknots,spar=spar)
    if(plotit==TRUE){
      plot(model,pch=19,type="l",lwd=2,col="red",xlab="mass",ylab="Resolution");
      points(resmass[,1],resmass[,2],pch=19,cex=0.7);
    }
    ############################################################################
    # (3) predict ##############################################################
    resolution<-c();
    for(i in 1:length(checked[,1])){
       it<-predict(model,checked[i,3])
       if(plotit==TRUE){
         points(checked[i,3],it$y,pch=19,cex=0.7,col="green");
       }
      resolution<-c(resolution,it$y)  
    }
    return(resolution);
    ############################################################################

}
