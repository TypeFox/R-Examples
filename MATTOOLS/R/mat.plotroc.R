 

 
mat.plotroc<-function (inRocObj) 
{
#$> names(inRocObj[[1]])
#[1] "zonezone" "outzone"  "ROC"     
#> names(mat.roc.out3[[1]]$ROC)
#[1] "TNF"         "TPF"         "FPE"         "slope"       "threasholds"
#[6] "optDissVal"  "maxPos"      "AUC"         "SEAUC"  

devn=vector("numeric")
classnames=names(inRocObj)
numzones=length(inRocObj)-1
ndevices=ceiling(numzones/4)
devcnt=1
for(k in 1:ndevices){
        windows(width = 8.5, height = 11, pointsize = 10)
        devn=c(devn,dev.cur())
        thelay=layout(matrix(1:12,ncol =3,nrow=4,byrow = TRUE))
        layout.show(thelay)
  }
dev.set(devn[1])

for(i in 1:numzones){
        
#  1. Get the current list from the ROC MATTOOLS object

        list1=inRocObj[[i]]

#  2. Get Optimal value of SCD and store FPE and TPF
 
        maxslopeloc=list1$ROC$maxPos[1]
        maxslope = list1$ROC$optDissVal[1]

#  3. Plot Histograms for SCD values within zone and superimpose on these the histogram outside of zone.  Calculate relative frequencies so that sum of bar tops sums to 1.

        hzone=hist(list1$zonezone,breaks=seq(0,2,0.05),freq=F,xlim=c(0,2),plot=F)
        houtzone=hist(list1$outzone,breaks=seq(0,2,0.05),freq=F,xlim=c(0,2),plot=F)
        hzone$density=0.05*hzone$density
        houtzone$density=0.05*houtzone$density
        plot(hzone,freq=F,xlim=c(0,2),ylim=c(0,1),add=F,col=1,xlab="SCD", ylab="Relative frequency",main=classnames[i])
        plot(houtzone,freq=F,xlim=c(0,2),add=T,density=0,border="grey",col=4)
        par(xpd=T)
        arrows(maxslope,-.1,maxslope,0,length=0.1,lwd=2,col="red")
        text(1.7,.75,paste("#Within = ",list1$within,"\n#Outside = ",list1$outside,sep=""))
        par(xpd=F)

        #hist(list1$zonezone,breaks=seq(0,2,0.05),freq=F,xlim=c(0,2),add=T,col=1,main=classnames[i])
        #hist(list1$outzone,breaks=seq(0,2,0.05),freq=F,xlim=c(0,2),add=T,col=4)

#  4. Plot the ROC curve and include on this the point and tangent line where optimal SCD is found.
        
        
        fpe=list1$ROC$FPE
        tpf=list1$ROC$TPF
        plot(c(0,fpe),c(0,tpf),type="l",xlim=c(0,1),ylim=c(0,1),xlab="1 - TNF (sensitivity)", ylab = "TPF (specificity)",main=classnames[i])
        abline(0,1,lty=2,col=4)
        a=tpf[maxslopeloc]-(fpe[maxslopeloc])
        abline(a,1)
        points(fpe[maxslopeloc],tpf[maxslopeloc],pch=16)
        text(.7,.2,paste("AUC = ",signif(list1$ROC$AUC,3),"\nSE = ",signif(list1$ROC$SEAUC,3),sep=""))
        
#  5. Create barplot of slope over range of SCD.

        barplot(list1$ROC$slope,names.arg= list1$ROC$threasholds,density=0,ylim=c(0,1),xlab = "Dissimilarity", ylab = "TPF - (1 - TNF)", main="TPF-FPE over range of dissimilarity")

        maxpoly=cbind(x=c(maxslopeloc,maxslopeloc+1,maxslopeloc+1,maxslopeloc),y=c(0,0,list1$ROC$slope[maxslopeloc],list1$ROC$slope[maxslopeloc]))
        polygon(maxpoly,density=-1,col="red")

        par(usr=c(0,1,0,1))
        scdopt=as.character(format(signif(maxslope,3),digits=3))
        text(.7,.9,paste("Optimal SCD = ",scdopt,"\n(value of joint minimization)",sep=""))

        #Change plot device after 4th plot
         
                if(i%%4==0){
                        devcnt=devcnt+1
                        dev.set(devn[devcnt])
                }
  }
}
