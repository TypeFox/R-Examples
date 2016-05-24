# plot tracer file from zelig

read.plot.zelig.out <- function(fileout,spp,grp1,grp2,pdfout=F){

if(pdfout==T) pdf(paste(fileout,".pdf",sep=""))
 mat<- matrix(1:2,2,1,byrow=T)
 layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")

nvar <- length(spp)+ 6 + 1

tra <- matrix(scan(fileout), byrow=T, ncol=nvar)
# add timezero to values
time <- c(0,tra[,1])
tot <- rbind(rep(0,4),tra[,c(2,5,6,7)])
tot[,1] <- tot[,1]/100 # per Km2; 

basp<- rbind(rep(0,8),tra[,8:nvar])

timeunit <- "Years"
variab <- c("Density (Stem/ha)/100", "BA (m2/ha)",
            "LAI", "Height (m)")

# plot bas by spp
matplot(time, basp[,grp1], type="l", col=1, xlab=timeunit, ylab="BasalArea (m2/ha)",ylim=c(0,25),
         xlim=c(0,1.5*max(time)),lty=c(1:length(grp1)))
legend(1.05*max(time),20, spp[grp1], lty=c(1:length(grp1)), col=1,cex=0.7)
matplot(time, basp[,grp2], type="l", xlab=timeunit, ylab="BasalArea (m2/ha)",ylim=c(0,25),
        xlim=c(0,1.5*max(time)),lty=c(1:length(grp1)), col=1)
legend(1.05*max(time),20, spp[grp2], lty=c(1:length(grp1)), col=1,cex=0.7)

# plot set of variables aggregated for stand
matplot(time, tot, type="l",col=1,ylab="Stand totals",xlab=timeunit,xlim=c(0,1.5*max(time)))
legend("topright",legend=variab,col=1,lty=1:5,cex=0.7)

if(pdfout==T) dev.off()

}
