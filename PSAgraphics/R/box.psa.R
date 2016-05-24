box.psa<-function(continuous, treatment = NULL, strata = NULL, 
 boxwex = 0.17, offset = 0.17,
 col = c("yellow","orange","black","red","darkorange3"), xlab="Stratum", 
 legend.xy=NULL, legend.labels = NULL,
 pts = TRUE, balance = FALSE, trim = 0, B=1000, ...){
 
#Plots side by side boxplot of C/T for each strata.  Means may be connected if desired.
#Point clouds for each boxplot may be plotted.  The nonparametric Kolmogorov-Smirnov test
#may be used to find p-values for the test of equivalence of distribution between C/T in
#each stratum.  In legend.lables, note that treatment levels A and B are actually taken from the treatment factor.
 
#If "continuous" has three columns, treat as m, t, s.
if(dim(as.data.frame(continuous))[2]==3){ treatment   <- continuous[,2]
                                           strata     <- continuous[,3]
                                           continuous <- continuous[,1]}
cts<-data.frame(continuous,treatment,strata)

#Sort the unique treatment levels.  To be used in legend as well.
sut<-sort(unique(treatment))

#Getting the legend labels sorted out.
leg.test<-is.null(legend.labels)
if(balance){
if(leg.test){
    legend.labels<-c(paste("Treatment",sut[1]),paste("Treatment",sut[2]),"Treatment Means Compared", 
    "KS p-values", "Strata-Treatment Size")}}else{legend.labels<-legend.labels}
if(!balance){
if(leg.test){
    legend.labels<-c(paste("Treatment",sut[1]),paste("Treatment",sut[2]),"Treatment Means Compared", 
    "Strata-Treatment Size")}}else{legend.labels<-legend.labels}

cs.0<-subset(cts,treatment==sut[1],select=c(continuous,strata))
cs.1<-subset(cts,treatment==sut[2],select=c(continuous,strata))  
size<-table(treatment,strata)
       
if(balance){bal.ms<-bal.ms.psa(continuous,treatment,strata,B,trim=trim)
            cat("Press <enter> for bar chart...")
            readline()} 

s.d<-dim(table(strata))
rgy<-range(continuous)
ht<-rgy[2]-rgy[1]
if(is.null(legend.xy)){legend.xy<-c(1,rgy[2]+.22*ht)}

boxplot(continuous ~ strata,
        boxwex = boxwex, at = 1:s.d - offset,axes=FALSE,ylim=c(rgy[1]-.13*ht,rgy[2]+.17*ht),
        subset= treatment == sut[1], col=col[1],
        xlab=xlab, ...)

boxplot(continuous ~ strata,  add = TRUE,axes=FALSE,
        boxwex = boxwex, at = 1:s.d + offset,
        subset= treatment == sut[2], col=col[2], ...)

#Add point clouds to boxplots
if(pts){
  for(i in 1:s.d){cs.0.strat<-subset(cs.0,cs.0[,2]==sort(unique(strata))[i])
                  points(jitter(rep(i,length(cs.0.strat[,1]))-offset,amount=.25*boxwex),
                  jitter(cs.0.strat[,1],amount=.10*boxwex),col="yellow3",cex=.8)}
  for(i in 1:s.d){cs.1.strat<-subset(cs.1,cs.1[,2]==sort(unique(strata))[i])
                  points(jitter(rep(i,length(cs.1.strat[,1])),amount=.25*boxwex)+offset,
                  jitter(cs.1.strat[,1],amount=.10*boxwex),col="orange3",cex=.8)}
       }
       
#Add Kolgomorov-Smirnov p-values     
if(balance){p.ks<-round(bal.ks.psa(continuous,treatment,strata),2)
       for (i in 1:s.d) {text  (i,rgy[1]-(rgy[2]-rgy[1])/10,p.ks[i],col=col[4],cex=.7)}   
           }

axis(1,1:s.d,labels=sort(unique(strata)))
axis(2,at=NULL)

#Size of strata
for (i in 1:s.d) {text (i-offset,rgy[1]-(rgy[2]-rgy[1])/7.4,size[1,i],col=col[5],cex=.7)}
for (i in 1:s.d) {text (i+offset,rgy[1]-(rgy[2]-rgy[1])/7.4,size[2,i],col=col[5],cex=.7)}

strata.factor<-as.factor(strata)
levels(strata.factor)<-1:s.d
strata.vector<-as.vector(strata.factor)

#Create lines matching means
x.t<-as.vector(rbind(1:s.d-offset+boxwex/2,1:s.d+offset-boxwex/2,rep("NA",s.d)))
last<-length(x.t)
x.coord<-x.t[-last]
y.t<-aggregate(cts[,1],by=list(treatment=treatment,strata=strata),FUN=mean,trim=trim)
y.tt<-NULL
for(i in 1:s.d){y.tt<-c(y.tt,y.t[(2*i-1),3],y.t[2*i,3],"NA")}
y.coord<-y.tt[-last]
lines(x.coord, y.coord, col=col[3],lty=1,type='o',lwd=2)

#Legend
if(balance){legend(x=legend.xy[1],y=legend.xy[2], legend = legend.labels,
           col=col, pch=c(15,15,-1,15,15),lty=c(-1,-1,1,-1,-1),
           lwd=c(2,2,2,2),cex=.8)}
           else{legend(x=legend.xy[1],y=legend.xy[2], legend = legend.labels,
                col=col[c(1,2,3,5)], pch=c(15,15,-1,15),lty=c(-1,-1,1,-1),lwd=c(2,2,2,2),cex=.8)
               }         
}



