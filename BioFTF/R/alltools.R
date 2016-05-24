alltools=function(x){

  #################################################################################
  ############################      DATA    #######################################
  #################################################################################
  library(BioFTF)
  for (i in 1:ncol(x)) colnames(x)=paste("species",c(1:ncol(x)),sep=" n.")
  for (i in 1:nrow(x)) rownames(x)=paste("community",c(1:nrow(x)),sep=" n.")
  x

  class(x)
  x=as.matrix(x)

  nsiti=nrow(x)
  nspecie=ncol(x)

  # matrix with relative abundance

  xrel=prop.table(x,1)
  xrel

  #################################################################################
  ###############################     DOMAIN    ###################################
  #################################################################################

  # fix points of domain

  b=seq(-1,1,0.1)
  b[11]=0.01    #avoid "shannon" jump close to 0
  b[1]=-0.9999  #avoid "richness" jump close to -1

  # temp matrix length(b)

  appo.beta=matrix(rep(b),length(b),nsiti)
  appo.beta

  #################################################################################
  ##############################     PROFILES    ##################################
  #################################################################################

  # temp matrix for computing beta profile

  appo.profile=matrix(NA,length(b),nsiti)
  appo.profile
  appo.siti=matrix(rep((xrel),each=length(b)),nsiti*length(b),nspecie)
  appo.siti

  # compute beta profile

  all=length(b)*nsiti

  for(i in 1:all) {appo.profile[i]=(1-sum(appo.siti[i,]^(appo.beta[i]+1)))/(appo.beta[i])}
  appo.profile[1,]=round(appo.profile[1,])  #the richness must be integer

  # plot beta profile

  par(mfrow=c(2,3))
  for (i in 1:nsiti) {plot(b,appo.profile[,i],type="l",xlim=c(-1,1),ylim=c(0,(nspecie)),lty=i,col=i,xlab="Beta",ylab="Diversity",main="Diversity Profiles")
    par(new=TRUE)}
legend(1,nspecie, paste("Com.",c(1:nsiti)),lty=c(1:nsiti),lwd=1,y.intersp=0.65,ncol=1,col=c(1:nsiti),xjust=1,merge=TRUE)


  #################################################################################
  #########################       FIRST           #################################
  #################################################################################

  # temp matrix for computing betaprofile first derivative

  appo.first=matrix(NA,length(b),nsiti)

  # compute beta profile first derivative without species with null frequencies

  appo.siti
  appo.siti2=log(appo.siti)
  appo.siti2[appo.siti2==-Inf]=0
  appo.siti2

  for (i in 1:all) appo.first[i]=(sum(appo.siti[i,]^(appo.beta[i]+1)*(1-appo.beta[i]*appo.siti2[i,]))-1)/((appo.beta[i])^2)
  appo.first

  # plot beta profile first derivative

  for (i in 1:nsiti) {plot(b,appo.first[,i],type="l",xlim=c(-1,1),ylim=c(min(appo.first),0),lty=i,lwd=1,col=i,xlab="Beta",ylab="First Derivative",main="First Derivative")
    par(new=TRUE)}
 legend(1,max(appo.first-1), paste("Com.",c(1:nsiti)),y.intersp=0.65,ncol=1,lty=c(1:nsiti),lwd=1,col=c(1:nsiti),xjust=1,merge=TRUE)

  #################################################################################
  #############################     SECOND     ####################################
  #################################################################################

  # temp matrix for computing betaprofile second derivative

  appo.second=matrix(NA,length(b),nsiti)

  # compute beta profile second derivative without species with null frequencies

  for (i in 1:all) appo.second[i]=(sum((appo.siti[i,]^(appo.beta[i]+1))*(-(appo.beta[i])^2*((appo.siti2[i,])^2)-2+2*appo.beta[i]*appo.siti2[i,]))+2)/((appo.beta[i])^3)
  appo.second

  # plot beta profile second derivative

  for (i in 1:nsiti)
  {plot(b,appo.second[,i],type="l",xlim=c(-1,1),ylim=c(0,max(appo.second)),lty=i,lwd=1,col=i,xlab="Beta",ylab="Second Derivative",main="Second Derivative")
    par(new=TRUE)}
  legend(1,max(appo.second), paste("Com.",c(1:nsiti)),y.intersp=0.65,ncol=1,lty=c(1:nsiti),lwd=1,col=c(1:nsiti),xjust=1,merge=TRUE)

  #################################################################################
  ##############################    CURVATURE     #################################
  #################################################################################

  # temp matrix for computing betaprofile curvature

  appo.curva=matrix(NA,length(b),nsiti)

  # compute beta profile curvature

  for (i in 1:all) appo.curva[i]=abs(appo.second[i])/(1+appo.first[i]^2)^(3/2)
  appo.curva

  # plot beta profile curvature

  for (i in 1:nsiti)
  {plot(b,appo.curva[,i],type="l",xlim=c(-1,1),ylim=c(min(appo.curva),max(appo.curva)),lty=i,lwd=1,col=i,xlab="Beta",ylab="Curvature",main="Curvature")
    par(new=TRUE)}
  legend(1,max(appo.curva), paste("Com.",c(1:nsiti)),lty=c(1:nsiti),y.intersp=0.65,ncol=1,lwd=1,col=c(1:nsiti),xjust=1,merge=TRUE)

  #################################################################################
  ##############################      RADIUS      #################################
  #################################################################################

  # temp matrix for computing beta profile radius

  appo.radius=matrix(NA,length(b),nsiti)

  # compute beta profile radius

  for (i in 1:all) appo.radius[i]=((1+appo.first[i]^2)^(3/2))/abs(appo.second[i])
  appo.radius

  # plot beta profile radius

  for (i in 1:nsiti)
  {plot(b,appo.radius[,i],type="l",xlim=c(-1,1),ylim=c(0,max(appo.radius)),lty=i,lwd=1,col=i,xlab="Beta",ylab="Radius",main="Radius of curvature")
    par(new=TRUE)}
   legend(1,max(appo.radius), paste("Com.",c(1:nsiti)),lty=c(1:nsiti),y.intersp=0.65,ncol=1,lwd=1,col=c(1:nsiti),xjust=1,merge=TRUE)

  # plot beta profile radius zoom

  for (i in 1:nsiti)
  {plot(b,appo.radius[,i],type="l",xlim=c(-0.5,0.5),ylim=c(0,mean(appo.radius)),lty=i,lwd=1,col=i,xlab="Beta",ylab="Radius",main="Radius of curvature")
    par(new=TRUE)}
  legend(0.5,mean(appo.radius), paste("Com.",c(1:nsiti)),lty=c(1:nsiti),y.intersp=0.65,ncol=1,lwd=1,col=c(1:nsiti),xjust=1,merge=TRUE)

  #################################################################################
  ##############################     ARC LENGTH   #################################
  #################################################################################

  # temp matrix for computing parts of arc length

  appo.tratti=matrix(NA,length(b),nsiti)

  # compute parts for arc

  for (i in 1:all) appo.tratti[i]=sqrt((appo.beta[i+1]-appo.beta[i])^2+(appo.profile[i]-appo.profile[i+1])^2)
  appo.tratti
  tratti=appo.tratti[-nrow(appo.beta),]
  tratti

  # temp vector for arc

  appo.arc=rep(NA,nsiti)

  # compute arc

  for (i in 1:nsiti) appo.arc[i]=sum(tratti[,i])
  appo.arc

  # table with arc

  tabarc=data.frame(appo.arc)
  for (i in 1:ncol(x)) colnames(tabarc)="arc length"
  for (i in 1:nrow(x)) rownames(tabarc)=paste("community",c(1:nrow(tabarc)),sep=" n.")
  tabarc

  #################################################################################
  ##############################       AREA       #################################
  #################################################################################

  # temp matrix and vector for computing surface area

  appo.aree=matrix(NA,length(b),nsiti)
  appo.surface=rep(NA,nsiti)

  # compute area

  for (i in 1:all) appo.aree[i]=(appo.beta[i+1]-appo.beta[i])*(appo.profile[i+1]+appo.profile[i])/2
  aree=appo.aree[-nrow(appo.beta),]
  for (i in 1:nsiti) appo.surface[i]=sum(aree[,i])
  appo.surface

  # table with area

  tabarea=data.frame(appo.surface)
  for (i in 1:ncol(x)) colnames(tabarea)="surface area"
  for (i in 1:nrow(x)) rownames(tabarea)=paste("community",c(1:nrow(tabarea)),sep=" n.")
  tabarea

  #################################################################################
  ##############################      TABLES     #################################
  #################################################################################

  # temp matrix for summary

  appo.sintesi=matrix(NA,nsiti,5)

  # summary table

  appo.sintesi[,1]=appo.profile[1,]
  appo.sintesi[,2]=appo.profile[11,]
  appo.sintesi[,3]=appo.profile[21,]
  appo.sintesi[,4]=appo.arc
  appo.sintesi[,5]=appo.surface
  appo.sintesi

  results=data.frame(appo.sintesi)
  for (i in 1:ncol(x)) colnames(results)=c("Richness","Shannon","Simpson","Arc","Area")
  for (i in 1:nrow(x)) rownames(results)=paste("community",c(1:nrow(results)),sep=" n.")
  results

  # summary table ordered

  ranking=results[with(results, order(-Area,-Richness,-Shannon,-Arc)), ]
  print("Communities Ranking")
  ranking

}
