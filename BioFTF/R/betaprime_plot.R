betaprime_plot=function(x){
  library(BioFTF)
  for (i in 1:ncol(x)) colnames(x)=paste("species",c(1:ncol(x)),sep=" n.")
  for (i in 1:nrow(x)) rownames(x)=paste("community",c(1:nrow(x)),sep=" n.")


  class(x)
  x=as.matrix(x)

  nsiti=nrow(x)
  nspecie=ncol(x)

  # matrix with relative abundance

  xrel=prop.table(x,1)


  #################################################################################
  ###############################     DOMAIN    ###################################
  #################################################################################

  # fix points of domain

  b=seq(-1,1,0.1)
  b[11]=0.01    #avoid "shannon" jump close to 0
  b[1]=-0.9999  #avoid "richness" jump close to -1

  # temp matrix length(b)

  appo.beta=matrix(rep(b),length(b),nsiti)

  #################################################################################
  ##############################     PROFILES    ##################################
  #################################################################################

  # temp matrix for computing beta profile

  appo.profile=matrix(NA,length(b),nsiti)

  appo.siti=matrix(rep((xrel),each=length(b)),nsiti*length(b),nspecie)


  # compute beta profile

  all=length(b)*nsiti

  for(i in 1:all) {appo.profile[i]=(1-sum(appo.siti[i,]^(appo.beta[i]+1)))/(appo.beta[i])}
  appo.profile[1,]=round(appo.profile[1,])  #the richness must be integer

  #################################################################################
  #########################       FIRST           #################################
  #################################################################################

  # temp matrix for computing betaprofile first derivative

  appo.first=matrix(NA,length(b),nsiti)

  # compute beta profile first derivative without species with null frequencies


  appo.siti2=log(appo.siti)
  appo.siti2[appo.siti2==-Inf]=0

  for (i in 1:all) appo.first[i]=(sum(appo.siti[i,]^(appo.beta[i]+1)*(1-appo.beta[i]*appo.siti2[i,]))-1)/((appo.beta[i])^2)

  for (i in 1:nsiti) {plot(b,appo.first[,i],type="l",xlim=c(-1,1),ylim=c(min(appo.first),0),lty=i,lwd=1,col=i,xlab="Beta",ylab="First Derivative",main="First Derivative")
    par(new=TRUE)}
  legend(1,max(appo.first-1), paste("Com.",c(1:nsiti)),lty=c(1:nsiti),y.intersp=0.65,ncol=2,lwd=1,col=c(1:nsiti),xjust=1,merge=TRUE)


}
