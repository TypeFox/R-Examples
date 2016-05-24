beta_plot=function(x){
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


  for (i in 1:nsiti) {plot(b,appo.profile[,i],type="l",xlim=c(-1,1),ylim=c(0,(nspecie)),lty=i,col=i,xlab="Beta",ylab="Diversity",main="Diversity Profiles")
    par(new=TRUE)}
  legend(1,nspecie, paste("Com.",c(1:nsiti)),lty=c(1:nsiti),lwd=1,y.intersp=0.65,ncol=2,col=c(1:nsiti),xjust=1,merge=TRUE)

}
