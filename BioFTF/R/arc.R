arc=function(x){

  #################################################################################
  ############################      DATA    #######################################
  #################################################################################
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
  ##############################     ARC LENGTH   #################################
  #################################################################################

  # temp matrix for computing parts of arc length

  appo.tratti=matrix(NA,length(b),nsiti)

  # compute parts for arc

  for (i in 1:all) appo.tratti[i]=sqrt((appo.beta[i+1]-appo.beta[i])^2+(appo.profile[i]-appo.profile[i+1])^2)

  tratti=appo.tratti[-nrow(appo.beta),]


  # temp vector for arc

  appo.arc=rep(NA,nsiti)

  # compute arc

  for (i in 1:nsiti) appo.arc[i]=sum(tratti[,i])


  # table with arc

  tabarc=data.frame(appo.arc)
  for (i in 1:ncol(x)) colnames(tabarc)="arc length"
  for (i in 1:nrow(x)) rownames(tabarc)=paste("community",c(1:nrow(tabarc)),sep=" n.")
  tabarc


}
