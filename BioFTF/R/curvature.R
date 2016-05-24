curvature=function(x){

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
  #########################       FIRST           #################################
  #################################################################################

  # temp matrix for computing betaprofile first derivative

  appo.first=matrix(NA,length(b),nsiti)

  # compute beta profile first derivative without species with null frequencies

  appo.siti2=log(appo.siti)
  appo.siti2[appo.siti2==-Inf]=0

  for (i in 1:all) appo.first[i]=(sum(appo.siti[i,]^(appo.beta[i]+1)*(1-appo.beta[i]*appo.siti2[i,]))-1)/((appo.beta[i])^2)

  #################################################################################
  #############################     SECOND     ####################################
  #################################################################################

  # temp matrix for computing betaprofile second derivative

  appo.second=matrix(NA,length(b),nsiti)

  # compute beta profile second derivative without species with null frequencies

  for (i in 1:all) appo.second[i]=(sum((appo.siti[i,]^(appo.beta[i]+1))*(-(appo.beta[i])^2*((appo.siti2[i,])^2)-2+2*appo.beta[i]*appo.siti2[i,]))+2)/((appo.beta[i])^3)

  #################################################################################
  ##############################    CURVATURE     #################################
  #################################################################################

  # temp matrix for computing betaprofile curvature

  appo.curva=matrix(NA,length(b),nsiti)

  # compute beta profile curvature

  for (i in 1:all) appo.curva[i]=abs(appo.second[i])/(1+appo.first[i]^2)^(3/2)
  appo.curva


}
