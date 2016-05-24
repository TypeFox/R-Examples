area=function(x){

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
  appo.profile[1,]=round(appo.profile[1,])

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


  # table with area

  tabarea=data.frame(appo.surface)
  for (i in 1:ncol(x)) colnames(tabarea)="surface area"
  for (i in 1:nrow(x)) rownames(tabarea)=paste("community",c(1:nrow(tabarea)),sep=" n.")
  tabarea


}
