ranking=function(x){


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
  appo.tratti
  tratti=appo.tratti[-nrow(appo.beta),]


  # temp vector for arc

  appo.arc=rep(NA,nsiti)

  # compute arc

  for (i in 1:nsiti) appo.arc[i]=sum(tratti[,i])


  # table with arc

  tabarc=data.frame(appo.arc)
  for (i in 1:ncol(x)) colnames(tabarc)="arc length"
  for (i in 1:nrow(x)) rownames(tabarc)=paste("community",c(1:nrow(tabarc)),sep=" n.")


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


  results=data.frame(appo.sintesi)
  for (i in 1:ncol(x)) colnames(results)=c("Richness","Shannon","Simpson","Arc","Area")
  for (i in 1:nrow(x)) rownames(results)=paste("community",c(1:nrow(results)),sep=" n.")


  # summary table ordred

  ranking=results[with(results, order(-Area)), ]
  print(ranking)




}
