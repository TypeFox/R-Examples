################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-03-25
# Chapter 10. Putting it all together in a case study
################################################################################
#' Calcul la vraisemblance d'un jeu de parametre
#' categorie : MCMC
#' lance la fonction moderato_AllCampagnesMCMC 
#' @title vraisemblance
#' @name vraisemblance
#' @param v.param 1 vecteur de parametres 
#' @param All.mesure data.frame de toutes les mesures experimentales
#' @param LAI.sigma variance du LAI fixee
#' @param B.sigma variance de la biomasse fixee
#' @param rdt0.sigma variance du rendement fixee
#' @param sites vecteur de tous les sites a etudier
#' @param All.mesure data.frame de toutes les mesures experimentales
#' @param id_campagne matrice 2 colonnes, nombredecampagne lignes, chaque nom de campagne
#' @param moderato.chemin chemin du dossier moderato
#' @param Gamma booleen : calcul de shape et scale si Gamma=TRUE  (pour Gibbs)
#' @return la vraisemblance du jeu de parametre par rapport aux donnees experimentales des sites choisis.
 
vraisemblance<-function(param, list_sy, data, LAI.sigma,B.sigma,sdate,ldate,Gamma=FALSE){
    #simulations et observations sur les site-year
    sim=maize.multisy(param, list_sy,sdate,ldate)[,c("sy","day","LAI","B")]
    data=subset(data,sy%in%list_sy)[,c("sy","day","LAIobs","Bobs")]
    simobs = merge(sim,data, by=c("sy","day"))

    LAI.err2<-(simobs$LAIobs-simobs$LAI)^2
    B.err2<-(simobs$Bobs-simobs$B)^2
    LAI.SSE<-sum(LAI.err2,na.rm=TRUE)
    B.SSE<-sum(B.err2,na.rm=TRUE)
    LAI.nobs<-sum(!is.na(LAI.err2))
    B.nobs<-sum(!is.na(B.err2))

    LAI.LogVrais<- - (LAI.nobs/2)*log(2*pi*LAI.sigma^2)- (1/(2*LAI.sigma^2))*(LAI.SSE)
    B.LogVrais<- - (B.nobs/2)*log(2*pi*B.sigma^2)- (1/(2*B.sigma^2))*(B.SSE)
    LogVrais<-LAI.LogVrais+B.LogVrais
    
    LAI.MSE<-LAI.SSE/LAI.nobs
    B.MSE<-B.SSE/B.nobs
    resultat<-list("LAI.LogVrais"=LAI.LogVrais, "B.LogVrais"=B.LogVrais, "LAI.MSE"=LAI.MSE, "B.MSE"=B.MSE)

    if(Gamma){
    #calcul des parametres shape et scale de la loi Gamma 
    LAI.shape<-(LAI.nobs/2)+2
    LAI.scale<-2/LAI.SSE
    B.shape<-(B.nobs/2)+2
    B.scale<-2/B.SSE

    LAI.variance.new<-(rgamma(1,shape=LAI.shape,scale=LAI.scale))^-1
    B.variance.new<-(rgamma(1,shape=B.shape,scale=B.scale))^-1   
    LAI.sigma.new<-sqrt(LAI.variance.new)
    B.sigma.new<-sqrt(B.variance.new)
    LAI.LogVrais.new<- - (LAI.nobs/2)*log(2*pi*LAI.sigma.new^2)- (1/(2*LAI.sigma.new^2))*(LAI.SSE)
    B.LogVrais.new<- - (B.nobs/2)*log(2*pi*B.sigma.new^2)- (1/(2*B.sigma.new^2))*(B.SSE)
    LogVrais.new<-LAI.LogVrais.new+B.LogVrais.new
    resultat<-list("LogVrais"=LogVrais,"LAI.MSE"=LAI.MSE, "B.MSE"=B.MSE, "LAI.SSE"=LAI.SSE, "B.SSE"=B.SSE, "LAI.sigma.new"=LAI.sigma.new,"B.sigma.new"=B.sigma.new,"LogVrais.new"=LogVrais.new,"LAI.shape"=LAI.shape,"B.shape"=B.shape,"LAI.scale"=LAI.scale,"B.scale"=B.scale,"LAI.nobs"=LAI.nobs,"B.nobs"=B.nobs)  }
return(resultat)
}


################################################################################
#' pour le cas ou teta1 est refuse : retire logvraisemblance ds la loi de teta0
Calcul.Gamma<-function(LAI.shape,B.shape,LAI.scale,B.scale,LAI.nobs,B.nobs,LAI.SSE,B.SSE)
{
  LAI.variance.new<-(rgamma(1,shape=LAI.shape,scale=LAI.scale))^-1 
  B.variance.new<-(rgamma(1,shape=B.shape,scale=B.scale))^-1   

  LAI.sigma.new<-sqrt(LAI.variance.new) 
  B.sigma.new<-sqrt(B.variance.new)

  LAI.LogVrais.new<- - (LAI.nobs/2)*log(2*pi*LAI.sigma.new^2)- (1/(2*LAI.sigma.new^2))*(LAI.SSE)
  B.LogVrais.new<- - (B.nobs/2)*log(2*pi*B.sigma.new^2)- (1/(2*B.sigma.new^2))*(B.SSE)

  LogVrais.new<-LAI.LogVrais.new+B.LogVrais.new
  resultat<-list("LogVrais.new"=LogVrais.new) 
  return(resultat) 
}


#' similaire a Algorithme de Metropolis-Hastings
#' mais la variance residuelle est calculee a chaque pas au lieu d'etre fixee
#' @name  MetropolisHastings_Gibbs
#' @title Algorithme de Metropolis-Hastings within Gibbs
#' @param  moderato.chemin
#' @param  All.mesure  toutes les mesures
#' @param  coeff.teta0  entier entre 0 et 1, pour definir un point de depart entre borne min et borne max
#' @param  Nb.iterations nombre d'iterations de la chaine de Markov
#' @param  LAI.sigma valeur initiale de LAI.sigma (est recalculee a chaque pas)
#' @param  B.sigma valeur initiale de B.sigma  (est recalculee a chaque pas)
#' @param  NomFichierSortie  nom du fichier ou sera enregistre la chaine
#' @param  LectureParam data.frame contenant les noms des parametres, bornes, parametres a estimer
#' @param  ReprisePointDepart booleen. Reprise d'une ancienne chaine ?
#' @param  NomfichierReprise   si ReprisePointDepart==TRUE, quelle est le nom de l'ancien NomFichierSortie
#' @param  cov  en cas de reprise, utilise t'on une autre matrice de var-covar ?
#MetropolisHastings_Gibbs<-function(vraisemblance, model,parametres.apriori,All.mesure,coeff.teta0,Nb.iterations,LAI.sigma,B.sigma,NomFichierSortie,LectureParam,ReprisePointDepart=FALSE,NomfichierReprise="NA",cov=FALSE)


##test
MetropolisHastings_Gibbs<-function(param.apriori=maize.define.param(),param.opti=c("RUE"),Nb.iterations,LAI.sigma,B.sigma,coeff.teta0=0.25,list_sy,data,sdate,ldate, NomFichierSortie="out")
 { 
 t0<-Sys.time() 
 NomFichierSortieExtension<-paste(NomFichierSortie,".csv",sep="")

 ####################chargement des parametres #################################
 #1/ valeurs des parametres init
 names.param.apriori<-names(param.apriori)

 #2/ les parametres qui vont etre estimes
 param.apriori.opti=as.matrix(param.apriori[,param.opti])
 colnames(param.apriori.opti)= param.opti

# cat("coeff multiplicateur des parametres a priori pour le calcul de teta0 :",coeff.teta0,"\n")

 #3/ Les bornes des parametres a estimer
 valeurMin = param.apriori.opti["binf",]
 valeurMax = param.apriori.opti["bsup",]

 #4/ Matrice des variances : var loi uniforme (b-a)^2/12, b=param*1.3, a=param*0.7
 VarcovPrior<-matrix(0,dim(param.apriori.opti)[2],dim(param.apriori.opti)[2])
 diag(VarcovPrior)<- as.numeric((valeurMax-valeurMin)^2/12)
 p.tune<-1
 ####################initialisation de la chaine################################
 rmq<-0
 teta0<-valeurMin+(valeurMax-valeurMin)*coeff.teta0   #coeff.teta0 entre 0 et 1

 if(any(teta0<valeurMin)|any(teta0>valeurMax)) {cat("Attention, teta0 n'est pas dans l'a priori")}
 param=param.apriori["nominal",]

 param[param.opti]=teta0

 resu_vrais<-vraisemblance(param=param,list_sy=list_sy ,data=data,LAI.sigma,B.sigma,sdate,ldate,Gamma=TRUE)
 LogVrais.teta0<- resu_vrais$LogVrais
 probapost0<-LogVrais.teta0+0

 MSE0<-c(LAI.MSE0=resu_vrais$LAI.MSE,B.MSE0=resu_vrais$B.MSE)
 SSE0<-c(LAI.SSE0=resu_vrais$LAI.SSE,B.SSE0=resu_vrais$B.SSE)
 shape0<-c(LAI.shape=resu_vrais$LAI.shape,B.shape=resu_vrais$B.shape)
 scale0<-c(LAI.scale=resu_vrais$LAI.scale,B.scale=resu_vrais$B.scale)

  ## teta0 est accepte ds la chaine
 Sigma0<-c(LAI.sigma=resu_vrais$LAI.sigma.new,B.sigma=resu_vrais$B.sigma.new) 
 LogVrais.new<- resu_vrais$LogVrais.new 
 probapost.new<- LogVrais.new+0
 
 alpha<-NA
 p.unif<-NA
 recap<-c(t(as.numeric(teta0)),LogVrais.teta0,as.numeric(MSE0),as.numeric(SSE0),t(as.numeric(teta0)),LogVrais.teta0,as.numeric(MSE0),as.numeric(SSE0),as.numeric(Sigma0),alpha,p.unif,probapost0,probapost.new,shape0,scale0,p.tune,rmq)
 names(recap)<-c(paste(param.opti,"P",sep=""),"LogVrais.tetaP","LAI.MSEP","B.MSEP","LAI.SSEP","B.SSEP",param.opti,"LogVrais.teta","LAI.MSE","B.MSE","LAI.SSE","B.SSE","LAI.sigma.new","B.sigma.new","alpha","p.unif","probapost","probapost.new","LAI.shape","B.shape","LAI.scale","B.scale","p.tune","rmq")
 name_recap<-names(recap)
 # Nombre d'observations nobs=SSE/MSE, valeur fixe 
 LAI.nobs<-as.numeric(SSE0["LAI.SSE0"]/MSE0["LAI.MSE0"])
 B.nobs<-as.numeric(SSE0["B.SSE0"]/MSE0["B.MSE0"])
 ### Ecriture des informations
 write(name_recap,NomFichierSortieExtension,ncolumns=length(recap),sep=";")
 write(recap,NomFichierSortieExtension,append=T,ncolumns=length(recap),sep=";")

################################################################################
 for(i in 1:Nb.iterations)
  {

  #recuperation de teta0, de la SEQ de teta0, de la variance des erreurs et de la log vraisemblance de teta0 : 2eme partie du vect recap

  teta0<-data.frame(t(recap[param.opti]))  
  LogVrais.teta0<-recap["LogVrais.teta"]
  probapost0 <-recap["probapost"]
  probapost.new <-recap["probapost.new"]
  MSE0<-c(recap["LAI.MSE"],recap["B.MSE"])
  SSE0<-c(recap["LAI.SSE"],recap["B.SSE"])
  shape0<-c(recap["LAI.shape"],recap["B.shape"])
  scale0<-c(recap["LAI.scale"],recap["B.scale"])

  Sigma0<-c(recap["LAI.sigma.new"],recap["B.sigma.new"])

  teta1<-teta0
  teta1<-as.vector(rmnorm(1,as.numeric(teta0),p.tune*VarcovPrior))
  param[param.opti]=teta1
  
  in.apriori<-(all(teta1>=valeurMin)& all(teta1<=valeurMax))

  if(!in.apriori) 
     {  # Cas 1
     rmq<-1
     LogVrais.new<-Calcul.Gamma(LAI.shape=as.numeric(shape0["LAI.shape"]),B.shape=as.numeric(shape0["B.shape"]),LAI.scale=as.numeric(scale0["LAI.scale"]),B.scale=as.numeric(scale0["B.scale"]),LAI.nobs=LAI.nobs,B.nobs=B.nobs,LAI.SSE=as.numeric(SSE0["LAI.SSE"]),B.SSE=as.numeric(SSE0["B.SSE"]))

     recap<-c(t(as.numeric(teta1)),as.numeric(LogVrais.teta0),as.numeric(MSE0),as.numeric(SSE0),t(as.numeric(teta0)),as.numeric(LogVrais.new),as.numeric(MSE0),as.numeric(SSE0),as.numeric(Sigma0),NA,NA,probapost0,probapost.new,shape0,scale0,p.tune,rmq)
     names(recap)<-name_recap
     }
  else
     {  #Cas 2      



     resu_vrais1<-vraisemblance(param=param,list_sy=list_sy ,data=data,LAI.sigma=as.numeric(Sigma0["LAI.sigma.new"]),B.sigma=as.numeric(Sigma0["B.sigma.new"]),sdate,ldate,Gamma=TRUE)

     MSE1<-c(resu_vrais1$LAI.MSE,resu_vrais1$B.MSE)
     SSE1<--c(resu_vrais1$LAI.SSE,resu_vrais1$B.SSE)
     Sigma1<-c(resu_vrais1$LAI.sigma.new,resu_vrais1$B.sigma.new)
     LogVrais.teta1<- resu_vrais1$LogVrais
     shape1<-c(resu_vrais1$LAI.shape,resu_vrais1$B.shape) 
     scale1<-c(resu_vrais1$LAI.scale,resu_vrais1$B.scale) 

     Logprior1<-0
     probapost1<-LogVrais.teta1 + Logprior1   #calcul de la proba a posteriori de teta1P
     #Calcul de alpha : rapport des proba a posteriori
     alpha<-as.numeric(exp(probapost1-probapost.new))

     if (alpha>1)           
          { #CAS 2.1 : alpha est superieur a 1 -> teta1P = teta1 (la proposition est retenue dans la chaine)
          rmq<-2.1
          Accept=TRUE
          }
     else 
          { #CAS 2.2 : alpha est inferieur a 1 -> On tire un nombre entre 0 et 1 selon une loi uniforme
          p.unif<-runif(1,0,1)      #p.unif = proba d'accepter le parametre

          if (p.unif>=alpha)
               { #CAS 2.2.1 : Si la p.unif > = alpha , teta1* = teta0
               rmq<-2.21
               Accept=FALSE
               }
          else 
               { #CAS 2.2.2 : Si la p.unif < alpha , teta1* = teta1
               rmq<-2.22
               Accept=TRUE
               }              
          }
      #on met les resultats dans un vecteur recapitulatif   
     if(Accept) 
          {
          LogVrais.new<- resu_vrais1$LogVrais.new
          probapost.new<- LogVrais.new+0
          p.unif<-NA
          recap<-c(t(as.numeric(teta1)),as.numeric(LogVrais.teta1),as.numeric(MSE1),as.numeric(SSE1),t(as.numeric(teta1)),LogVrais.teta1,as.numeric(MSE1),as.numeric(SSE1),as.numeric(Sigma1),alpha,p.unif,probapost1,probapost.new,shape1,scale1,p.tune,rmq)
          names(recap)<-name_recap

          }
     else
          { 
           LogVrais.new<-Calcul.Gamma(LAI.shape=as.numeric(shape0["LAI.shape"]),B.shape=as.numeric(shape0["B.shape"]),LAI.scale=as.numeric(scale0["LAI.scale"]),B.scale=as.numeric(scale0["B.scale"]),LAI.nobs=LAI.nobs,B.nobs=B.nobs,LAI.SSE=as.numeric(SSE0["LAI.SSE"]),B.SSE=as.numeric(SSE0["B.SSE"]))
          recap<-c(t(as.numeric(teta1)),as.numeric(LogVrais.teta1),as.numeric(MSE1),as.numeric(SSE1),t(as.numeric(teta0)),as.numeric(LogVrais.new),as.numeric(MSE0),as.numeric(SSE0),as.numeric(Sigma0),alpha,p.unif,probapost0,probapost.new,shape0,scale0,p.tune,rmq)
          names(recap)<-name_recap
          }
     } # FIN CAS 2
  write(recap,NomFichierSortieExtension,ncolumns=length(recap), append=T,sep=";")   
  # tunning 
  param_tun<-100
  if (i%%param_tun==0 & i>=50)
   {
   MCMC.result<-read.csv(NomFichierSortieExtension,sep=";")
   MCMC.result<-MCMC.result[(i-param_tun+1):i,]
   nb_rejet<-length(MCMC.result[MCMC.result$rmq==2.21,1])+length(MCMC.result[MCMC.result$rmq==1,1])
   taux_rejet<-nb_rejet/param_tun
   if (taux_rejet<0.65) {p.tune<-p.tune*2}
   if (taux_rejet>0.85) {p.tune<-p.tune/2}
    }
  #sauvegarde des resultats  
  sauv<-10000    
  if (i%%sauv==0)
     {
     print(i) 
     print(Sys.time()-t0)
     sauv.fichier<-read.csv(NomFichierSortieExtension,sep=";")
     write.table(sauv.fichier,paste(NomFichierSortie,"_",i,".csv",sep=""),sep=";",row.names=FALSE)
     if(i!=10000){file.remove(paste(NomFichierSortie,"_",i-10000,".csv",sep=""))}
     }
  } # fin for(i in 1:Nb.iterations)   
 } #fin de la fonction




#' Lit les chaines et les fixe a la meme longueur
#' @name Lecture.chaines
#' @title Lecture des chaines
#' @param Nomchaines  vecteur contenant le nom des chaines
#' @param CheminChaines chemin des chaines
#' @return les chaines
Lecture.chaines<-function(NomChaines,CheminChaines=NA)
{
nb_chaines<-length(NomChaines)

for(i in 1:nb_chaines){
    nom=paste("MCMC.result",i,sep="")
    if(is.na(CheminChaines)) nom_fichier=NomChaines[i]
    else nom_fichier=paste(CheminChaines,NomChaines[i],sep="//")
    assign(nom,read.csv(nom_fichier,sep=";"))}

  if(nb_chaines==2){
if(dim(MCMC.result1)[1]!=dim(MCMC.result2)[1])
{
  dim.min<-min(dim(MCMC.result1)[1],dim(MCMC.result2)[1])
  print(dim.min)
  MCMC.result1<-MCMC.result1[1:dim.min,]
  MCMC.result2<-MCMC.result2[1:dim.min,]
}
return(list(MCMC.result1=MCMC.result1,MCMC.result2=MCMC.result2))}
  if(nb_chaines==4){
if(dim(MCMC.result1)[1]!=dim(MCMC.result2)[1] | dim(MCMC.result2)[1]!=dim(MCMC.result3)[1] | dim(MCMC.result3)[1]!=dim(MCMC.result4)[1] )
{
  dim.min<-min(dim(MCMC.result1)[1],dim(MCMC.result2)[1],dim(MCMC.result3)[1],dim(MCMC.result4)[1])
  print(dim.min)
  MCMC.result1<-MCMC.result1[1:dim.min,]
  MCMC.result2<-MCMC.result2[1:dim.min,]
  MCMC.result3<-MCMC.result3[1:dim.min,]
  MCMC.result4<-MCMC.result4[1:dim.min,]
}
return(list(MCMC.result1=MCMC.result1,MCMC.result2=MCMC.result2,MCMC.result3=MCMC.result3,MCMC.result4=MCMC.result4))}
}