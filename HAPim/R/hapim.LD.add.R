`hapim.LD.add` <-
function(hap.trans.pere,hap.trans.mere,perf,CD,map,position,temps.depart,perfectLD,marq.hap.left) {

# vérification des dimensions:
#############################
     if(dim(hap.trans.mere)[2] != dim(hap.trans.pere)[2]  | dim(hap.trans.pere)[2] != (length(map)+1) ) stop("the numbers of genotype information are not consistent in  hap.trans.mere, hap.trans.pere or map ",call.=FALSE)

     if(length(CD) != length(perf) | length(perf) != dim(hap.trans.mere)[1] | dim(hap.trans.mere)[1] != dim(hap.trans.pere)[1]  ) stop("the numbers of individual information are not consistent in CD, perf, hap.trans.pere or hap.trans.mere",call.=FALSE)

#Recodage des haplotypes
########################
	hap.pop       = 	rbind(hap.trans.pere,hap.trans.mere)
	all.marq      = 	allele.marq(hap.pop)
	hap.trans.pere= 	recode.hap(hap.trans.pere,all.marq)
	hap.trans.mere= 	recode.hap(hap.trans.mere,all.marq)
	hap.pop       = 	rbind(hap.trans.pere,hap.trans.mere)

  	
#Calcul du nombre d'intervalles
###############################
	nbre.int	= 	length(map) 


# calcul du nombre de marques de l'haplotype
# même nbre de marques à droite et à gauche de la position de test
##################################################################
	 marq.hap	=	2*marq.hap.left


# distance des marques par rapport à l'origine
##############################################
	 dist.marq 	= 	distance.marqueurs(map)




#########################################################
#  vérification: les positions rentrées par l'utilisateur doivent être 
#  comprises entre marq.hap.left et nbre.marque-marq.hap.left+1
#########################################
         position=sort(position) #rangement 
         position=round(position,5) #Il faut arrondir!

    dist.marq1=dist.marq[marq.hap.left:(length(dist.marq)-marq.hap.left+1)]

         borne.inf=round(dist.marq1[1],5) #Il faut arrondir!
         borne.sup=round(dist.marq1[length(dist.marq1)],5) #Il faut arrondir!


         diff.left=position-borne.inf 
         diff.right=position-borne.sup 

         which=(1:length(position))[(diff.left>=0)&(diff.right<=0)]
         position.new=round(sort(position[which]),5)
         if ((length(position.new))!=length(position)) stop("error in test positions",call.=FALSE) 





#le nombre de positions de tests
################################
        num.pos	   =    length(position.new)


# frequence allelique
#####################
	freq.a	=	freq.all(hap.pop)


# distance des positions de tests par rapport à l'origine
##########################################################
	dist.test   = 	distance.test(position.new,dist.marq)



# calcul de la vraisemblance sous H0 (Pas de QTL)
#################################################
	ML.H0		= 	obj.LD.add.H0(perf,CD)


# Boucle
########

# Initialisation des vecteurs/matrices résultats:LRT,pos.test,assoc.est,hap.ass.est,param.est
	 LRT       	= 	rep(NA,num.pos)
	 pos.test  	= 	rep(NA,num.pos)
	 assoc.est	= 	rep(NA,num.pos)
	 hap.ass.est 	= 	rep(NA,num.pos)
	 if(perfectLD==TRUE){
	      param.est = matrix(NA,ncol=num.pos,nrow=5)
	 }
	 if(perfectLD==FALSE){
	      param.est = matrix(NA,ncol=num.pos,nrow=6)
	 }

# Initialisation de l'index comptant les positions de test:i.pos
         i.pos=0 



   for (i in marq.hap.left:(nbre.int-marq.hap.left+1)) { # Début de la boucle sur les intervalles (FOR1)

   #extraction des informations pour l'intervalle i
   nbre.all.marq	=	rep(0,marq.hap)
   freq.marq		=	list()
   all.marq.int		=	list()

       for(im in 1:(marq.hap))  { # Début FOR2
          nbre.all.marq[im]	=	length(all.marq[[(i-marq.hap.left)+im]])
          freq.marq[[im]]	=	freq.a[[(i-marq.hap.left)+im]]
          all.marq.int[[im]]	=	all.marq[[(i-marq.hap.left)+im]]
       } # Fin FOR2

   #calcul des structures pour l'intervalle i
   res.structure	=	structure.hap(marq.hap,nbre.all.marq)
   pi.hap		=	pi.hap(freq.marq,res.structure)
   cor.pere		=	corresp(hap.trans.pere[,((i-marq.hap.left)+1):((i-marq.hap.left)+marq.hap)],res.structure)
   cor.mere		=	corresp(hap.trans.mere[,((i-marq.hap.left)+1):((i-marq.hap.left)+marq.hap)],res.structure)
   hap.assoc		=	unique(c(cor.pere$assoc,cor.mere$assoc))   
   nbre.ass		= 	length(hap.assoc) 
   nbre.pas		=	length(dist.test[[i]])
   param		=	depart.LD(perf,CD)


          for(h in 1:nbre.pas){ # Début de la boucle sur les pas des intervalles (FOR3)

              for (z in 1:length(position.new)) { # début FOR4

                if (((h==1)&(position.new[z]==dist.test[[i]][h]))|((h>=2)&(position.new[z]==dist.test[[i]][h]))) { # début FOR5

                   i.pos	=	i.pos+1
                   i.dist	=	c(dist.marq[(i-marq.hap.left+1):i],dist.test[[i]][h],dist.marq[(i+1):(i+marq.hap.left)])
                   poids.D	=	poids.D(i.dist,marq.hap,res.structure)

                   #initialisation des vecteurs nécessaires dépendant du nombre d'associations
                   ML	=	rep(NA,nbre.ass)
                   if(perfectLD==TRUE){
                       val.par	=	matrix(NA,nrow=nbre.ass,ncol=5)
                   }
                   if(perfectLD==FALSE){
                       val.par	=	matrix(NA,nrow=nbre.ass,ncol=6)
                   }

    
                    for (j in 1:nbre.ass){# Début de la boucle sur toutes les associations possibles (FOR6)
                         assoc 	= 	hap.assoc[j]
                         don 	= 	list(perf,CD,assoc,cor.pere,cor.mere,pi.hap,res.structure,poids.D)

                         # recherche du minimum de la fonction obj.LDL.add (inverse de la vraisemblance)
                         if(perfectLD==TRUE){
                              start    	=  	c(param[2],0,temps.depart,0.25,param[1])
                              op  	= 	optim(start,obj.LD.add,NULL,method="BFGS",lower=-Inf,upper=Inf,control=list(),hessian=FALSE,don)
                         }

                         if(perfectLD==FALSE){
                              start    	= 	c(param[2],0,temps.depart,0.25,0.5,param[1])
#                              op  	= 	optim(start,obj.LD.add.alpha,NULL,method="BFGS",lower=-Inf,upper=Inf,control=list(),hessian=FALSE,don)
                         }
   
                         val.par[j,]    =        op$par  
                         ML[j] 	        =       -op$value

                   }# Fin de la boucle sur les associations (FOR6)

               # Recherche de la valeur du maximum de vraisemblance  pour chaque position
               # sur les association possibles.
               # Calcul des valeurs du test et stockage des paramètres estimés

               pos.test[i.pos] 	= 	dist.test[[i]][h]
               LRT[i.pos] 	= 	-2*(ML.H0- max(ML))
               i.val 	        =       which.max(ML)
               assoc.est[i.pos] = 	hap.assoc[i.val]
               param.est[,i.pos]=	val.par[i.val,]

               #retrouver les allèles de l'haplotype associé
               hap.ass.est[i.pos]	=	retrouve.all(assoc.est[i.pos],res.structure,all.marq.int)

              } # Fin FOR5
           } # Fin FOR4
        }# Fin de la boucle sur les pas des intervalles (FOR3)
     }# Fin de la boucle sur les intervalles (FOR1)

# On met à 0 les LRT négatifs
  LRT[LRT<0]	=	0

# On regroupe les différents résultats dans un tableau
  res		=	data.frame(pos.test,round(LRT,4),hap.ass.est,round(t(param.est),4))

# On nomme les différentes colonnes du tableau
	# cas1: perfectLD=true

          if(perfectLD==TRUE){
            dimnames(res)[[2]]	=	c("position","LRT","haplotype","variance","effect.Q","temps","pi.Q(0)","mean")
          }

	# cas2: perfectLD=false

          if(perfectLD==FALSE){
            dimnames(res)[[2]]	=	c("position","LRT","haplotype","variance","effect.Q","temps","alpha","pi.Q(0)","mean")
          } 
  

res

}

