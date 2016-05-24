`hapim.LDL.add` <-
function(hap.chrom1.pere,hap.chrom2.pere,hap.trans.mere,perf,CD,genea,PLA,map,position,temps.depart,perfectLD,marq.hap.left){


# vérification des dimensions:
############################
     if(dim(hap.trans.mere)[2] != dim(hap.chrom1.pere)[2] |  dim(hap.chrom1.pere)[2] != dim(hap.chrom2.pere)[2] | dim(hap.chrom2.pere)[2] != (length(map)+1) ) stop("the numbers of genotype information are not consistent in  hap.trans.mere, hap.chrom1.pere, hap.chrom2.pere or map ",call.=FALSE)

     if(length(CD) != length(perf) | length(perf) != dim(genea)[1] | dim(genea)[1] != dim(PLA)[1] | dim(PLA)[1] != dim(hap.trans.mere)[1]) stop("the numbers of half-sib information are not consistent in CD, perf, genea, PLA or hap.trans.mere",call.=FALSE)

     if(dim(hap.chrom1.pere)[1] != dim(hap.chrom2.pere)[1] | dim(hap.chrom2.pere)[1] != length(unique(genea[,2])))  stop("the numbers of sire information are not consistent in hap.chrom1.pere, hap.chrom2.pere or unique(genea[,2]) ",call.=FALSE)

     if  (dim(PLA)[2]!=length(position)) stop("the number of columns in PLA is not consistent with the number of test positions",call.=FALSE)


# rangement de perf, CD, PLA, hap.trans.mere, hap.chrom1.pere, hap.chrom2.pere en fonction de l'index des pères (2nde colonne de genea)
##################################################################################################################################
#on range perf, CD, PLA et hap.trans.mere en fonction de l'index des pères
	ord.fils=order(genea[,2])
        perf=perf[ord.fils]
        CD=CD[ord.fils]
        hap.trans.mere=hap.trans.mere[ord.fils,]
        PLA=PLA[ord.fils,]
#on range hap.chrom1.pere et hap.chrom2.pere en fonction de l'index des pères
	ord.pere=order(unique(genea[,2]))
        hap.chrom1.pere=hap.chrom1.pere[ord.pere,]
	hap.chrom2.pere=hap.chrom2.pere[ord.pere,]
#on range genea
        genea=genea[ord.fils,]

# recodage des haplotypes
#########################
	 hap.pop	=	rbind(hap.chrom1.pere,hap.chrom2.pere,hap.trans.mere)
	 all.marq 	= 	allele.marq(hap.pop)
	 hap.chrom1.pere= 	recode.hap(hap.chrom1.pere,all.marq)
	 hap.chrom2.pere= 	recode.hap(hap.chrom2.pere,all.marq)
	 hap.trans.mere = 	recode.hap(hap.trans.mere,all.marq)
	 hap.pop	=	rbind(hap.chrom1.pere,hap.chrom2.pere,hap.trans.mere)
 
# calcul du nombre d'intervalles
################################
	 nbre.int	=	length(map) 


# calcul du nombre de marques de l'haplotype: même nbre de marques à droite et à gauche de la position de test
##############################################################################################################
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




# calcul du nombre de positions de tests
########################################
	num.pos		=	length(position.new)

# distance des position de tests par rapport à l'origine
########################################################
        dist.test 	= 	distance.test(position.new,dist.marq)


# frequence allelique
#####################
	freq.a        	= freq.all(hap.pop)


# calcul de la vraisemblance sous H0 (pas de QTL)
#################################################
	 desc.pere     = 	descendant.pere(genea)
	 moyenne.pere  = 	moyenne.pere(perf,CD,desc.pere)
	 ML.H0         = 	obj.LDL.add.H0(moyenne.pere,perf,CD,desc.pere)


# Boucle
########
# Initialisation des vecteurs/matrices résultats:LRT,pos.test,assoc.est,hap.ass.est,param.est
	 LRT		= 	rep(NA,num.pos)
	 pos.test  	= 	rep(NA,num.pos)
	 assoc.est 	= 	rep(NA,num.pos)
	 hap.ass.est	= 	rep(NA,num.pos)
	 if(perfectLD==TRUE){
	      param.est = matrix(NA,ncol=num.pos,nrow=5)
	 }
	 if(perfectLD==FALSE){
	      param.est = matrix(NA,ncol=num.pos,nrow=6)
	 }

# Initialisation de l'index comptant les positions de test
 	i.pos	=	0

# Initialisation de l'index comptant le nbre de colonnes à passer dans proba_de_liaison
 	i.pos.PLA	=	0


#  if(marq.hap.left>1) {
#      for(i in 1:(marq.hap.left-1)) {   
#         nbre.pas	=	length(dist.test[[i]])
#         i.pos.PLA	=	i.pos.PLA    +   nbre.pas
#
#      }
#   } 


    for (i in marq.hap.left:(nbre.int-marq.hap.left+1)) { # Début de la boucle sur les intervalles (FOR1)

    #extraction des informations pour l'intervalle i
    nbre.all.marq	=	rep(0,marq.hap)
    freq.marq		=	list()
    all.marq.int	=	list()

        for(im in 1:(marq.hap)){ # Début FOR2
             nbre.all.marq[im]	=	length(all.marq[[(i-marq.hap.left)+im]])
             freq.marq[[im]]	=	freq.a[[(i-marq.hap.left)+im]]
             all.marq.int[[im]]	=	all.marq[[(i-marq.hap.left)+im]]
        }# Fin FOR2


    #calcul des structures pour l'intervalle i
    res.structure	=	structure.hap(marq.hap,nbre.all.marq)
    pi.hap		=	pi.hap(freq.marq,res.structure)
    cor.chrm1.pere	=	corresp(hap.chrom1.pere[,((i-marq.hap.left)+1):((i-marq.hap.left)+marq.hap)],res.structure)
    cor.chrm2.pere	=	corresp(hap.chrom2.pere[,((i-marq.hap.left)+1):((i-marq.hap.left)+marq.hap)],res.structure)
    cor.mere		=	corresp(hap.trans.mere[,((i-marq.hap.left)+1):((i-marq.hap.left)+marq.hap)],res.structure)
    hap.assoc		=	unique(c(cor.chrm1.pere$assoc,cor.chrm2.pere$assoc,cor.mere$assoc))   
    nbre.ass		= 	length(hap.assoc) 
    nbre.pas		=	length(dist.test[[i]])


        for(h in 1:nbre.pas){  # Début de la boucle sur les pas des intervalles (FOR3)

             for (z in 1:length(position.new)){# début FOR4

                 if (((h==1)&(position.new[z]==dist.test[[i]][h]))|((h>=2)&(position.new[z]==dist.test[[i]][h]))){# début FOR5

                     i.pos	=	i.pos+1
                     depart   	= 	depart.LDL(moyenne.pere,perf,CD,PLA[,(i.pos+i.pos.PLA)],desc.pere)
                     i.dist	=	c(dist.marq[(i-marq.hap.left+1):i],dist.test[[i]][h],dist.marq[(i+1):(i+marq.hap.left)])
                     poids.D	=	poids.D(i.dist,marq.hap,res.structure)


		     # initialisation des vecteurs nécessaires dépendant du nombre d'associations
                     ML		=	rep(NA,nbre.ass)

                     if(perfectLD==TRUE){
                         val.par	=	matrix(NA,nrow=nbre.ass,ncol=5)
                     }

                     if(perfectLD==FALSE){
                         val.par	=	matrix(NA,nrow=nbre.ass,ncol=6)
                     }

    
                    for (j in 1:nbre.ass){ # Début de la boucle sur toutes les associations possibles (FOR6)
                         assoc 	= 	hap.assoc[j]
                         don 	= 	list(perf,CD,PLA[,(i.pos+i.pos.PLA)],desc.pere,moyenne.pere,assoc,cor.chrm1.pere,cor.chrm2.pere,cor.mere,pi.hap,res.structure,poids.D)

                         # recherche du minimum de la fonction obj.LDL.add (inverse de la vraisemblance)
                         if(perfectLD==TRUE){
                              start    	= 	c(depart,temps.depart,0.25,0)
                              op  	= 	optim(start,obj.LDL.add,NULL,method="BFGS",lower=-Inf,upper=Inf,control=list(),hessian=FALSE,don)
                         }

                         if(perfectLD==FALSE){
                              start    	= 	c(depart,temps.depart,0.9,0.25,0)
#                              op  	= 	optim(start,obj.LDL.add.alpha,NULL,method="BFGS",lower=-Inf,upper=Inf,control=list(),hessian=FALSE,don)
                         }
    

                         val.par[j,]   =        op$par  
                         ML[j]         =        -op$value

                     }# Fin de la boucle sur les associations (FOR6)


                   # Recherche de la valeur du maximum de vraisemblance  pour chaque position
                   # sur les association possibles.
                   # Calcul des valeurs du test et stockage des paramètres estimés

                   pos.test[i.pos] 	= 	dist.test[[i]][h]
                   LRT[i.pos]           = -2*(ML.H0- max(ML))
                   i.val                =	which.max(ML)
                   assoc.est[i.pos]     = 	hap.assoc[i.val]
                   param.est[,i.pos]	=	val.par[i.val,]

		   # retrouver les allèles de l'haplotype associé
                    hap.ass.est[i.pos]=retrouve.all(assoc.est[i.pos],res.structure,all.marq.int)

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

