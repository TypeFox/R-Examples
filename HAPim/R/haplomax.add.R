`haplomax.add` <-
function(hap.trans.pere,hap.trans.mere,perf,CD,map,marq.hap){

# vérification des dimensions:
#############################
     if(dim(hap.trans.mere)[2] != dim(hap.trans.pere)[2]  | dim(hap.trans.pere)[2] != (length(map)+1) ) stop("the numbers of genotype information are not consistent in  hap.trans.mere, hap.trans.pere or map ",call.=FALSE)

     if(length(CD) != length(perf) | length(perf) != dim(hap.trans.mere)[1] | dim(hap.trans.mere)[1] != dim(hap.trans.pere)[1]  ) stop("the numbers of individual information are not consistent in CD, perf, hap.trans.pere or hap.trans.mere",call.=FALSE)




# recodage des haplotypes et des genotypes
##########################################
	 hap.pop   	=	rbind(hap.trans.pere,hap.trans.mere)
	 all.marq  	=	allele.marq(hap.pop)
	 hap.trans.mere = 	recode.hap(hap.trans.mere,all.marq)
	 hap.trans.pere = 	recode.hap(hap.trans.pere,all.marq)
	 hap.pop  	=	rbind(hap.trans.pere,hap.trans.mere)

# On calcule la distance des marques par rapport à l'origine
############################################################
	 dist.marq	=	distance.marqueurs(map)

# On calcule le nombre d'intervalles entre les marques
######################################################
	 nbre.int	=	length(map) 


# On calcule du nombre de positions de tests
############################################
	  num.pos	= 	nbre.int-marq.hap+2




#initialisation des vecteurs résultats
######################################
	 Fisher		=	rep(NA,num.pos)
	 assoc.est 	= 	rep(NA,num.pos) 
	 hap.ass.est	= 	rep(NA,num.pos)
	 param.est 	= 	matrix(NA,ncol=num.pos,nrow=2)

for (i in 1:num.pos)

        {  # Début de la boucle sur les positions de tests

              #extraction des informations pour la position  i
              nbre.all.marq	=	rep(0,marq.hap)
              all.marq.int	=	list()
              for(im in 1:(marq.hap)){
                 nbre.all.marq[im]	=	length(all.marq[[(i+im-1)]])
                 all.marq.int[[im]]	=	all.marq[[(i+im-1)]]
              }

              #calcul des structures pour la position i
               res.structure	=	structure.hap(marq.hap,nbre.all.marq)
               cor.hap		=	corresp(hap.pop[,i:(i+marq.hap-1)],res.structure)
               pi.hap		=	pi.hap.NI(res.structure,cor.hap)
               cor.pere		=	corresp(hap.trans.pere[,i:(i+marq.hap-1)],res.structure)
               cor.mere		=	corresp(hap.trans.mere[,i:(i+marq.hap-1)],res.structure)  
               hap.assoc	=	unique(c(cor.pere$assoc,cor.mere$assoc))   
               nbre.ass		=	length(hap.assoc) 

              # initialisation des vecteurs nécessaires dépendant du nombre d'associations
              F.assoc	=	rep(NA,nbre.ass)
              val.par	=	matrix(NA,nrow=nbre.ass,ncol=2)
 
              # Début de la boucle sur toutes les associations possibles:
              if(nbre.ass>1) { # début du if sur nbre.ass
                   for (j in 1:nbre.ass){  # Début de la boucle sur toutes les associations possibles
                      assoc		=	hap.assoc[j]  
                      res		=	obj.haplomax.add(perf,CD,assoc,res.structure,pi.hap,cor.pere,cor.mere)
                      sum.res	=	summary(res)
                      temp		=	sum.res[[1]]
                      F.assoc[j]	=	temp$F[1]
                      val.par[j,2]	=	res$coefficient[2]
                      val.par[j,1]	=	temp$Mean[length(temp$Mean)]
               }  # Fin de la boucle sur les associations

              i.val=which.max(F.assoc)
              if(!is.na(i.val)){
              assoc.est[i]	=	hap.assoc[i.val]
              Fisher[i]	=	F.assoc[i.val]
              param.est[,i]	=	val.par[i.val,]

	      #retrouver les allèles de l'haplotype associé
              hap.ass.est[i]	=	retrouve.all(assoc.est[i],res.structure,all.marq.int)
              } 

             }# fin du if sur nbre.assoc

         }# fin de la boucle sur les positions de tests

#On calcul les positions de test (milieu de l'haplotype)
###################################
    nbre.marq 	= 	length(dist.marq)
    temp1=dist.marq[1:(nbre.marq-marq.hap+1)]
    temp2=dist.marq[marq.hap:nbre.marq]
    pos.test=temp1+(temp2-temp1)/2


#Regroupement des résultats sous forme de tableau
#################################################
    res		=	data.frame(pos.test,Fisher,hap.ass.est,t(param.est)) 

#On nomme les colonnes du data.frame
dimnames(res)[[2]]=c("position","F","haplotype","variance","effect.Q")

res

}

