`descendant.pere` <-
function(genea){

  iden.pere = unique(genea[,2])
  nbre.pere = length(iden.pere)

# Initialisation des vecteur/matrice de résultat
  desc.pere 	 = 	matrix(NA,ncol=2,nrow=nbre.pere)
  nbre.desc.pere = 	rep(NA,nbre.pere)


  for(ip in 1:nbre.pere){

    nbre.desc.pere[ip] = sum(genea[,2]==iden.pere[ip])

  }


  desc.pere[1,1] = 1
  desc.pere[1,2] = nbre.desc.pere[1]
 
  if(nbre.pere>1){

      for(ip in 2:nbre.pere){

          desc.pere[ip,1]=desc.pere[(ip-1),2]+1
          desc.pere[ip,2]=desc.pere[ip,1]+nbre.desc.pere[ip]-1 

      }
  }
  
  desc.pere


}

