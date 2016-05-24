`moyenne.pere` <-
function(perf,CD,desc.pere){

  nbre.pere	=	length(desc.pere[,1])

  #Initialisation du vecteur de résultat
  moyenne.pere	=	rep(0,nbre.pere)


  Y 		= 	perf*CD*CD
  CDsquare	=	CD*CD


  for(ip in 1:nbre.pere){

    deb   =  desc.pere[ip,1]
    fin   =  desc.pere[ip,2]

    if(sum(CDsquare[deb:fin])>0.00000001){

         moyenne.pere[ip]=sum(Y[deb:fin])/sum(CDsquare[deb:fin])
    }

  }

  moyenne.pere


}

