dependogram <- function(X,vecd.ou.p,N=10,B=2000,alpha=0.05,display=TRUE,graphics=TRUE,nbclus=1) {


  if (nbclus>1) {

    suppressWarnings(pkg.present <- require(snow))
    if (!pkg.present) stop("Package snow is not installed!")
    suppressWarnings(pkg.present <- require(rsprng))
    if (!pkg.present) stop("Package rsprng is not installed!")
    suppressWarnings(pkg.present <- require(Rmpi))
    if (!pkg.present) stop("Package Rmpi is not installed!")
    
  }


  X <- as.matrix(X)

# si length(vecd.ou.p)>1 alors cas non sériel sinon cas sériel

  if (length(vecd.ou.p) > 1) {
# on fait le cas non sériel
    seriel <- 0
    vecd <- vecd.ou.p
    p <- length(vecd)
    taille <- 2^p-p-1
#


    
    RnAs <- rep(0,2^p-p-1)
    Rn <- 0

# On charge la fonction C dans la mémoire
#	     dyn.load(paste("dependogram", .Platform$dynlib.ext, sep=""))

# Remarque: quand on passe une matrice dans la fonction .C elle est reçue comme un vecteur obtenu en concaténant les colonnes de cette matrice
# et le résultat est lui aussi renvoyé sous la forme d'un tel vecteur




# On démarre la grappe de calculs
    if (nbclus > 1) {

      B <- round(B/nbclus)*nbclus
      

      cl <- makeCluster(nbclus, type = "MPI") 
      clusterSetupSPRNG(cl)
                                        
      
      myfunc <- function(B,p) {
        RnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
        Rnstar <- rep(0,B)
        require(IndependenceTests)
# On appelle la fonction C dependogram
        .C("dependogram",
           as.integer(N),
           as.integer(vecd),
           as.integer(length(vecd)),
           as.integer(p),
           as.numeric(X),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           B=as.integer(B),
           as.numeric(alpha),
           RnAs=as.numeric(RnAs),
           RnAsstar=as.numeric(RnAsstar),
           Rn=as.numeric(Rn),
           Rnstar=as.numeric(Rnstar),
           as.integer(seriel),PACKAGE="IndependenceTests")
        
      }
      
      out2 <- clusterCall(cl, myfunc, round(B/nbclus),p)
      
                   
      
      # On arrete la grappe de calcul
      stopCluster(cl)

      out <- list(RnAs=c(),RnAsstar=c(),Rn=c(),Rnstar=c())


      for (clus in 1:nbclus) {
        
        out$Rnstar <- c(out$Rnstar,out2[[clus]]$Rnstar)
        out$RnAsstar <- c(out$RnAsstar,out2[[clus]]$RnAsstar)
        
      }

      out$Rn <- out2[[clus]]$Rn
      out$RnAs <- out2[[clus]]$RnAs
            
      
    } else {

      RnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
      Rnstar <- rep(0,B)

# On appelle la fonction C dependogram
      out <- .C("dependogram",
		as.integer(N),
		as.integer(vecd),
		as.integer(length(vecd)),
		as.integer(p),
		as.numeric(X),
		as.integer(nrow(X)),
		as.integer(ncol(X)),
		as.integer(B),
		as.numeric(alpha),
		RnAs=as.numeric(RnAs),
		RnAsstar=as.numeric(RnAsstar),
		Rn=as.numeric(Rn),
		Rnstar=as.numeric(Rnstar),
		as.integer(seriel),PACKAGE="IndependenceTests")
      
    }




# On décharge la fonction C de la mémoire
#		dyn.unload(paste("dependogram", .Platform$dynlib.ext, sep=""))



    if (display) {
# require(combinat)
      RES <- as.list(1:(p-1))
      for (cardA in 2:p) {RES[[cardA]] <- as.matrix(combn(p,cardA))}
      nb <- 0
      for (cardA in 2:p) {
        for (j in 1:(choose(p,cardA))) {
          nb <- nb+1
          
          cat(c(nb,": A=",RES[[cardA]][,j],": ||RnA||=",round(out$RnAs[nb],3),"\n"))
          
          
        }
      }
    }
    

# Ordonne les éléments de chaque ligne
    RnAsstar <- matrix(out$RnAsstar,nrow=(2^p-p-1),ncol=B,byrow=FALSE)
    matseuils <- t(apply(RnAsstar,FUN=sort,MARGIN=1))


    beta <- (1-alpha)^(1/taille)
# Contient les beta-quantiles des stats RnA pour chacun des des 2^p-p-1 ensembles A
    AllThresholds <- matseuils[,round(beta*B)]

    if (graphics) {
# On trace une barre verticale pour chaque A de hauteur ||RnA||
      plot(out$RnAs,type="h",ylim=c(0,max(c(max(AllThresholds),max(out$RnAs),max(out$Rnstar)))+0.1),xlim=c(0,2^p-p),main="Dependogram",xlab="Subsets",ylab="||RnA||")
    }

    Rnstar <- sort(out$Rnstar)
    GlobalThreshold <- Rnstar[round((1-alpha)*B)]
# abline(h=GlobalThreshold,col="red")


    if (graphics) { 
# On met une étoile pour chaque beta-quantile de ||R_A||
      points((1:(2^p-p-1)),AllThresholds,pch="*")
    }

    res <- list(norm.RnA=out$RnAs,Rn=out$Rn,rA=AllThresholds,r=GlobalThreshold,RnAsstar=RnAsstar)
    return(res)


  }



  if (length(vecd.ou.p) == 1) {
# on fait le cas sériel
    seriel <- 1
    p <- vecd.ou.p
    vecd <- rep(ncol(X),p)
    taille <- 2^(p-1)-1




    SnAs <- rep(0,taille)
    Sn <- 0

# On charge la fonction C dans la mémoire
#	     dyn.load(paste("dependogram", .Platform$dynlib.ext, sep=""))

# Remarque: quand on passe une matrice dans la fonction .C elle est reçue comme un vecteur obtenu en concaténant les colonnes de cette matrice
# et le résultat est lui aussi renvoyé sous la forme d'un tel vecteur




# On démarre la grappe de calculs
    if (nbclus > 1) {
      
      B <- round(B/nbclus)*nbclus

      SnAsstar <- matrix(0,nrow=taille,ncol=B)
      Snstar <- rep(0,B)
      
      cl <- makeCluster(nbclus, type = "MPI") 
      clusterSetupSPRNG(cl)
                                        
      
      myfunc <- function(B,p) {
        SnAsstar <- matrix(0,nrow=(2^p-p-1),ncol=B)
        Snstar <- rep(0,B)
        require(IndependenceTests)
# On appelle la fonction C dependogram
        .C("dependogram",
           as.integer(N),
           as.integer(vecd),
           as.integer(length(vecd)),
           as.integer(p),
           as.numeric(X),
           as.integer(nrow(X)),
           as.integer(ncol(X)),
           B=as.integer(B),
           as.numeric(alpha),
           SnAs=as.numeric(SnAs),
           SnAsstar=as.numeric(SnAsstar),
           Sn=as.numeric(Sn),
           Snstar=as.numeric(Snstar),
           as.integer(seriel),PACKAGE="IndependenceTests")
      }
      
      out2 <- clusterCall(cl, myfunc, round(B/nbclus),p)
      
                                        
      
      # On arrete la grappe de calcul
      stopCluster(cl)

      out <- list(SnAs=c(),SnAsstar=c(),Sn=c(),Snstar=c())
      

      for (clus in 1:nbclus) {
        
        out$Snstar <- c(out$Snstar,out2[[clus]]$Snstar)
        out$SnAsstar <- c(out$SnAsstar,out2[[clus]]$SnAsstar)
        
      }

      out$Sn <- out2[[clus]]$Sn
      out$SnAs <- out2[[clus]]$SnAs

      
    } else {

      SnAsstar <- matrix(0,nrow=taille,ncol=B)
      Snstar <- rep(0,B)
# On appelle la fonction C dependogram
      out <- .C("dependogram",
		as.integer(N),
		as.integer(vecd),
		as.integer(length(vecd)),
		as.integer(p),
		as.numeric(X),
		as.integer(nrow(X)),
		as.integer(ncol(X)),
		as.integer(B),
		as.numeric(alpha),
		SnAs=as.numeric(SnAs),
		SnAsstar=as.numeric(SnAsstar),
		Sn=as.numeric(Sn),
		Snstar=as.numeric(Snstar),
		as.integer(seriel),PACKAGE="IndependenceTests")
      
    }





# On décharge la fonction C de la mémoire
#		dyn.unload(paste("dependogram", .Platform$dynlib.ext, sep=""))


    if (display) {

# require(combinat)
      RES <- as.list(1:(p-1))
      for (cardA in 2:p) {RES[[cardA]] <- as.matrix(rbind(rep(1,choose(p-1,cardA-1)),as.matrix(combn(p-1,cardA-1)+1)))}
      nb <- 0
      for (cardA in 2:p) {
        for (j in 1:(choose(p-1,cardA-1))) {
          nb <- nb+1

          cat(c(nb, ": A=",RES[[cardA]][,j],": ||SnA||=",round(out$SnAs[nb],3),"\n"))

        }
      }
    }


    SnAsstar <- matrix(out$SnAsstar,nrow=taille,ncol=B,byrow=FALSE)
    Sn <- max(out$SnAs)


    beta <- (1-alpha)^(1/taille)


# Le beta-quantile de S_n,A est calculé en amalgamant toutes les
# valeurs  S_n,A^* avec |A|=k comme dans l'article
# Il y a choose(p-1,|A|-1) ensembles A de taille |A| qui contiennent 1
# Il faut donc prendre dans la matrice SnAsstar des paquets de choose(p-1,|A|-1) lignes, pour |A|=2 to p.
# Il y aura p-1 paquets.
# Pour chacun de ces paquets, on crée un vecteur vecA en prenant tous les éléments du paquet, 
# puis on calcule AllThresholds[|A|]<-vecA[round(beta*B*choose(p-1,|A|-1))] pour |A|=2 to p
# Ce vecteur AllThresholds contiendra donc les hauteurs des barres horizontales à placer sur le dependogram (une barre horizontale 
# pour chaque |A|)


    AllThresholds <- rep(0,p-1)
    begin <- 1
    end <- 0
    for (cardA in 2:p) {
      end <- end+choose(p-1,cardA-1)
      vecA <- as.vector(SnAsstar[begin:end,])
      vecA <- sort(vecA)
      AllThresholds[cardA-1] <- vecA[round(beta*B*choose(p-1,cardA-1))]
      begin <- end+1
    }

    if (graphics) {
# On trace une barre verticale pour chaque A de hauteur ||SnA||
      plot(out$SnAs,type="h",ylim=c(0,max(c(max(AllThresholds),max(out$SnAs),max(out$Snstar)))+0.1),xlim=c(0,2^(p-1)),main="Dependogram",xlab="Subsets",ylab="||SnA||")
    }
    
    Snstar <- sort(out$Snstar)
    GlobalThreshold <- Snstar[round((1-alpha)*B)]
# abline(h=GlobalThreshold,col="red")


# Il reste à placer les lignes horizontales des seuils critiques pour chaque |A|


    if (graphics) {
      begin <- 1
      end <- 0
      for (cardA in 2:p) {
        end <- end+choose(p-1,cardA-1)

        segments(begin-0.5,AllThresholds[cardA-1],end+0.5,AllThresholds[cardA-1],lty=4)

        begin <- end+1
      }
    }





    res <- list(norm.SnA=out$SnAs,Sn=out$Sn,sA=AllThresholds,s=GlobalThreshold,SnAsstar=SnAsstar)
    return(res)
  }
  


}



