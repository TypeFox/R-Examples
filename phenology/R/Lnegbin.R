
#' .Lnegbin estimate a negative binomial likelihood.
#' @title The function ".Lnegbin"
#' @author Marc Girondot
#' @return Return the likelihood
#' @param x Set of parameters to be fitted
#' @param pt Transfer parameters
#' @description Function of the package phenology

.Lnegbin <- function(x, pt) {
  
  # 19/3/2016: je rajoute cofactors et add.cofactors dans pt
  
  #  if (length(pt)==1 & names(pt)[1]=="pt") pt <- pt$pt
  # .phenology.env<- NULL
  # rm(.phenology.env)
  
  # pt=list(data=data, fixed=parametersfixed, incertitude=method_incertitude, zerocounts=zero_counts)
  
  sum=0
  # je mets tous les paramètres dans xpar
  xpar <- c(x, pt$fixed)
  out <- pt$out
  
  datatot <- pt$data
  infinite <- pt$infinite
  daily_count <- getFromNamespace(".daily_count", ns="phenology")
  format_par <- getFromNamespace(".format_par", ns="phenology")
  
  for(k in seq_along(datatot)) {
    
    # print(paste("pop", k))
    data <- datatot[[k]]
    # quel est le nom de la série en cours
    nmser <- names(datatot)[k]
    # Prend en compte les 0 ou non 5/2/2012
    zero <- pt$zerocounts[k]
    deb <- ifelse(zero, 0, 1)
    
    # je fais une fonction où j'envoie les paramètres et la série et il renvoie ceux à utiliser directement
    xparec <- format_par(xpar, nmser)
    #  if (xparec["MinE"]<0) print(xpar, xparec)
    th <- xparec["Theta"]
    
    for (i in 1:dim(data)[1]) 
    {

      # est ce que j'ai une observation ?
      # transformé en 0 ou non avec zero TRUE ou FALSE: 5/2/2012
      if ((data$nombre[i]!=0)||zero) {
        
        if (is.na(data$Date2[i])) {
          #######################
          # on a une seule date #
          #######################
          
          # ici je devrais aussi faire le calcul en conditionnel
          sumnbcount <- daily_count(data$ordinal[i], xparec, print=FALSE, zero=pt$zero)
          
          # 19/3/2016: je dois inclure les cofacteurs
          if (!is.null(pt$add.cofactors)) {
 #           print("Nom des cofacteurs")
#            print(pt$add.cofactor)
#            print("Date")
#            print(data$Date[i])
#            print("Paramètres")
#            print(xparec[pt$add.cofactor])
#            print("Cofacteurs")
#            print(pt$cofactors[ pt$cofactors$Date==data$Date[i], pt$add.cofactor] * xparec[pt$add.cofactor])
#            print("Nombre avant cofacteurs")
#            print(sumnbcount)
            sumnbcount <- sumnbcount + sum(pt$cofactors[ pt$cofactors$Date==data$Date[i], pt$add.cofactor] * xparec[pt$add.cofactor])
#            print("Nombre après cofacteurs")
#            print(sumnbcount)
            if (sumnbcount <= pt$zero) sumnbcount <- pt$zero
          }
          
          # dans zero j'ai l'information si je prends les 0 ou non pour cette série
          if (!zero) {
            lnli2 <- -log(dnbinom(data$nombre[i], size=th, mu=sumnbcount, log=FALSE)/(1-dnbinom(0, size=th, mu=sumnbcount, log=FALSE)))
          } else {
            lnli2 <- -dnbinom(data$nombre[i], size=th, mu=sumnbcount, log=TRUE)
          }
          
        } else {
          ###################
          # on a deux dates #
          ###################
          
          # nombre de jours. On met des NA  
          nbjour <- data$ordinal2[i]-data$ordinal[i]+1
          #je met le nombre théorique de ce jour dans nbcount
          nbcount <- daily_count((1:nbjour)+data$ordinal[i]-1, xparec, print=FALSE)
          
          #je somme
          sumnbcount <- sum(nbcount)
          
          if (pt$incertitude==1 | pt$incertitude=="convolution") {
            
            # je suis sur la méthode 1 d'incertitude
            if (!zero) {
              lnli2 <- -log(dSnbinom(data$nombre[i], size=th, mu=nbcount, log=FALSE, infinite=infinite)/(1-dSnbinom(0, size=th, mu=nbcount, log=FALSE, infinite=infinite)))
            } else {
              lnli2 <- -dSnbinom(data$nombre[i], size=th, mu=nbcount, log=TRUE, infinite=infinite)
              # if (!is.finite(lnli2)) print(dput(list(x=data$nombre[i], size=th, mu=nbcount)))
            }
            
            
          } 
          
          if (pt$incertitude==2 | pt$incertitude=="combinatory") {
            
            # je suis sur la méthode 2 d'incertitude
            # nbcount est le nombre théorique
            nbcountrel <- nbcount/sumnbcount
            nbcountrel[nbcountrel==0] <- pt$zero
            nbcountrel[nbcountrel==1] <- 1-pt$zero
            
            #dans nbcountrel j'ai les probabilités
            #dans nbjour j'ai le nombre de jours
            #dans data$nombre[i] j'ai le nombre de nids observés
            
            #Indistinguishable Objects to Distinguishable Boxes
            # http://2000clicks.com/mathhelp/CountingObjectsInBoxes.aspx
            # The number of different ways to distribute n indistinguishable balls into
            # k distinguishable boxes is C(n+k-1,k-1).
            # For example, 5 balls into 3 boxes can be done in these C(7,2) = 21 ways:
            
            N <- data$nombre[i]
            
            # The number of different ways to distribute n indistinguishable balls into
            # k distinguishable boxes is C(n+k-1,k-1).
            # nb<-choose(N+nbjour-1,nbjour-1)=dim(tb)[1]
            # divers<-matrix(rep(0, nbjour*nb), ncol=nbjour)
            
            # generate all possible positions of the boundaries
            xx <- combn(N+nbjour-1, nbjour-1)
            # compute the number of balls in each box
            a <- cbind(0, diag(nbjour)) - cbind(diag(nbjour), 0)
            tb <- t(a %*% rbind(0, xx, N+nbjour) - 1)
            
            
            # je calcule déjà les dbnbinom pour toutes les solutions: 4/2/2012
            # il me faut un tableau de 1:nbjour et de 0:N
            a <- try(matrix(rep(0, (N+1)*nbjour), nrow=N+1, ncol=nbjour), silent=TRUE)
            
            if (class(a)=="try-error") {
              stop("Too many incertitudes on the days. Use the other method named 'convolution'.")
            }
            
            # Je calcule la matrice des vraisemblances
            # de toutes les solutions pour éviter de les recalculer à chaque fois
            
            for(ii in deb:N) {
              for(countday in 1:nbjour) {
                if (!zero) {
                  a[ii+1, countday] <- log(dnbinom(ii, size=th, mu=nbcount[countday], log=FALSE)/(1-dnbinom(0, size=th, mu=nbcount[countday], log=FALSE)))
                } else {
                  a[ii+1, countday] <- dnbinom(ii, size=th, mu=nbcount[countday], log=TRUE)
                }
              }
            }
            
            
            sump <- 0
            for(ii in 1:dim(tb)[1]) {
              p <- dmultinom(tb[ii,1:nbjour], prob=nbcountrel, log=TRUE)
              for(countday in 1:nbjour) {
                p <- p + a[tb[ii,countday]+1, countday]
              }
              sump <- sump + exp(p)
            }
            lnli2 <- -log(sump)
          }
          
          
          # fin du test if ((is.na(data$Date2[i]))) {
        }
        
        # transformé en test de zero
      } else {
        # 19/3/2016: Si data$nombre[i]==0 alors lnli2 gardait la valeur d'avant
        lnli2 <- NA
        sumnbcount <- NA
      }
      
      datatot[[k]]$LnL[i] <- lnli2
      datatot[[k]]$Modeled[i] <- sumnbcount
      # sum <- sum + lnli2
      
      # fin de la boucle des jours
    }
    
    sum <- sum + sum(datatot[[k]]$LnL, na.rm = TRUE)
    
    # fin de la boucle des séries
  }
  if (out) {
    return(sum)
    } else {
    return(datatot)
  }
  # fin de la fonction
}
