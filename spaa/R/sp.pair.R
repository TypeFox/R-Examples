sp.pair <-
function(matr)
{  
   if(any(is.na(matr))){
      matr <- na.omit(matr)
      cat("The data matrix contains NA, and have been removed.\n")
   }
   pearson  <- cor(matr, method = "pearson")
   spearman <- cor(matr, method = "spearman")
   matr1 = matr
   N <- nrow(matr1) ##number of plots
   matr1[matr1>1] <- 1
   n <- apply(matr1, 2, sum) ##number of species for each plots
   
   tmatr <- t(matr1)
   df <- as.matrix(tmatr)
   
   a <- df %*% t(df)          ## Number of plots Presence for both speciesA and  speciesA
   b <- df %*% (1 - t(df))    ## Number of plots Presence for speciesA
   c <- (1 - df) %*% t(df)    ## Number of plots Presence for speciesB
   d <- ncol(df) - a - b - c  ## Number of plots absence for both speciesA and  speciesA
   ##Yate's correction for chisq
   chisq <- data.frame()
   V <- data.frame()
   Ochiai <- data.frame()
   Dice <- data.frame()
   Jaccard <- data.frame()
   PCC <- data.frame() 
   for (i in 1:nrow(a)){
       for (j in 1:ncol(a)){
           chisq[i,j]   <- ((((abs(a[i,j] * d[i,j] - b[i,j]*c[i,j]) 
                           - 0.5*N)^2)*N )/((a[i,j]+b[i,j])*(a[i,j]+c[i,j])*(b[i,j]
                           +d[i,j])*(c[i,j]+d[i,j])))
           V[i,j]       <- ((a[i,j]+d[i,j])-(b[i,j]+c[i,j]))/(a[i,j]+b[i,j]+c[i,j]+d[i,j]) # V>0,positive, V<0,negative # Thanks for Ms Xueni Zhang for pointing out the error.
           Ochiai[i,j]  <- a[i,j]/sqrt((a[i,j]+b[i,j])*(a[i,j]+c[i,j])) #Ochiai index
           Dice[i,j]    <- 2*a[i,j] /(2*a[i,j] + b[i,j] + c[i,j]) #Dice index
           Jaccard[i,j] <- a[i,j]/(a[i,j]+b[i,j]+c[i,j]) #Jaccard index
           PCC[i,j]     <- (a[i,j] * d[i,j] - b[i,j] * c[i,j])/((a[i,j] + b[i,j])*(a[i,j] 
                            + c[i,j])*(c[i,j] + d[i,j])*(b[i,j] + d[i,j])) ##Percentage cooccurance 
       }
    }

   chisq <- data.frame(chisq)
   dd <- data.frame() 
   for (i in 1:nrow(chisq)){
       for (j in 1:ncol(chisq)){ 
            if (chisq[i, j] > 6.635){ 
		        chisq[i, j] <- paste(chisq[i, j], "**") 
				}
            if (chisq[i, j] < 6.635& chisq[i, j] > 3.841){ 
		        chisq[i, j] <- paste(chisq[i, j], "*")
				}
        dd[i,j] <- a[i,j] * d[i,j] - b[i,j] * c[i,j]
        }
   }
   AC = data.frame() ## Association Coefficients
   for (i in 1: nrow(a)){
       for (j in 1: ncol(a)){
           if (a[i,j] * d[i,j] >= b[i,j] * c[i,j]){ 
		       AC[i,j] <- (a[i,j] * d[i,j] - b[i,j] * c[i,j])/
			              ((a[i,j] + b[i,j])*(b[i,j] + d[i,j])) 
			   }
           if ( a[i,j]*d[i,j] < b[i,j]*c[i,j] & d[i,j] >= a[i,j]){ 
		      AC[i,j] <- (a[i,j] * d[i,j] - b[i,j] * c[i,j])/
			             ((a[i,j] + b[i,j])*(a[i,j] + c[i,j])) 
			  }
    
           if (b[i,j]*c[i,j] > a[i,j] * d[i,j] & d[i,j] < a[i,j] ){ 
		      AC[i,j] <- (a[i,j] * d[i,j] - b[i,j] * c[i,j])/
			             ((b[i,j] + d[i,j])*(d[i,j] + c[i,j]))}
       }
   }
   result <- list(chisq=chisq, chisqass = dd, V=V, Ochiai=Ochiai, 
             Dice=Dice, Jaccard=Jaccard, Pearson = pearson, 
			 Spearman = spearman, PCC=PCC, AC=AC )
   return(result)
}

