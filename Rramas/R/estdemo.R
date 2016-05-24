estdemo <-
function(v0,mat,stmat=NULL, fecundity1=TRUE){
   # fecundity (B) is picked from a poisson distribution
   # survival (S) is picked from a binomial distribution
   v <- v0 # TODO: change v for v0 in all the function
   B <- NULL  
   S <- NULL  # survival
   
   # option 1: all data in first row are assumed to be fecundities; data in the
   #   rest of rows are assumed to represent survival probabilities if <=1 and 
   # fecundities if >1
   
   if (is.null(stmat) & fecundity1==TRUE){
      for (i in 1: (dim(mat)[1])){
        # for each individual in each stage (vi) asign descendence 
        # based in a poisson distribution with mean mat[1,i]
        # (only if mat[1,i] >0)
        Bi <- 0
        if(mat[1,i] >0) Bi <- rpois(v[i],mat[1,i])
        B <- c(B, Bi) 
        Sj <- NULL
        for (j in 1:(dim(mat)[1]-1)){
          # for each individual in each stage (v[i]) asign survival or descendence
          # if mat[j+1,i]>1 we assume that this is reproduction and asign
          # descendence based in a poisson distribution with mean mat[j+1,i].
          # if mat[j+1,i]<=1 we assume that this is survival and assign to each 
          # individual in each stage (v[i])reproduction and asign survival based in a
          # binomila distribution with probability  mat[j+1,i]
    
          Sij <- 0 
          #only sample rpois or rbinom if transition-probabilities >0
          if(mat[j+1,i] > 0){
              if(mat[j+1,i]>1) Sij <- sum(rpois(v[i],mat[j+1,i])) 
                else Sij <- sum(rbinom(v[i],size=1,mat[j+1,i]))
          }   
          Sj <- c(Sj,Sij)  
        }
        S <- cbind(S,Sj)     
       }
       B <- sum(B)
       #B0 <- sum(B)
       S <- apply(S,1,sum)
       #S0<- apply(S,1,sum)
       result <- c(B,S)
       #result <- c(B0,S0)
   }
   
   # option 2: all data > 1 are assumed to represent
   # fecundities; all data <=1 are assumed to represent survival
   
   if (is.null(stmat) & fecundity1!=TRUE){
      for (i in 1: (dim(mat)[1])){
         Sj <- NULL
        for (j in 1:(dim(mat)[1])){
          # for each individual in each stage (v[i]) asign survival or descendence
          # if mat[j,i]>1 we assume that this is reproduction and asign
          # descendence based in a poisson distribution with mean mat[j,i].
          # if mat[j,i]<=1 we assume that this is survival and assign to each 
          # individual in each stage (v[i])reproduction and asign survival based in a
          # binomila distribution with probability  mat[j+1,i]
          Sij <- 0
          #only sample rpois or rbinom if transition-probabilities >0
          if(mat[j,i] > 0){ 
            if(mat[j,i]>1) Sij <- sum(rpois(v[i],mat[j,i])) 
               else Sij <- sum(rbinom(v[i],size=1,mat[j,i]))
          }
          Sj <- c(Sj,Sij)  
        }
        S <- cbind(S,Sj)     
       }
       S <- apply(S,1,sum)
       result <- S
   }
   
   # option 3: accept matrix "stmat" that gives proportion of probability
   # transitions that must be considered fecundities
   
   if (!is.null(stmat) ){
    #real fecundity probabilities
    matF <- mat * stmat
    #real survival probabilities  
    matS <- mat-matF
    for (i in 1: (dim(mat)[1])){
      BSj <- NULL
      for (j in 1:(dim(mat)[1])){
        # for each individual in each stage (v[i]) asign survival or descendence
        # if mat[j,i]>1 we assume that this is reproduction and asign
        # descendence based in a poisson distribution with mean mat[j,i].
        # if mat[j,i]<=1 we assume that this is survival and assign to each 
        # individual in each stage (v[i])reproduction and asign survival based in a
        # binomila distribution with probability  mat[j+1,i]
        Bij <- ifelse(matF[j,i] >0, sum(rbinom(v[i],size=1,matF[j,i])), 0)
        Sij <- ifelse(matS[j,i] >0, sum(rpois(v[i], matS[j,i])), 0) 
        BSij <-  Bij+Sij #sum the contributions of fertility and survival
        BSj <- c(BSj, BSij)  
      }
      S <- cbind(S,BSj)     
    }
    S <- apply(S,1,sum)
    result <- S
   }
   # TODO: option 4: as option 3 but fecundities weighted by survivals, i.e.,
   # those individuals that are asigned a 0 from rbinomial lose all offspring 
   # produced by rpois(or a percentage of these offspring)
   
   return(result)
}
