multiresultm <- function(n,T,F,varF=NULL)
{
  # First, determine numbers from survival and growth
          clas  <- length(n)                            # number of classes
          death <-  1 - colSums(T)                      # vector of death rates
          T     <- rbind(T,death)                       # append death rates to matrix T
          outcome <- matrix(0,nrow(T),clas)             # initialize matrix "outcome"
## see Box 8.11 in Morris and Doak
          for (j in 1:clas) {
            Tj  <- T[,j]                        # extract column j of in matrix T
            ni  <- n[j]                         # extract entry j of vector n
            ind <- matrix(1,1,ni)
            pp  <- cumsum(Tj/sum(Tj))           # find cumulative probabilities
            rnd <- runif(ni,min=0,max=1)		    # make a uniform random for each individual
            for (ii in 1:length(Tj)) {          # find each individual's random fate
                ind <- ind + (rnd > pp[ii])
                } # end for ii
            for (ii in 1:length(Tj)) {          # add up the individuals in each fate
                outcome[ii,j] <- sum(ind == ii)
                } # end for ii
          } #end for j

    # Second, determine numbers of offspring
      if (length(varF) != 0) { # if there is a matrix of inter-individual
                               # variance of fertilities
    # This routine generates lognormal numbers from mean fertilities and their variance.

          offspring <- matrix(0,clas,clas)       # initialize matrix "offspring"
          for (j in 1:clas) {
            fj  <- F[,j]                          # extract column j of matrix F
            if (max(fj) > 0) {                    # skip loop if there is no fertility
              ni  <- n[j]                         # extract entry j of vector n
              for (i in 1:length(fj)) {
                if (F[i,j] > 0) {          # skip if fertility is null
                 
                    # rndfert  <- lnorms(F[i,j],varF[i,j],rnorm(ni,0,1)) # make lognormal random fertilities
                      rndfert  <- lnorms(ni, F[i,j],varF[i,j])   # updated lnorms in version 2.0
                  offspring[i,j] <- sum(rndfert)   # computes number of offsprings from fertilities
                } # end if
              } # end for i
            } # end if
          } # end for j
      } # end if      
      else {
      # This routine generates binomial random births from mean fertilities,
      # for population with clutch size = 1.
          offspring <- matrix(0,clas,clas)        # initialize matrix "offspring"
          for (j in 1:clas) {
            fj  <- F[,j]                          # extract column j of matrix F
            if (max(fj) > 0) {                    # skip loop if there is no fertility
              ni  <- n[j]                         # extract entry j of vector n
              for (i in 1:length(fj)) {
                if (F[i,j] > 0) {                  # skip if fertility is null
                  rndbirth  <- rbinom(ni,1,F[i,j]) # make binomial random births
                  offspring[i,j] <- sum(rndbirth)  # computes number of offspring from random births
                } # end if
              } # end for i
            } # end if
          } # end for j
      } # end else

      multiresultm <- matrix(rowSums(offspring) + rowSums(outcome[1:clas,]),clas,1, dimnames=list(rownames(F), "t+1"))
       
      multiresultm

 } # end of function

