# Version: 30-11-2012, Daniel Fischer

# Changes:
# Added the CSubmat option here , 29-11-2012, DF
# Changed the stop messages, 29-06-2013, DF

  PE1 <- function(X,g,t,goi,alg){

      result <- 0 

    # Elements in group t
      Nt <- sum(g==t)
    # Total amount of elements
      N <- length(g)
    # Set of remaining group labels   
      tprime <- unique(g[g!=t])
    # Go through remaining group labels and sum up over P_tt' then, see Chapter 2.1 in gMWT manuscript
      if(alg=="Cnaive")
      {
	for(i in 1:length(tprime))
	{
	    result <- result + sum(g==tprime[i])*getP.Cnaive(X[g==t],X[g==tprime[i]])
	}
      } else if(alg=="Rgrid"){
	for(i in 1:length(tprime))
	{
	    result <- result + sum(g==tprime[i])*getP.grid(X[g==t],X[g==tprime[i]])
	}
      } else if(alg=="Rnaive"){
	for(i in 1:length(tprime))
	{
	    result <- result + sum(g==tprime[i])*getP.Rnaive(X[g==t],X[g==tprime[i]])
	}
      } else {
	stop("For the single probabilistic indices are the only options: 'Cnaive', 'Rgrid' and 'Rnaive'\n")
      }
      result/(N-Nt)
  } # Ned of function PE1
#----------------------------------------------------------------------------------------------------------------------------------------------

# This function is needed for the P_tt' case .
# Daniel, Tampere 20-10-2012
# Function tests, 20-10.2012

PE2 <- function(X,g,comb,alg){
  if(alg=="Cnaive"){
    result <- getP.Cnaive(X[g==comb[1]],X[g==comb[2]])
  } else if(alg=="Rgrid"){
    result <- getP.grid(X[g==comb[1]],X[g==comb[2]])
  } else if(alg=="Rnaive"){
    result <- getP.Rnaive(X[g==comb[1]],X[g==comb[2]])
  } else {
     stop("For the probabilistic indices based on pairs are the only options: 'Cnaive', 'Rgrid' and 'Rnaive'\n")
  }
  result
}
#----------------------------------------------------------------------------------------------------------------------------------------------


# This function is needed for the P_tt't'' case .
# Daniel, Tampere 20-10-2012
# Function tests, 20-10.2012

PE3 <- function(X,g,comb,alg,nper=1){
  if(alg=="Cnaive"){
    result <- getP.Cnaive(X[g==comb[1]],X[g==comb[2]],X[g==comb[3]])
  } else if(alg=="Rsubmat"){
    result <- getP.Rsub(X[g==comb[1]],X[g==comb[2]],X[g==comb[3]])
  } else if(alg=="Rnaive"){
    result <- getP.Rnaive(X[g==comb[1]],X[g==comb[2]],X[g==comb[3]])
  } else if(alg=="Csubmat"){
    result <- getP.Csub(X[g==comb[1]],X[g==comb[2]],X[g==comb[3]],nper)
  } else {
     stop("For the probabilistic indices based on triples are the only options: 'Cnaive', 'Rsubmat', 'Rnaive' and 'Csubmat'\n")
  }
  result
}
#----------------------------------------------------------------------------------------------------------------------------------------------