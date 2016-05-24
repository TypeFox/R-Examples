################################################################################
## General
################################################################################

#' Scramble a vector
#' 
#' This function takes a vector and returns it in a random order
#' 
#' @param x a vector
#' @return a vector, now in random order
#' @export
#' @examples
#' scramble(1:10)
scramble <- function(x) {
  if(length(x) == 1) x
  else sample(x, length(x))
}

#' Converts scores to a ranking
#' 
#' takes in vector of scores (with the largest score being the one most preferred)
#' and returns back a vector of WINNER, SECOND PLACE, ... LAST PLACE
#' 
#' @param scores the scores (e.g. means) of a set of alternatives
#' @return an ordering of the index of the winner, second place, etc.
#' @export
#' @examples
#' scores <- Generate.RUM.Parameters(10, "exponential")$Mean
#' scores.to.order(scores)
scores.to.order <- function(scores) (1:length(scores))[order(-scores)]

#' Helper function for the graphing interface
#'
#' As named, this function takes a vector where each element is a mean, 
#' then returns back a list, with each list item having the mean
#'
#' @param Parameters a vector of parameters
#' @param name Name of the parameter
#' @return a list, where each element represents an alternative and has a Mean value
#' @export
convert.vector.to.list <- function(Parameters, name = "Mean") {
  m <- length(Parameters)
  List <- rep(list(NA), m)
  for(i in 1:m) {
  	List[[i]] <- list(placeholdername = Parameters[i])
  	names(List[[i]])[1] <- name
  }
  List
}

################################################################################
## Breaking
################################################################################

#' Generate a matrix of pairwise wins
#' 
#' This function takes in data that has been broken up into pair format.
#' The user is given a matrix C, where element C[i, j] represents 
#' (if normalized is FALSE) exactly how many times alternative i has beaten alternative j
#' (if normalized is TRUE)  the observed probability that alternative i beats j
#' 
#' @param Data.pairs the data broken up into pairs
#' @param m the tot
#' al number of alternatives
#' @param weighted whether or not this generateC should use the third column of Data.pairs as the weights
#' @param prior the initial "fake data" that you want to include in C. A prior 
#' of 1 would mean that you initially "observe" that all alternatives beat all
#' other alternatives exactly once.
#' @param normalized if TRUE, then normalizes entries to probabilities
#' @return a Count matrix of how many times alternative i has beat alternative j
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' generateC(Data.Test.pairs, 5)
generateC <- function(Data.pairs, m, weighted = FALSE, prior = 0, normalized = TRUE) {
  
  # Converts a pairwise count matrix into a probability matrix
  # 
  # @param C original pairwise count matrix
  # @return a pairwise probability matrix
  # @export
  # @examples
  # C= matrix(c(2,4,3,5),2,2)
  # normalizeC(C)
  normalizeC <- function(C){
    normalized <- C / (C + t(C))
    diag(normalized) <- 0
    normalized
  }
     
  
  # C is the transition matrix, where C[i, j] denotes the number of times that
  C <- matrix(data = prior, nrow = m, ncol = m) - prior * m * diag(m)
  
  pairs.df <- data.frame(Data.pairs)
  
  if(ncol(pairs.df) > 2) {
    names(pairs.df) <- c("i", "j", "r")
    C.wide <- ddply(pairs.df, .(i, j), summarize, n = length(i), r = sum(r))
  }
  else {
    names(pairs.df) <- c("i", "j")
    C.wide <- ddply(pairs.df, .(i, j), summarize, n = length(i))
  }
  
  for(l in 1:nrow(C.wide)) {
    # i wins, j loses
    i <- C.wide[l, "i"]
    j <- C.wide[l, "j"]
    if(ncol(C.wide > 2) & weighted) 
      C[i, j] <- C[i, j] + C.wide[l, "r"]
    else
      C[i, j] <- C[i, j] + C.wide[l, "n"]
  }
  
  if(normalized) normalizeC(C)
  else C
}

#' Turns inference object into modeled C matrix.
#' 
#' For parametric models, plug in a pairwise function for get.pairwise.prob.
#' For nonparametric models, set nonparametric = TRUE
#'
#' @param Estimate inference object with a Parameter element, with a list of parameters
#' for each alternative
#' @param get.pairwise.prob (use this if its a parametric model) function that takes in two lists of parameters and 
#' computes the probability that the first is ranked higher than the second
#' @param nonparametric set this flag to TRUE if this is a non-parametric model
#' @param ... additional arguments passed to generateC.model.Nonparametric (bandwidth)
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' Estimate <- Estimation.Normal.GMM(Data.Test.pairs, 5)
#' generateC.model(Estimate, Normal.Pairwise.Prob)
generateC.model <- function(Estimate, get.pairwise.prob = NA, nonparametric = FALSE, ...) {
  if(nonparametric)
    generateC.model.Nonparametric(Estimate, ...)
  else {
    m <- length(Estimate$Parameters)
    C.matrix.model <- matrix(0, nrow = m, ncol = m)
    for(i in 1:m) for(j in 1:m) if(i != j) C.matrix.model[i, j] <- get.pairwise.prob(Estimate$Parameters[[i]], Estimate$Parameters[[j]])
    C.matrix.model
  }
}

#' Breaks full or partial orderings into pairwise comparisons
#' 
#' Given full or partial orderings, this function will generate pairwise comparison
#' Options
#' 1. full - All available pairwise comparisons. This is used for partial
#' rank data where the ranked objects are a random subset of all objects
#' 2. adjacent - Only adjacent pairwise breakings
#' 3. top - also takes in k, will break within top k
#' and will also generate pairwise comparisons comparing the
#' top k with the rest of the data
#' 4. top.partial - This is used for partial rank data where the ranked 
#' alternatives are preferred over the non-ranked alternatives
#' 
#' The first column is the preferred alternative, and the second column is the
#' less preferred alternative. The third column gives the rank distance between
#' the two alternatives, and the fourth column enumerates the agent that the
#' pairwise comparison came from.
#' 
#' @param Data data in either full or partial ranking format
#' @param method - can be full, adjacent, top or top.partial
#' @param k This applies to the top method, choose which top k to focus on
#' @return Pairwise breakings, where the three columns are winner, loser and rank distance (latter used for Zemel)
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
Breaking <- function(Data, method = "full", k = NULL) {
  m <- ncol(Data)
  
  pair.full <- function(rankings) pair.top.k(rankings, length(rankings))
  
  pair.top.k <- function(rankings, k) {
    pair.top.helper <- function(first, rest) {
      pc <- c(); 
      z <- length(rest)
      for(i in 1:z) pc <- rbind(pc, array(as.numeric(c(first, rest[i], i))))
      pc
    }
    if(length(rankings) <= 1 | k <= 0) c()
    else rbind(pair.top.helper(rankings[1], rankings[-1]), pair.top.k(rankings[-1], k - 1))
  }
  
  pair.adj <- function(rankings) {
    if(length(rankings) <= 1) c()
    else rbind(c(rankings[1], rankings[2]), pair.adj(rankings[-1]))
  }
  
  pair.top.partial <- function(rankings, m) {
    # this is used in the case when we have missing ranks that we can
    # fill in at the end of the ranking. We can assume here that all 
    # ranked items have higher preferences than non-ranked items
    # (e.g. election data)
    
    # the number of alternatives that are not missing
    k <- length(rankings)
    
    # these are the missing rankings
    missing <- Filter(function(x) !(x %in% rankings), 1:m)
    
    # if there is more than one item missing, scramble the rest and place them in the ranking
    if(m - k > 1) missing <- scramble(missing)
    
    # now just apply the top k breaking, with the missing elements
    # inserted at the end
    pair.top.k(c(rankings, missing), k)
  }
  
  break.into <- function(Data, breakfunction, ...) {
    n <- nrow(Data)
    # applying a Filter(identity..., ) to each row removes all of the missing data
    # this is used in the case that only a partial ordering is provided
    tmp <- do.call(rbind, lapply(1:n, function(row) cbind(breakfunction(Filter(identity, Data[row, ]), ...), row)))
    colnames(tmp) <- c("V1", "V2", "distance", "agent")
    tmp
  }
  
  if(method == "full") Data.pairs <- break.into(Data, pair.full)
  if(method == "adjacent") Data.pairs <- break.into(Data, pair.adj)
  if(method == "top") Data.pairs <- break.into(Data, pair.top.k, k = k)
  if(method == "top.partial") Data.pairs <- break.into(Data, pair.top.partial, m = m)
  
  Data.pairs
}

##################################################Sampling#############################################

# Conditional truncated sampler for normal
# 
# @param X sample values for alternative, m dimensional
# @param pi ranks pi[1] is the best pi[m] is the worst sampling X[i]
# @param theta m dimentional parameter values
samplen=function(dist, parameter, i, X, pi, race=FALSE)
{
  m <- length(X)
  
  #case for full orders
  if(!sum(pi==0)){
    rank <- which(pi==i)

    if(rank>1 & rank<m ){
      lower <- X[as.numeric(pi[rank+1])]
      upper <- X[as.numeric(pi[rank-1])]
    }
    
    if(rank==m){
      lower <- -Inf
      upper <- X[as.numeric(pi[rank-1])]
    }
    
    if(rank==1){
      lower <- X[as.numeric(pi[rank+1])]
      upper <- Inf
    }
    
  }
  
  #case for nonfull orders
  if(sum(pi==0)){
    listemp <- as.numeric(pi[which(pi>0)])
    listempz <- setdiff(1:m,listemp)
    
    if(!sum(pi==i)){
      if(pi[m]==0){
        lower <- -Inf
        upper <- min(X[listemp])
        if(race){
          upper <- Inf
        }
      }
      if(pi[1]==0){
        upper <- Inf
        lower <- max(X[listemp])
        if(race){
          lower <- -Inf
        }
      }
    }
    
    if(sum(pi==i)){
      rank <- which(pi==i)
    
      if(rank>1 & rank<m ){
        if(pi[rank+1]>0 & pi[rank-1]>0){  
          lower <- X[as.numeric(pi[rank+1])]
          upper <- X[as.numeric(pi[rank-1])]
        }
        if(pi[rank+1]==0 & pi[rank-1]>0){  
          lower <- max(X[listempz])
          upper <- X[as.numeric(pi[rank-1])]
          if(race){
            lower <- -Inf
          }
        }
        if(pi[rank+1]>0 & pi[rank-1]==0)
        {  
          
          lower <- X[as.numeric(pi[rank+1])]
          upper <- min(X[listempz])
          if(race){
            upper <- Inf
          }
        }
        
      }
      
      if(rank==m){
        if(pi[rank-1]>0){
          lower <- -Inf
          upper <- X[as.numeric(pi[rank-1])]
        }
        if(pi[rank-1]==0){
          lower <- -Inf
          upper <- min(X[listempz])
          if(race){
            upper <- Inf
          }
        }
        
      }
      
      if(rank==1){
        if(pi[rank+1]>0){
          lower <- X[as.numeric(pi[rank+1])]
          upper <- Inf
        }
        if(pi[rank+1]==0){
          lower <- max(X[listempz])
          upper=Inf
          if(race){
            lower <- -Inf
          }
        }
      }
    }
    
  }
  ##############
  if(upper < lower){
    upper <- abs(lower)+.01
  }
  
  
  if(dist=="norm" || dist=="norm.fixedvariance"){
    OldSample <- X[i]
    X[i] <- rtrunc(1, spec=dist, a = lower, b = upper,mean=parameter$Mu[i],sd=parameter$Var[i]^.5)
    
    if( X[i]==Inf | X[i]==-Inf | is.na(X[i]) ){
      X[i] <- OldSample
    }
    
    if( (X[i]<lower) | (X[i]>upper) ){
      X[i] <- OldSample
    }
  }
  
  if(dist=="exp"){ 
    OldSample <- X[i]
    X[i] <- rtrunc(1, spec=dist, a = lower, b = upper,rate=parameter$Mu[i]^-1)
    
    if(X[i]==Inf | X[i]==-Inf | is.na(X[i])){
      X[i] <- OldSample
    }
    
    if( (X[i]<lower) | (X[i]>upper) ){
      X[i] <- OldSample
    }
    
  }
  
  X 
}

# Gibbs sampler
#
# Samples S samples Using a Gibbs sampler from the joint distribution of normal or exponential with a constraint on the
# order of samples given by pi as the rank. 
#   
# @param S number of samples
# @param X sample values for alternative, m dimensional
# @param pi ranks pi[1] is the best pi[m] is the worst sampling X[i]
# @param dist distribution of utilities
# @param parameter m dimensional parameter values
GibbsSampler=function(S, X, pi, dist, parameter, race=FALSE)
{
  m <- length(X)
  out <- matrix(0,S,m)
  out[1,] <- X  
  
  #gcon=mean(-log(rexp(100000,1)))
  
  for(s in 2:S){
    i <- sample.int(m, size = 1, replace = FALSE, prob = NULL)
    out[s,] <- samplen(dist,parameter,i,out[s-1,],pi,race)
    #print(i)
  }
  
  list(M1=colMeans(out[(round(S/10)):S,]),M2=colMeans(out[(round(S/10)):S,]^2), X=out[S,])
}

#################################Helpers for GMM for Normal
#

f= function(Mu, SD=0, Var = FALSE)
{
  if(Var==FALSE)
  {
    m=length(Mu)
    A=matrix(0,m,m)
    for(i in 1:m)
    {
      for(j in 1:m)
      {
        A[i,j]=pnorm(Mu[i]-Mu[j],0,sqrt(2))
      }
    }
    diag(A)=1-colSums(A);
  }
  
  if(Var==TRUE)
  {
    m=length(Mu)
    A=matrix(0,m,m)
    for(i in 1:m)
    {
      for(j in 1:m)
      {
        A[i,j]=pnorm(Mu[i]-Mu[j],0,sd=sqrt(SD[i]^2+SD[j]^2))
      }
    }
    diag(A)=1-colSums(A);
  }
  diag(A)=0
  ;A
}

VarMatrix = function(SD)
{
  m=length(SD)
  Var=matrix(0,m,m)
  for(i in 1:m)
  {
    for(j in 1:m)
    {
      Var[i,j]=(SD[i]^2+SD[j]^2)
    }
  }
  ;Var
}

delta= function(Mu, SD=0, Var = FALSE)
{
  if(Var==FALSE)
  {
    m=length(Mu)
    A=matrix(0,m,m)
    for(i in 1:m)
    {
      for(j in 1:m)
      {
        A[i,j]=Mu[i]-Mu[j]
      }
    }
  }
  if(Var==TRUE)
  {
    m=length(Mu)
    A=matrix(0,m,m)
    for(i in 1:m)
    {
      for(j in 1:m)
      {
        A[i,j]=2^.5*(Mu[i]-Mu[j])/(SD[i]^2+SD[j]^2)^.5
      }
    }
  }
  ;A
}

#full break
fullbreak= function(data)
{
  n=dim(data)[1]
  m=dim(data)[2]
  A=matrix(0,m,m)
  for(i in 1:n)
  {
    for(j in 1:(m-1))
    {
      for(k in (j+1):m)
      {
        A[data[i,j],data[i,k]]=A[data[i,j],data[i,k]]+1
      }
    }
  }
  A=A/n;
  diag(A)=.5-colSums(A);
  ;A
}

#top-c break
topCbreak= function(data,c)
{
  n=dim(data)[1]
  m=dim(data)[2]
  A=matrix(0,m,m)
  for(i in 1:n)
  {
    for(j in 1:(c-1))
    {
      for(k in (j+1):c)
      {
        A[data[i,j],data[i,k]]=A[data[i,j],data[i,k]]+1;
      }
    }
  }
  
  An=matrix(0,m,m)
  for(i in 1:m)
  {
    for(j in 1:m)
    {
      if((A[i,j]+A[j,i])!=0)
      {
        An[i,j]=A[i,j]/(A[i,j]+A[j,i])
      }
    }
  }
  diag(An)=.5-colSums(An);
  ;An
}


Analyze=function(brokendata, m)
{
  e=eigen(brokendata)
  ;sort(e$vectors[,m],index=TRUE,decreasing=TRUE)$x
}

###

#' Converts a matrix into a table
#' 
#' takes a matrix and returns a data frame with the columns being row, column, entry
#' 
#' @param A matrix to be converted
#' @param uppertriangle if true, then will only convert the upper right triangle of matrix
#' @return a table with the entries being the row, column, and matrix entry
#' @export
turn_matrix_into_table <- function(A, uppertriangle = FALSE) {
  m <- dim(A)[1]
  transcribed.table <- data.frame()
  if(uppertriangle) for(i in 1:(m-1)) for(j in (i+1):m) transcribed.table <- rbind(transcribed.table, c(i, j, A[i, j]))
  else for(i in 1:m) for(j in (1:m)) if(i != j) transcribed.table <- rbind(transcribed.table, c(i, j, A[i, j]))
  names(transcribed.table) <- c("row", "column", "prob")
  transcribed.table
}

################################################################################
# Calculating Pairwise Probabilities
################################################################################


#' Pairwise Probability for PL Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
PL.Pairwise.Prob <- function(a, b) {
  # this is flipped because in the PL case, smaller utilities are ranked higher than
  # lower utilities
  1 - a$Mean / (a$Mean + b$Mean)
}

#' Pairwise Probability for Zemel
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Zemel.Pairwise.Prob <- function(a, b){
  exp(a$Score - b$Score) / (exp(a$Score - b$Score) + exp(b$Score - a$Score))
}

#' Pairwise Probability for Normal Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Normal.Pairwise.Prob <- function(a, b) {
  # Let W = X - Y
  if(is.null(a[["Variance"]])) a$Variance <- 1
  if(is.null(b[["Variance"]])) b$Variance <- 1
  
  mu <- a$Mean - b$Mean
  sigma <- sqrt(a$Variance + b$Variance)
  #P(X - Y > 0)
  1 - pnorm(-mu/sigma)
}

#' Pairwise Probability for Normal Multitype Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Normal.MultiType.Pairwise.Prob <- function(a, b) {
  mean.grid <- expand.grid(a$Mean, b$Mean)
  variance.grid <- expand.grid(a$Variance, b$Variance)
  gamma.grid <- expand.grid(a$Gamma, b$Gamma)
  n <- nrow(mean.grid)
  
  probs <- rep(NA, n)
  for(i in 1:n) {
    probs[i] <- Normal.Pairwise.Prob(list(Mean = mean.grid[i, 1], Variance = variance.grid[i, 1]), list(Mean = mean.grid[i, 2], Variance = variance.grid[i, 2]))
    probs[i] <- probs[i] * prod(gamma.grid[i, ])
  }
  sum(probs)
}

#' Pairwise Probability for PL Multitype Model
#' 
#' Given alternatives a and b (both items from the inference object)
#' what is the probability that a beats b?
#' 
#' @param a list containing parameters for a
#' @param b list containing parameters for b
#' @return probability that a beats b
#' @export
Expo.MultiType.Pairwise.Prob <- function(a, b) {
  mean.grid <- expand.grid(a$Mean, b$Mean)
  n <- nrow(mean.grid)
  
  probs <- rep(NA, n)
  for(i in 1:n) {
    probs[i] <- PL.Pairwise.Prob(list(Mean = mean.grid[i, 1]), list(Mean = mean.grid[i, 2]))
    probs[i] <- probs[i] * prod(mean.grid[i, ])
  }
  sum(probs)
}

dgumbel <- function (x, scale = 1, location = 0, log = FALSE) 
{
  fx <- 1/scale * exp(-(x - location)/scale - exp(-(x - location)/scale))
  if (log) 
    return(log(fx))
  else return(fx)
}