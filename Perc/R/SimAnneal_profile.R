
SimAnneal <- function(PMat, kmax=1000) {

#PMat = ret

## Part I: set up calculation

  N <- nrow(PMat)   # number of rows
  ColSum <- sapply(1:N, function(n) sum(PMat[,n])) # ColSum <- colSums(PMat) sum of prob of loses.
  
  # order by sum dom Prob of wins. 
  Order <- order(ColSum)  # order by sum prob of loses in incrasing order. That is order by sum prob of wins in decreasing order.
  # idea behind: When Matrix being re-organized by sum Prob wins/loses, this resulted in upper triangle with most domProb value > 0.5.
  C <- Cost(Reorder(PMat, Order)) # cost of the original domProb
  C1 <- C                      # assign the cost to C1 and Cb
  Cb <- C
  Ordb <- Order                # assign the order to Ordb
  profCn = numeric(kmax * 10)  # create a numeric vector
  profC = numeric(kmax * 10)
  profCn[1] = C                # assign cost to profCn[1] and profC[1]
  profC[1] = C

## Part II
  k <- 0                       # looping index used in while(k < kmax) statement.
  l <- 1
  while(k < kmax) {   # repeat this for kmax * 10 times.
    # idea: each time a new cost found, update C and ordb.
    # profC keep track of all C.
    # profCn keep track of all Cn found.
    # Cb: a single value representing the best cost.
    
    NewOrder <- FindNeighbor(PMat,Order) # determine a new order using find neighbor (by domProb)
    Cn <- Cost(Reorder(PMat, NewOrder))  # find the cost of the new order (Cn).
    if(Cn < C) {                         # if a smaller cost found, adopt the smaller cost and the new order
      C <- Cn
      Order <- NewOrder
      if(Cn < Cb) {                     # update Cb each time a better Cn found.
        Cb <- Cn
        Ordb <- NewOrder
      }
    } else
    if( AP(C, Cn, Temp(k/kmax)) > runif(1) ) {     # until order with smallest cost found,
      C <- Cn
      Order <- NewOrder
    }
    l <- l+1
    profCn[l] = Cn
    profC[l] = C
    if((l %% 10) == 0) k <- k+1
  }
  return(list(C1=C1, Cb=Cb, Ordb=Ordb, profC = profC, profCn = profCn))
}


Cost <- function(Mat) {
  # cost function: search lower triangle for Mat[i, j] > 0.5
  # if Mat[i, j] in lower triangle is greater than 0.5, 
  # then add cost log(2*(1-Mat[i, j])) * exp((N + 1 - j)*(i - j)*2/N^2).
  N <- nrow(Mat)
  Cost <- 0
  for(i in 2:N) for(j in 1:(i-1)) Cost <- Cost + max(0, -log(2*(1-Mat[i,j]))) * 
     exp((N+1-j)*(i-j)*2/N^2)
  return(Cost)
}

# This is way faster than the old Reorder function!
Reorder = function(Mat, Ord){
  Mat[Ord, Ord]
}
 

FindNeighbor <- function(Mat, Ord) {
  # determine a new order.
  # It transforms the lower triangle of the domProb matrix into a vector 
  #    and find the 80% quantile of this vector (Q). 
  # For each domProb value in the lower triangle, 
  #    if this value is greater than Q, save the row index (i) and the column index (j). 
  # All the row and column indexes for cells in lower triangle which have domProb value 
  #    greater than Q are composed into two vectors respectively, I and J. 
  # All these domProb values which are greater than Q is saved into a vector Pr. 
  # Then, sample one index (S) which indexing a combination I, J, 
  #     and Pr according to the value of Pr 
  #         (probability for each index of being sampled = Pr/sum(Pr)). 
  # Then switch the position of I[S] and J[S] in the order (ord) which is the new order.
  #Mat = ret
  #Ord = Order
  CMat <- Reorder(Mat, Ord)
  N <- nrow(Mat)
  Vals <- vector("numeric")
  for(i in 1:(N-1))  Vals <- c(Vals, CMat[(i+1):N,i]) # build a vector of lower triangle of domProb
  Q <- quantile(Vals, 0.8) # find the 80% quantile of the vector
  I <- vector("numeric")
  J <- vector("numeric")
  Pr <- vector("numeric")
  for(i in 2:N) for(j in 1:(i-1)) if(CMat[i,j] >= Q) { 
    # find column and row indexes where their corresponding domProb is greater than Q (the 80% quantile of the vector)
    I <- c(I, i)
    J <- c(J, j)
    Pr <- c(Pr, CMat[i,j]) # a vector of domProb.
  }
  Pr <- Pr/sum(Pr) # Proportion of domProb 
  S <- sample(x=1:length(Pr), size=1, prob=Pr) # sample one of the index
  Seq <- 1:N
  Seq[c(I[S],J[S])] <- Seq[c(J[S],I[S])]
  Ord <- Ord[Seq]
  return(Ord)
}


  

AP <- function(e1, e2, T) {
  if(e2 < e1) return(1)
  return( exp((e1 - e2)/T/0.1) )
}

Temp <- function(num) {
  return( 1 - exp((num-1)*5) )
}
