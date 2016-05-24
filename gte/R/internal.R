# This function was extracted from the function 'Aintmap' in the package 'interval' 
# version 2011-04-07 : R/icfit.R

get.intmap <- function(L, R){
  n<-length(L)
  
  # Lin and Rin indicates which of '(' or '[' is going to be used in the output 
  Lin<-rep(FALSE,n)
  Rin<-rep(TRUE,n)
  Lin[L==R]<-TRUE     ## left because it was in the original code, but these 2 lines are useless here
  Rin[R==Inf]<-FALSE  ## since we got rid of L==R and of Inf values before calling get.intmap.
  
  # calculate a small number, eps, to differentiate between e.g.,  [L,R] and (L,R]
  # we will treat, (L,R] as [L+eps,R], and [L,R) as [L,R-eps] 
  # since eps is only 
  # used in ranking we do not need to make it super small
  # just smaller than the smallest difference
  LRvalues<-sort(unique(c(0,L,R,Inf)))
  eps<- min(diff(LRvalues))/2
  Le<-L
  Re<-R
  Le[!Lin]<-L[!Lin]+eps
  Re[!Rin]<-R[!Rin]-eps
  # let s be the vector of ordered L and R values with 
  # R values later when there are ties
  # then intmap are values s[i] and s[i+1] where s[i] is 
  # associated with L and s[i+1] is associated with R
  oLR<-order(c(Le,Re+eps/2) )
  # find the Turnbull intervals, or innermost intervals
  # this is the same as the primary reduction of 
  ### Aragon and Eberly (1992) J of Computational and Graphical
  ###     Statistics 1:129-140
  # label L=1 and R=2
  Leq1.Req2<-c(rep(1,n),rep(2,n))
  # order and see if an R is followed by an L
  # take difference of Leq1.Req2 after putting them in 
  # order, then if the difference is 1 then the R=2 is followed by L=1 
  flag<- c(0,diff( Leq1.Req2[oLR] ))
  R.right.of.L<- (1:(2*n))[flag==1]
  intmapR<- c(L,R)[oLR][R.right.of.L]
  intmapL<- c(L,R)[oLR][R.right.of.L - 1]
  intmapRin<- c(Lin,Rin)[oLR][R.right.of.L]
  intmapLin<- c(Lin,Rin)[oLR][R.right.of.L - 1]
  intmap <- matrix(c(intmapL,intmapR),byrow=TRUE,nrow=2)
  attr(intmap,"LRin") <- matrix(c(intmapLin,intmapRin),byrow=TRUE,nrow=2)
  
  intmap
}
