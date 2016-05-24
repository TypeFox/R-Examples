`TotalVar` <-
function(x,DiscretePrevision,n,s){

   x.nodevalues <- DiscretePrevision[[1]]
   omega.nodevalues <- DiscretePrevision[[2]]
   UpperPrevisions <- DiscretePrevision[[3]]

   m <- ncol(omega.nodevalues)     # Number of additional supporting nodes

  # Loosing some more or less unneccessary nodes gives:
   ListRedNnodev <- FctReducedNodevalues(n,m,s,x.nodevalues,omega.nodevalues)
   x.nodevalues <- ListRedNnodev[[1]]
   x.frequency <- ListRedNnodev[[2]]
   omega.nodevalues <- ListRedNnodev[[3]]
   mr <- ListRedNnodev[[4]]

  # The number of observed nodes (which may differ from n now)
   nr <- length(x.frequency)    # 'nr' appreviates: n-reduced

  # Updating the number of additional supporting nodes
   m <- mr

  # all values of the functions f at all remaining nodes
   nodevalues <- cbind(x.nodevalues,omega.nodevalues)


  # Buiding the linear function 'a' which has to be optimized:
   a <- BuildOptVec(nr,m)

  # Building the matrix A:
   A0 <- BuildMatrix1(nodevalues,nr)
   A1234 <- BuildMatrix2(nr,m)
   A <- BuildMatrix(A0,A1234)

  # Building the bounds b:
   b0 <- BuildBounds1(UpperPrevisions)
   b1234 <- BuildBounds2(n,nr,x.frequency,m)
   b <- BuildBounds(b0,b1234)

  # Solving the maximization problem by use of LpSolve:
   library(lpSolve)

   l <- rep("<=",length(A[,1]))
   ergebnis <- lp(direction="max", objective.in=a, const.mat=A, const.dir=l, const.rhs=b)
   list(2*(1-ergebnis$objval), ergebnis$status)

}

