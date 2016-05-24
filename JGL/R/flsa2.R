flsa2 <-
function(A,L,lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
{
# 1st apply fused penalty:
S1 = abs(A[[1]]-A[[2]])<=2*lam2/L
X1 = (A[[1]]+A[[2]])/2
Y1 = X1

S2 = (A[[1]] > A[[2]]+2*lam2/L)
X2 = A[[1]] - lam2/L
Y2 = A[[2]] + lam2/L

S3 = (A[[2]] > A[[1]]+2*lam2/L)
X3 = A[[1]] + lam2/L
Y3 = A[[2]] - lam2/L

X = soft(a = S1*X1 + S2*X2 + S3*X3, lam = lam1/L, penalize.diagonal=penalize.diagonal)
Y = soft(a = S1*Y1 + S2*Y2 + S3*Y3, lam = lam1/L, penalize.diagonal=penalize.diagonal)

return(list(X,Y))
}

