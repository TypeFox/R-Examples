Sub.r <-
function(z, t, Mat, i, j){
M <- sum(z)
D <- M - mean(sapply(1:t, function(j) Sub.S(z, t, Mat, j)))
Chat <- 1 - mean(sapply(1:t, function(j) Sub.S(z, t, Mat, j)) / sapply(1:t, function(j) Sub.n(z, t, Mat, j)))
Bij <- Sub.B(z, t, Mat, i, j)
ni <- Sub.n(z, t, Mat, i)
nj <- Sub.n(z, t, Mat, j)


MatA <- matrix(0, t, t)
MatB <- matrix(0, t, t)
MatC <- matrix(0, t, t)

for(r in 1:(t - 1)){
nr <- Sub.n(z, t, Mat, r)
for(s in (r + 1):t){
ns <- Sub.n(z, t, Mat, s)
MatA[r, s] = Sub.A(z, t, Mat, r, s)
MatB[r, s] = Sub.B(z, t, Mat, r, s)
MatC[r, s] = MatA[r, s] * (D / Chat * MatB[r, s] /(nr * ns) - 1)
}
}

rij <- Bij / (ni * nj) * (D / Chat + 1 / (t * Chat) * sum(MatC)) - 1
names(rij) <- paste("r", i, j, sep="")
return(rij)
}
