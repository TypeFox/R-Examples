Chisqstat3 <-
function(U, V, counts)

{

k = length(U)

n = sum(counts)

N = apply(counts, 2, sum)


if(N[1] == 0) 
{
chisq = t(U[1:(k-1)]) %*% solve(V[1:(k-1), 1:(k-1)]) %*% U[1:(k-1)] / n
df = k - 1
}
else
{
equal = 1

ratios = t(apply(counts, 1, "/", N))
r = 2

while(equal == 1 && r <= k)

{

if(any(ratios[r,] != ratios[1,]))

equal = 0

r = r + 1

}

if(equal == 1)

{

chisq = t(U[1:(k-1)]) %*% solve(V[1:(k-1), 1:(k-1)]) %*% U[1:(k-1)] / n
df = k - 1

}

else

{

chisq = t(U) %*% solve(V) %*% U / n

df = k

}
}

c(chisq, df)

}

