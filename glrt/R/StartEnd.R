StartEnd <-
function(alpha)
{

n = nrow(alpha)

startend = matrix(0, n, 2)

cusum = apply(alpha, 1, cumsum)  

for(i in 1:n)

{

startend[i,1] = which(cusum[,i] == 1, arr.ind=T)[1]

startend[i,2] = which(cusum[,i] == max(cusum[,i]), arr.ind=T)[1]
}

startend

}

