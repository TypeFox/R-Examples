Multinomial <-
function(Nrs, probs){
if(missing(Nrs) || missing(probs))
stop("Nrs and/or probs missing.")

for(n in Nrs){
if(any(n!=n[1])){
warning("Unequal number of reads across samples.")
break
}
}

Sample.counts <- matrix(0, length(Nrs), length(probs))
for(i in 1:length(Nrs))
Sample.counts[i,] <- rmultinom(n=1, size=Nrs[i], prob=probs)
data <- Sample.counts[, colSums(Sample.counts) != 0]

return(data)
}
