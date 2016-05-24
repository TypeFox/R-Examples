Dirichlet.multinomial <-
function(Nrs, shape){
if(missing(Nrs) || missing(shape))
stop("Nrs and/or shape missing.")

if(Nrs <= 0 || shape <= 0)
stop("Nrs and shape must be positive.")

for(n in Nrs){
if(any(n != n[1])){
warning("Unequal number of reads across samples.")
break
}
}

Sample.counts <- matrix(0, length(Nrs), length(shape))
for(i in 1:length(Nrs))
Sample.counts[i,] <- rmultinom(n=1, size=Nrs[i], prob=dirmult::rdirichlet(1, shape))

return(Sample.counts)
}
