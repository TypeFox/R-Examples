loglikDM <-
function(data, alphap){
data <- data[,colSums(data)!=0]
alphap <- alphap[alphap!=0]

ll <- sum(lgamma(rowSums(data)+1) + lgamma(sum(alphap)) - lgamma(sum(alphap)+rowSums(data))) + 
sum(rowSums(lgamma(sweep(data, 2, alphap, "+")) - lgamma(data+1) - lgamma(t(replicate(nrow(data), alphap)))))

return(ll)
}
