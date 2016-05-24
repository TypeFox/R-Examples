#questa funzione, dato un vettore, trova il punto che può essere considerato come minimo locale. Rispetto alla definzione data da luck. Vale a dire un punto intorno massimo circondato da punti più alti (n.points, nella funzione).
# di default n.points=2

localmin <- function(x, n.points=2){
	vet=x
	dat=data.frame(index=1:length(vet))
	# creo un data.frame con gli indici. Mi serve per utilizzare la funzione apply.
	
	
	test=function(i,vet,n.points) {
	indices=c((i-n.points):(i-1), (i+1):(i+n.points)) # select observations n.points-far from the i (i)
	indices=indices[indices>0&indices<(length(vet))] #esclude 0 or negative indices and indices outside the vector length
	response=((!any(vet[indices]<=vet[i]))&length(indices)==(n.points*2)) # include only values in which i < surronding n.points points
																   # and there are all n.points*2 points
	return(response)
	}
	
	candidates.indices=apply(dat, 1, function(k){test(i=k, vet, n.points)})
	
	if(any(candidates.indices)){
	candidates=vet[candidates.indices]
	return(min(candidates)) } else {
		return(NA)
	}
	
	}
