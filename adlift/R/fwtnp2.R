fwtnp2 <-
function (input, f, nkeep = 2, intercept = TRUE, initboundhandl = 0, 
    neighbours = 1, closest = FALSE, LocalPred = 1, do.W = FALSE, 
    varonly = FALSE) 
{
n <- length(f)

v <- 0
W <- 0
Wnew<-vnew<-NULL
  
if((do.W==1)&(varonly==1)){
	varonly<-FALSE
}

ex<-do.W+varonly  

if (ex==1) {
	W<-matrix(0,n,n)
}

if (varonly) {
       	v <- rep(0, times = n)
}

if ((n - nkeep) > 0) {
	coeff <- rep(0, times = n)
        pointsin <- coeff
        lengthsremove <- rep(0, times = n - nkeep)
        lca <- matrix(0, n - nkeep, 6 * neighbours + 5)
        ans <- .C("fwtnp", as.double(input), as.double(f), as.integer(nkeep), 
            as.integer(intercept), as.integer(initboundhandl), 
            as.integer(neighbours), as.integer(closest), as.integer(LocalPred), 
            as.integer(n), coeff = as.double(coeff), lr = as.double(lengthsremove), 
            lengths = as.double(coeff), lca = as.double(t(lca)), 
            po = as.integer(pointsin), nc = as.integer(0), as.integer(do.W), 
            W = as.double(t(W)), as.integer(varonly), v = as.double(v), 
            PACKAGE = "adlift")

        ans$po <- ans$po[1:nkeep]
        ans$lengths <- ans$lengths[1:nkeep]
        ans$lca <- matrix(ans$lca[1:((n - nkeep) * ans$nc)], ncol = ans$nc, byrow = TRUE)
        coeff <- ans$coeff
        lengthsremove <- ans$lr
        lca <- ans$lca
        pointsin <- ans$po
        lengths <- ans$lengths

	if(ex==1){
	        if(do.W){
		        Wnew <- matrix(ans$W, n, n, byrow = TRUE)
			vnew<-NULL
		}
		else{        
	        	vnew <- ans$v
			Wnew<-NULL
		}
	}
}
else {
	coeff <- f
        pointsin <- order(input)
        lengths <- lengthsremove <- lca <- NULL
}

return(list(x=input,coeff = coeff, lengthsremove = lengthsremove, 
        lca = lca, pointsin = pointsin, lengths = lengths, W = Wnew, 
        v = vnew))

}

