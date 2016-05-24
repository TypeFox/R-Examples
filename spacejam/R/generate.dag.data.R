generate.dag.data <-function(g,n,basesd = 1, basemean = 0, bfuns = function(x){cbind(x,x^2,x^3)}, funclist = NULL, usenorm=T){
	p <- vcount(g)
	X <- matrix(nrow=n,ncol=p)
	ord <- topological.sort(g)
	makeFunctions <- is.null(funclist)
	if(makeFunctions){
		funclist <- replicate(p,replicate(p,function(x){0},simplify=F),simplify=F)
	}
	for(i in ord){
		parents <- neighbors(g,i,mode="in")
		v <- length(parents)
		if(usenorm){
			X[,i] <-  rnorm(n=n,1)
		} else {
			X[,i] <- runif(n,min=-sqrt(12),max=sqrt(12))
		}
		if(v > 0){
			for(j in 1:v){
					if(makeFunctions){
						x <- X[,parents[j] ]
						xx <- as.matrix(bfuns(x))
						beta <- rnorm(ncol(xx),mean=basemean, sd= basesd)
						#beta <- beta/sd(xx%*%beta)
						funclist[[i]][[parents[j]]] <- function(x)(bfuns(x)%*%beta)	
						environment(funclist[[i]][[parents[j]]]) <- new.env()
						assign("beta", beta, envir = environment(funclist[[i]][[parents[j]]]))
						X[,i] <- X[,i] + scale(xx%*%beta)
					} else {
						X[,i] <- X[,i] + scale(funclist[[i]][[parents[j]]](X[,parents[j] ]))
					}
				}
			}
		X[,i] <- scale(X[,i])
		}
	return(list("X" = X, "funclist" = funclist))
}
