"oncotree.fit" <-
function(dataset, error.fun=function(x,y){sum((x-y)^2)}){
    cur.matrix <- formatData(dataset)
    parent <- build.parent(cur.matrix)
 #   plotinfo <- build.plot(parent)
    
    estimate.error.rates <- function(otree, error.fun){
		mse <- function(eps, observed){
			err <- try({error.rates(otree) <- eps}, silent=TRUE)
			if (is.character(err)) return(otree$nmut)
			md.eps <- marginal.distr(otree, with.errors=TRUE)
			error.fun(observed[-1],md.eps[-1])
		  }
		
		obs <- colMeans(otree$data)
		min.p <- min(marginal.distr(otree, with.errors=FALSE))
		res <- optim(c(0,0), mse, method="L-BFGS-B", lower=c(0,0), upper=c(min.p-0.001, 1-min.p), observed=obs)
		res$par
	  }

    otree <- list(data=cur.matrix, nmut=ncol(cur.matrix), parent=parent)
		          # ,level=plotinfo$level, numchild=plotinfo$numchild, 
              #  levelnodes=plotinfo$levelnodes, levelgrp=plotinfo$levelgrp)
    if (!is.null(error.fun)){              
	    eps.hat <- estimate.error.rates(otree, error.fun=error.fun)
	    error.rates(otree) <- eps.hat
    }
    class(otree)<-"oncotree"
    return(otree)
}

