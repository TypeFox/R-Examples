S.m <- function(gam.fit){
  Ss <- list()
  off <- rank <- 0 
  jj <- 1
	for(j in 1:length(gam.fit$sp) ){
		Ss[[j]] <- gam.fit$smooth[[j]]$S[[1]]
            	rank[j] <- gam.fit$smooth[[j]]$rank
                off[j]  <- gam.fit$smooth[[j]]$first.para                                               
       }
list(rank=rank,off=off,Ss=Ss)
}


