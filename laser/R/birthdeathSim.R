birthdeathSim <- function(b=.1, d=0, CladeSize, NumberMissing, NumberOfReps){
		v <- vector("list", NumberOfReps);
# Klaus 15.4.2013
#		for (i in 1:NumberOfReps) v[[i]] <- birthdeath.tree(b = 0.1, d = 0, taxa.stop = CladeSize);
#		v <- lapply(v, old2new.phylo);
#		v <- lapply(v, drop.tip, as.character(sample(1:(CladeSize),  NumberMissing)));	
                for (i in 1:NumberOfReps) v[[i]] <- sim.bdtree(b = 0.1, d = 0, stop = "taxa", n = CladeSize);
	        v <- lapply(v, drop.tip, sample(CladeSize,  NumberMissing));
		fx <- function(x)return(rev(sort(as.numeric(branching.times(x)))));
		v <- lapply(v, fx);
		v <- as.matrix(t(as.data.frame(v)));
		
		return(v);
}




