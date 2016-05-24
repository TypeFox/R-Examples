"bootstrap.oncotree" <-
function(otree, R, type=c("nonparametric","parametric")){
		  type <- match.arg(type)
		
	    parentlist <- function(x,i){
		    new.d <- x[i,]
		    new.tree <- oncotree.fit(new.d, error.fun=NULL)
		    result <- numeric(ncol(x))
		    names(result) <- colnames(x)
			result[match(new.tree$parent$child, names(result))] <- new.tree$parent$parent.num
			result
		  }		
	    
#	   	original.parent <- otree$parent
	    nmut <- otree$nmut

	    if (type=="nonparametric"){
		    test<-boot(otree$data, parentlist, R)    
		    result<-test$t
	    }
	    else {	#type=="parametric"
	    	if (is.null(otree$eps)) stop("Need false positive and negative rates")
				result <- matrix(NA, nrow=R, ncol=nmut)
				N <- nrow(otree$data)
	      distr <- distribution.oncotree(otree, with.probs=TRUE, with.errors=TRUE,
		                     edge.weights="estimated")
				for (i in 1:R){
					ran.idx <- sample(1:nrow(distr), size=N, prob=distr$Prob, replace=TRUE)
					new.d <- distr[ran.idx, 2:otree$nmut]
					new.tree <- oncotree.fit(new.d, error.fun=NULL)
					result[i,match(new.tree$parent$child, otree$parent$child)] <- new.tree$parent$parent.num
				}
			}
	
	#frequency of each possible child-parent edge
		parent.freq <- apply(result, 2, function(x)table(factor(x,levels=0:nmut)))
		dimnames(parent.freq) <- list(Parent=c("", otree$parent$child),
		                              Child=otree$parent$child)
		#mostcommonparent <- apply(result, 2, 
		#       function(x){tb <- table(x); as.numeric(names(tb)[which.max(tb)])})
		mostcommonparent <- apply(parent.freq, 2, which.max)-1
		
		#What are the frequencies of the other trees?
		tree.list <- apply(result, 1, function(x) paste(x, sep="", collapse="."))
		tree.freq1 <- as.data.frame(table(tree.list))
		Freq <- tree.freq1[,2]
		o <- order(-Freq)
		tree.freq <- tree.freq1[o,1:2]
		colnames(tree.freq) <- c("Tree","Freq")
		rownames(tree.freq) <- 1:nrow(tree.freq)
		bootout <- list(original=otree$parent, consensus=mostcommonparent, 
		              parent.freq=parent.freq,
		              tree.list=tree.freq, boot.type=type)
		class(bootout) <- "boottree"
		return(bootout)
}

