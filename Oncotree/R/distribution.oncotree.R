"distribution.oncotree" <-
function(otree, with.probs=TRUE, with.errors=FALSE,
          edge.weights=if (with.errors) "estimated" else "observed"){
	edge.weights <- match.arg(edge.weights, c("observed","estimated"))
	
	get.to.N.events <- function(atree, N, with.probs=TRUE,  edge.weights=c("observed","expected")){
		nevents <- atree$nmut-1
		res <- data.frame()
		if (N==0){
			root.outcome <- c(1,rep(0,nevents)) #the root has always occurred
			if (with.probs) p.root.only <- p.no.more(1, root.outcome, atree, edge.weights)
			res <- rbind(res, rbind(c(root.outcome, 0, if (with.probs) p.root.only)) )
			colnames(res) <- c(otree$parent$child, "NEvents", if (with.probs) "Prob")
		}
		else {
			res <- get.to.N.events(atree, N-1, with.probs=with.probs, edge.weights=edge.weights)
			for (i in which(res$NEvents==N-1)){
				outcome <- res[i,1:atree$nmut]
				dsc <- descendants(outcome, atree, with.probs, edge.weights)
				if (nrow(dsc[[1]])>0)
				  if (with.probs){
					  res <- rbind(res, data.frame(dsc[[1]], NEvents=N, Prob=dsc[[2]]*res$Prob[i],
					                               check.names=FALSE))
				  }
				  else {
					  res <- rbind(res, data.frame(dsc[[1]], NEvents=N, check.names=FALSE))
				  }
			}
		}
		res
	}
	
	descendants <- function(outcome, atree, with.probs, edge.weights=c("observed","expected")){
		desc <- matrix(0, nrow=0, ncol=atree$nmut)
		prob.adj <- numeric()
		level.empty <- tapply(as.numeric(outcome), atree$level, function(x)all(x==0))
		level.empty <- c(level.empty, TRUE, TRUE)  #add two artificial levels
		for (j in which(outcome==1)){  #breadth first enlargement
			lev.j <- atree$level[j]
			if (level.empty[lev.j+2]){
				j.children <- which(atree$parent$parent.num==j)
				max1 <- max(0,which((outcome==1)&(atree$level==lev.j+1))) #left-to-right within level
				j.children <- j.children[j.children>max1]
				for (ch in j.children){
					new.outcome <- outcome
					new.outcome[ch] <- 1
					desc <- rbind(desc, new.outcome)
					if (with.probs){
						inprob <- ifelse(edge.weights=="observed", 
						     atree$parent$obs.weight[ch], atree$parent$est.weight[ch])
						p.stop <- p.no.more(ch, new.outcome, atree, edge.weights)
						prob.adj <- c(prob.adj, p.stop*inprob/(1-inprob))
					}
				}
			}
		}
		list(desc, prob.adj)
	}
		
	p.no.more <- function(event, outcome, atree, edge.weights=c("observed","expected")){
		children <- which(atree$parent$parent.num==event)
		if (length(children)>0) {
			children <- children[outcome[children] == 0]
			if (edge.weights=="observed"){
				p <- prod(1-atree$parent$obs.weight[children])
			}
			else {
				p <- prod(1-atree$parent$est.weight[children])
			}
		}
		else {p <- 1}
		p			
	}
	
	plotinfo <- build.plot(otree$parent) #need level info
	otree$level <- plotinfo$level
	distr <- get.to.N.events(otree, N=otree$nmut-1, with.probs=with.probs, edge.weights=edge.weights)
	distr$NEvents <- NULL
	
	if (with.errors){
	    if (is.null(otree$eps)) stop("Need false positive and negative rates")
     	epos <- otree$eps["epos"]
     	eneg <- otree$eps["eneg"]

		val.list <- list()
		for (vx in otree$parent$child[-1]){
			val.list.add <- list(c(0,1))
			names(val.list.add) <- vx
			val.list <- c(val.list, val.list.add)
		}
		res <- do.call(expand.grid, val.list)
		attr(res,"out.attrs") <- NULL
		attr(res,"colnames") <- NULL
		trans.mat <- matrix(0, nrow=nrow(distr), ncol=nrow(res))
		
		for (i in 1:nrow(distr)){
			for (j in 1:nrow(res)){
				x <- as.numeric(distr[i,2:otree$nmut])
				y <- as.numeric(res[j,])
				trans.mat[i,j] <- (1-epos)^((1-x)%*%(1-y)) * epos^((1-x)%*%y) *
				     (1-eneg)^(x%*%y) * eneg^(x%*%(1-y))
			}
		}
		
		res$Prob <- as.numeric(distr$Prob %*% trans.mat)
		distr <- cbind(Root=1, res)
	}
	distr
}

