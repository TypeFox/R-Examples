"generate.data" <-
function(N, otree, with.errors=TRUE,
         edge.weights=if (with.errors) "estimated" else "observed",
         method=c("S","D1","D2")){
  
  method <- match.arg(method)
  edge.weights <- match.arg(edge.weights, c("observed","estimated"))
  if (!with.errors && method=="D2") method <- "D1"  #D2 only makes sense with errors
  
  if (method=="D1"){
    #calculate full joint distribution and sample from it
    distr <- distribution.oncotree(otree, with.probs=TRUE, with.errors=with.errors,
                                   edge.weights=edge.weights)
    ran.idx <- sample(1:nrow(distr), size=N, prob=distr$Prob, replace=TRUE)
    ran.data <- distr[ran.idx, 2:otree$nmut]
    rownames(ran.data) <- 1:N    
  }
  
  if (method=="D2"){
    #calculate joint distribution without errors, then "corrupt" the sample
    distr <- distribution.oncotree(otree, with.probs=TRUE, with.errors=FALSE,
                                   edge.weights=edge.weights)
    ran.idx <- sample(1:nrow(distr), size=N, prob=distr$Prob, replace=TRUE)
    ran.data0 <- distr[ran.idx, 2:otree$nmut]
    rownames(ran.data0) <- 1:N    
    
    epos <- otree$eps["epos"]
    eneg <- otree$eps["eneg"]
    ran.data <- matrix(rbinom(prod(dim(ran.data0)), size=1, prob=ifelse(ran.data0==0,epos,1-eneg)),
                       nrow=nrow(ran.data0), ncol=ncol(ran.data0), 
                       dimnames=dimnames(ran.data0))
  }
  
  if (method=="S"){
    #simulate based on definition, then "corrupt" the sample
    plotinfo <- build.plot(otree$parent) #need level info
    maxLev <- max(plotinfo$level)
    level.list <- split(otree$parent$child, plotinfo$level)
    ran.data0 <- matrix(NA, nrow=N, ncol=otree$nmut)
    colnames(ran.data0) <- otree$parent$child
    w <- if (edge.weights=="observed") otree$parent$obs.weight else otree$parent$est.weight
    names(w) <- otree$parent$child
    ran.data0[,1] <- 1  #Root is always present
    
    for (lev in 2:maxLev){
      events <- level.list[[lev]]
      for (ev in events){
        ev.parent <- otree$parent$parent[otree$parent$child==ev]
        transition <- rbinom(N, size=1, prob=w[ev])
        ran.data0[,ev] <- ran.data0[,ev.parent] * transition #only happens, if parent already occurred
      }    
    }   
    ran.data0 <- ran.data0[,-1]  #remove root
    
    epos <- otree$eps["epos"]
    eneg <- otree$eps["eneg"]
    ran.data <- matrix(rbinom(prod(dim(ran.data0)), size=1, prob=ifelse(ran.data0==0,epos,1-eneg)),
                       nrow=nrow(ran.data0), ncol=ncol(ran.data0), 
                       dimnames=dimnames(ran.data0))
    
  }
	ran.data
}

