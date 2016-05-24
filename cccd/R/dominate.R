dominate.greedy <- function(g,weight=NULL,proportion=1.0)
{
   od <- degree(g,mode="out")+1
   S <- NULL
   A <- get.adjacency(g)
   diag(A) <- 0
   n <- nrow(A)
   covered <- rep(0,n)
   while(sum(covered)<n*proportion){
      i <- which.max(od)
		if(!is.null(weight)){
		   qq <- weight
			qq[od != od[i]] <- -Inf
			i <- which.max(qq)
		}
      covered[A[i,]>0] <- 1
		covered[i] <- 1
      S <- c(S,i)
      A[,covered>0] <- 0
		h <- graph.adjacency(A,mode="directed")
      od <- degree(h,mode="out")+1-covered
   }
   S
}

dominate.byR <- function(g,proportion=1.0)
{
   dominate.greedy(g,weight=g$R,proportion=proportion)
}

dominate.random.sample <- function(g,proportion=1.0)
{
	n <- vcount(g)
   dominate.greedy(g,weight=runif(n),proportion=proportion)
}

dominate <- function(g,method="greedy",proportion=1.0)
{
   METHODS = c("greedy","sample","byradius")
   method <- pmatch(tolower(method),METHODS)
   if(is.na(method)){
      stop("invalid method")
   }
	S <- NA
   if(method==1) S <- dominate.greedy(g,proportion=proportion)
   else if(method==2) S <- dominate.random.sample(g,proportion=proportion)
	else if(method==3) S <- dominate.byR(g,proportion=proportion)
   S
}

