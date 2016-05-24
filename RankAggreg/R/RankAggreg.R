`RankAggreg` <-
function(x, k, weights=NULL, method=c("CE", "GA"), 
		distance=c("Spearman", "Kendall"), seed=NULL, maxIter = 1000, 
		convIn=ifelse(method=="CE", 7, 30), importance=rep(1,nrow(x)),
            rho=.1, weight=.25, N=10*k^2, v1=NULL,
            popSize=100, CP=.4, MP=.01, verbose=TRUE, ...)
{
    if(!is.null(seed))    
    	  set.seed(seed)	   
    method <- match.arg(method, c("CE", "GA"))
    distance <- match.arg(distance, c("Spearman", "Kendall"))
    argss <- list(...)
    x <- x[,1:k]
    orig.x <- x

    orig.imp <- importance
    importance <- importance/sum(importance) #rescale importance weights
	
    distinct <- apply(x, 1, function(y) ifelse(length(unique(y)) < k, 1, 0))
    if(sum(distinct) >= 1)
        stop("Elements of Each Row Must Be Unique")
    if(nrow(x)<2)
        stop("X must have more than 1 row")
    if(CP == 0 | MP == 0)
	  stop("Neither CP nor MP can be 0") 

    compr.list <- unique(sort(as.vector(x)))
    n <- length(compr.list)
    #compr.list <- rev(compr.list)
	#cat(compr.list)

    comp.list <- 1:n
    x <- t(apply(x,1, function(xx) match(xx,compr.list)))
		
    if(!is.null(weights)){
        weights <- weights[,1:k]
        #standardize weights:
        weights <- t(apply(weights,1,function(z){if(max(z)==min(z)) 
			rep(0, length(z)) else (z-min(z))/(max(z)-min(z))}))
	  for(i in 1:nrow(weights)) # make sure 1 is the best score for all lists
        	if(weights[i,k]!=0)
            	weights[i,] <- 1-weights[i,]
        if(dim(x)[1] != dim(weights)[1] || dim(x)[2] != dim(weights)[2])
        stop("Dimensions of x and weight matrices have to be the same")
    }

    if (k > n)
        stop("k must be smaller or equal to n") 
	  
    fyRes <- matrix(0,1,2)
    colnames(fyRes) <- c("Minimums", "Medians")	
   
    if(method=="CE"){
        v <- matrix(1/n,n,k)
        if(!is.null(v1))
            v <- as.matrix(v1)
        y <- vector("numeric")
    
        Nhat <- round(rho*N)
        if (Nhat < 5)
            stop("rho is too small")

        t <- 1
		resN <- matrix(0,N,k)
		iter <- 0
		
        repeat
        {
            cands <- matrix(.C("sampling", as.integer(N), as.double(v), as.integer(n),
						as.integer(k), as.integer(resN), PACKAGE="RankAggreg")[[5]], N)

			#mcmcProc(v, N, argss$thin, argss$burn.in, comp.list, verbose)
            
			minf <- ifelse(t!=1, min(f.y), 0)
		        
            if(distance=="Spearman")
                f.y <- spearman(x, cands, importance, weights)
            else
                f.y <- kendall(x, cands, importance, weights)

            fy <- sort(f.y, ind=TRUE)
            y[t] <- fy$x[Nhat]
            good.cand <- cands[f.y <= y[t],]
            
		if(t==1)
			fyRes[1,] <- c(min(f.y), median(f.y))
		else
			fyRes <- rbind(fyRes,c(min(f.y), median(f.y)))
            
            v <- upd.prob(good.cand, v, weight, comp.list)
            
            best.cand <- compr.list[cands[fy$ix[1],]]
            rm(cands) # clean up
            
            y.l <- paste(best.cand, sep="", collapse=",")
            
            if(verbose){     
                cat("\n", "Iteration", t, ": ",  c("Optimal value: ", min(f.y),
                        "\n Optimal List:  ", y.l, "\n"))  
		    plotUpdate(f.y, fyRes, N, method)}	

            if(minf == min(f.y))
                iter <- iter+1
            else
                iter <- 1
            
            if(iter == convIn)
                break

            t <- t + 1
							
            if (t > maxIter){
            	cat("Did not converge after ", maxIter, " iterations. Please increase sample size N\n")
			break}
        }
    } else{
        #generate initial population randomly
        cands <- matrix(0, popSize, k)
        for(i in 1:popSize)
            cands[i,] <- sample(comp.list, k)
        
        #calculate obj. fn
            if(distance=="Spearman")
                f.y <- spearman(x, cands, importance, weights)
            else
                f.y <- kendall(x, cands, importance, weights)
        
	    fyRes[1,] <- c(min(f.y), median(f.y))	
        best.cand <- compr.list[cands[which.min(f.y),]]
        bestevery <- min(f.y)    
        
        conv=FALSE
        t <- 1
        iter <- 0
        while(!conv){
            #selection probability
            minf <- min(f.y)
            p.y <- (max(f.y)+1-f.y)/sum((max(f.y)+1-f.y))
            cpy <- cumsum(p.y)
            
            #select cands for the next generation
            ind <- runif(popSize)
            ind2 <- rep(0, popSize)
            for(i in 1:popSize)
                ind2[i] <- sum(ind[i] > cpy)+1
            cands <- cands[ind2,]

            # cross-over
            pairstocross <- floor(popSize*CP/2)
            samp <- sample(1:popSize, pairstocross*2)
            pointsofcross <- sample(2:k, pairstocross, replace=TRUE)
            for(i in 1:pairstocross){
		    swap <- ifelse(pointsofcross[i] < k/2, 1:pointsofcross[i], pointsofcross[i]:k)
                for(j in swap){
                # this loop performs partially matched crossover (PMX) described in Section 10.5 of 
                # Data Mining: Concepts, Models, Methods, and Algorithms by Mehmed Kantardzic (2003)

                    t1 <- cands[samp[i],j]
                    t2 <- cands[samp[i+pairstocross],j]
                    
                    if(!is.na(t3 <- match(t2, cands[samp[i],])))
                        cands[samp[i], t3] <- t1
                    if(!is.na(t3 <- match(t1, cands[samp[i+pairstocross],])))
                        cands[samp[i+pairstocross], t3] <- t2
                                                
                    cands[samp[i], j] <- t2
                    cands[samp[i+pairstocross], j] <- t1
                }
            }              
                          
            # random mutations with probability MP
            mutations <- round(popSize*k*MP)
		
            rows <- sample(1:popSize, mutations, replace=TRUE)
            cols <- sample(1:k, mutations, replace=TRUE)
		switchWith <- sample(comp.list, mutations, replace=TRUE)

            for(i in 1:mutations){	
		    tempI <- cands[rows[i], cols[j]]
		    if(switchWith[i] %in% cands[rows[i],])
		        cands[rows[i],which(switchWith[i]==cands[rows[i],])] <- tempI	
		    cands[rows[i], cols[j]] <- switchWith[i]
		}         
            
            #calculate obj. fn
            if(distance=="Spearman")
                f.y <- spearman(x, cands, importance, weights)
            else
                f.y <- kendall(x, cands, importance, weights)

 		fyRes <- rbind(fyRes,c(min(f.y), median(f.y)))

            y.l <- paste(compr.list[cands[which.min(f.y),]], sep="", collapse=",")
            if(verbose){     
                cat("\n", "Iteration", t, ": ",  c("Optimal value: ", min(f.y),
                        "\n Optimal List:  ", y.l, "\n"))
		    plotUpdate(f.y, fyRes, popSize, method)}
            
            if(minf == min(f.y))
                iter <- iter+1
            else
                iter <- 1
            
            if(iter == convIn)
                conv=TRUE
            
            if(min(f.y) < bestevery){
                best.cand <- compr.list[cands[which.min(f.y),]]
                bestevery <- min(f.y)
            }
            t <- t+1
							
		if(t > maxIter){
			cat("Did not converge after ", maxIter, " iterations.\n")
			break}
        }
    }
    rownames(fyRes) <- paste("Iter", 1:nrow(fyRes))
    res <- list(top.list=best.cand, optimal.value=ifelse(method=="CE", fy$x[1], bestevery), 
    sample.size = ifelse(method=="CE", N, popSize), num.iter=t, method=method, distance=distance,
    importance=orig.imp, lists = orig.x, weights = weights, sample=f.y, summary = fyRes)
    class(res) <- "raggr"
    res
}



