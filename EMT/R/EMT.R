.packageName <- "EMT"




"multinomial.test" <- 
function(observed, prob, useChisq = FALSE, MonteCarlo = FALSE, ntrial = 100000, atOnce = 1000000) 
{
    if(! is.vector(observed,mode="numeric")) stop(" Observations have to be stored in a vector, e.g.  'observed <- c(5,2,1)'")
    if(! is.vector(prob,mode="numeric")) stop(" Probabilities have to be stored in a vector, e.g.  'prob <- c(0.25, 0.5, 0.25)'")
    if(round(sum(prob),digits=1) != 1) stop("Wrong input: sum of probabilities must not deviate from 1.")
    if(length(observed) != length(prob)) stop(" Observations and probabilities must have same dimensions.")

    size <- sum(observed)
    groups <- length(observed)
    numEvents <- choose(size + groups - 1, groups - 1)  

    if ( MonteCarlo == FALSE ) {
        if (useChisq == FALSE) {
            res <- ExactMultinomialTest(observed, prob, size, groups, numEvents)
        } else {
            res <- ExactMultinomialTestChisquare(observed, prob, size, groups, numEvents)
        }
    } else {
        if ( ntrial < numEvents ) {
            cat(" \n WARNING: Number of simulated withdrawels is lower than the number of possible outcomes. 
                This might yield unreliable results!\n\n")}
            flush.console()
        
        if (useChisq == FALSE) {
            res <- MonteCarloMultinomialTest(observed, prob, size, groups, numEvents, ntrial, atOnce)
        } else {
            res <- MonteCarloMultinomialTestChisquare(observed, prob, size, groups, numEvents, ntrial, atOnce)
        }
    }
    invisible(res)
} 




"plotMultinom" <- 
function(listMultinom, showmax = 50) 
{
    if ((listMultinom$stat == "lowF") | (listMultinom$stat == "lowP")) {
        xlab <- "Probability of events (sorted)"
        ylab <- "Probability"
    } else {
        xlab <- "Chisquare of events (sorted)"
        ylab <- "Relative frequency"
    }

    nmax = min(showmax, length(listMultinom$allProb))
    h <- listMultinom$allProb[1:nmax]
    
    barplot(h, main = listMultinom$id, xlab = xlab, ylab = ylab, space = 1)
    mtext(paste("Size:",listMultinom$size,"   Groups:",listMultinom$groups, 
                "   p.value =",listMultinom$p.value), side = 3, col = "blue")
    if (sum(grep("Carlo", listMultinom$id))) mtext(paste(" Trials: ", listMultinom$ntrial), side = 4, col = "blue")  

    invisible(listMultinom)  			
}





"ExactMultinomialTest" <- 
function(observed, prob, size, groups, numEvents) 
{
    pObs    = dmultinom(observed, size=size, prob) 	
    eventMat <- findVectors(groups,size)    		
    if( nrow(eventMat) != numEvents ) stop("Wrong number of events calculated. \n This is probably a bug.")

    eventProb <- apply(eventMat, 1, function(x) dmultinom(x, size=size, prob=prob))  
    p.value = sum(eventProb[eventProb <= pObs])

    if(round(sum(eventProb),digits=2) != 1) stop("Wrong values for probabilities. \n This is probably a bug.")

    head <- paste("\n Exact Multinomial Test, distance measure: p\n\n")
    tab <- as.data.frame(cbind(numEvents, round(pObs, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   pObs","   p.value")
    cat(head); print(tab,row.names = FALSE)

    invisible(list(id = "Exact Multinomial Test", size = size, groups = groups, 
         stat = "lowP", allProb = sort(eventProb, decreasing = TRUE), 
         ntrial = NULL, p.value = round(p.value, digits = 4)))
} 





"ExactMultinomialTestChisquare" <- 
function(observed, prob, size, groups, numEvents) 
{
    expectedFreq  = size * prob                     	
    chi2Obs = chisqStat(observed,expectedFreq)   			

    eventMat <- findVectors(groups,size) 			
    if( nrow(eventMat) != numEvents ) stop("Wrong number of events calculated. \n This is probably a bug.")

    eventProb <- apply(eventMat, 1, function(x) dmultinom(x, size=size, prob=prob))   	
    eventChi2 <- apply(eventMat, 1, function(x) chisqStat(x,expectedFreq))   			
    eventPandChi2 <- cbind(eventProb,eventChi2)								

    if(round(sum(eventProb),digits=2) != 1) stop("Wrong values for probabilities. \n This is probably a bug.")

    p.value <- sum((eventPandChi2[eventPandChi2[,2] >= chi2Obs,])[,1])

    head <- paste("\n Exact Multinomial Test, distance measure: chisquare\n\n")
    tab <- as.data.frame(cbind(numEvents, round(chi2Obs, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   chi2Obs","   p.value")
    cat(head); print(tab,row.names = FALSE)

    invisible(list(id = "Exact Multinomial Test, Chisquare", size = size, groups = groups, 
         stat = "highChisq", allProb = sort(eventProb, decreasing = TRUE), 
         ntrial = NULL, p.value = round(p.value, digits = 4)))
} 



 
"MonteCarloMultinomialTest" <- 
function(observed, prob, size, groups, numEvents, ntrial, atOnce) 
{
    # run rmultinom for not more than "atOnce" trials at once to avoid huge arrays:
    IDofObs <- paste(observed, sep="", collapse="")
    sumLowFreq = 0
    loops = floor(ntrial/atOnce)
    rest  = ntrial%%atOnce
    if (loops > 0) {
        for (i in 1:loops) {
            res <- rmultinom(n = atOnce, size = size, prob = prob) 			
            vec <- apply(res,2,function(x) paste(x,sep="",collapse=""))    		
            frequencyTable <- table(vec)    							    				
            if (sum(rownames(frequencyTable) == IDofObs) == 0) {
                freqObs = 0   										  
            } else {
                freqObs <- frequencyTable[rownames(frequencyTable) == IDofObs][[1]]  	
            } 
            sumLowFreq  <- sumLowFreq + sum(frequencyTable[as.vector(frequencyTable) <= freqObs])
            cat(" Number of withdrawals accomplished: ",prettyNum(i*atOnce, scientific = FALSE, big.mark = "."),"\n")
            flush.console() 
        }
    }
    if (rest > 0 ) {
        res <- rmultinom(n = rest, size = size, prob = prob) 			
        vec <- apply(res,2,function(x) paste(x,sep="",collapse=""))    		
        frequencyTable <- table(vec)    							     				
        if (sum(rownames(frequencyTable) == IDofObs) == 0) {
            freqObs = 0   										  
        } else {
            freqObs <- frequencyTable[rownames(frequencyTable) == IDofObs][[1]]  	
        }   
        sumLowFreq  <- sumLowFreq + sum(frequencyTable[as.vector(frequencyTable) <= freqObs])
    }

    p.value = (sumLowFreq + 1)/(ntrial + 1) 

    head <- paste("\n Monte Carlo Multinomial Test, distance measure: f\n\n")
    tab <- as.data.frame(cbind(numEvents, round(freqObs/ntrial, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   fObs","   p.value")
    cat(head); print(tab,row.names = FALSE)

    invisible(list(id = "Monte Carlo Multinomial Test", size = size, groups = groups, 
         stat = "lowF", allProb = sort(as.vector(frequencyTable), decreasing = TRUE)/ntrial, 
         ntrial = ntrial, p.value = round(p.value, digits = 4)))
} 





"MonteCarloMultinomialTestChisquare" <- 
function(observed, prob, size, groups, numEvents, ntrial, atOnce) 
{ 
    expectedFreq  = size * prob                     			
    chi2Obs = chisqStat(observed,expectedFreq)

    # run rmultinom for not more than "atOnce" trials at once to avoid huge arrays:
    bigChis = 0
    loops = floor(ntrial/atOnce)
    rest  = ntrial%%atOnce
    if (loops > 0) {
        for (i in 1:loops) {
            res <- rmultinom(n = atOnce, size = size, prob = prob) 		
            chi2all <- apply(res,2,function(x) chisqStat(x,expectedFreq))		
            bigChis <- bigChis + length(chi2all[chi2all >= chi2Obs])     			 
            cat(" Number of withdrawals accomplished: ",prettyNum(i*atOnce, scientific = FALSE, big.mark = "."),"\n")
            flush.console()
        }
    }
    if (rest > 0 ) {
        res <- rmultinom(n = rest, size = size, prob = prob) 		
        chi2all <- apply(res,2,function(x) chisqStat(x,expectedFreq))		
        bigChis <- bigChis + length(chi2all[chi2all >= chi2Obs])     			 
    }
    if((atOnce * loops + rest) != ntrial) stop(" Number of withdrawals made incorrect.\n This is probably a bug.") 

    p.value = (bigChis + 1)/(ntrial + 1) 

    head <- paste("\n Monte Carlo Multinomial Test, distance measure: chisquare\n\n")
    tab <- as.data.frame(cbind(numEvents, round(chi2Obs, digits = 4), round(p.value, digits = 4)))
    colnames(tab) <- c("   Events","   chi2Obs","   p.value")
    cat(head); print(tab,row.names = FALSE)

    invisible(list(id = "Monte Carlo Multinomial Test, Chisquare", size = size, groups = groups, 
         stat = "highChi2", allProb = sort(as.vector(table(chi2all))/ntrial, decreasing = TRUE), 
         ntrial = ntrial, p.value = round(p.value, digits = 4)))			   
} 





"chisqStat" <- 
function(observed,expected) 
{
    chisq = sum((observed-expected)^2/expected)
    invisible(chisq)   	
}



"findVectors" <- 			
function(groups,size) 
{  
    if (groups == 1) {
        mat = size    
    } else { 
        mat <- matrix(rep(0,groups-1),nrow=1)  
        for (i in 1:size) {
            mat <- rbind(mat,findVectors(groups-1,i))  
        } 
        mat <- cbind(mat,size - rowSums(mat))   
    } 
    invisible(mat)
}



## 2010; uwe.menzel@math.uu.se   uwemenzel@gmail.com


