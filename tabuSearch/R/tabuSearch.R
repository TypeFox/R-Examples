tabuSearch <-
function(size = 10, iters = 100, objFunc = NULL, config = NULL, neigh = size, listSize = 9, nRestarts = 10,
 repeatAll = 1, verbose = FALSE){
    if (size < 2){
        stop("error: config too short!")
    }
    if (iters < 2){
        stop("error: not enough iterations!")
    }
    if (listSize >= size){
        stop("error: listSize too big!")
    }
    if (neigh > size){
        stop("error: too many neighbours!")
    }
    if (is.null(objFunc)) {
        stop("A evaluation function must be provided. See the objFunc parameter.")
    }
    if (is.null(config)) {
        config <- matrix(0, 1, size) 
        config[sample(1:size, sample(1:size, 1))] <- 1
    }
    else if (size != length(config)){
        stop("Length of the starting configuration != size")
    }
    if (repeatAll < 1){
        stop("error: repeatAll must be > 0")
    }

    iter <- 1
    configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), size)
    eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3))
    
    for(j in 1:repeatAll){
        
        if (j > 1) {
            config <- matrix(0, 1, size) 
            config[sample(1:size, sample(1:size, 1))] <- 1
        }
        
        tabuList <- matrix(0, 1, size)
        listOrder <- matrix(0, 1, listSize) #remembers order of tabu moves

        eUtility <- objFunc(config)
        aspiration <- eUtility    #highest utility found so far
        
    
        preliminarySearch<- function(){
    
            configKeep[iter, ] <- config
            eUtilityKeep[iter] <- eUtility
            iter <- iter + 1
    
            for (i in 2:iters){
                
                neighboursEUtility <- matrix(0, 1, size)  
                configTemp <- t(matrix(config, size, neigh))
                randomNeighbours <- sample(size, neigh) #pick random neighbours
                diag(configTemp[, randomNeighbours]) <- abs(diag(configTemp[, randomNeighbours]) - 1)#flip
                neighboursEUtility[randomNeighbours] <- apply(configTemp, 1, objFunc) 
                maxNontaboo <- max(neighboursEUtility[tabuList == 0])
                maxTaboo <- max(neighboursEUtility[tabuList == 1], 0)
    
                #find new move                
                move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration, 
                ifelse(length(which(neighboursEUtility == maxTaboo)) == 1, 
                       which(neighboursEUtility == maxTaboo), sample(which(neighboursEUtility == maxTaboo), 1)),  
                ifelse(length(which(neighboursEUtility == maxNontaboo & tabuList == 0)) == 1, 
                       which(neighboursEUtility == maxNontaboo & tabuList == 0), sample(which(neighboursEUtility == maxNontaboo & tabuList == 0), 1)))

                #if new utility is lower than old, add move to tabu list, else adjust aspiration, if necessary               
                if (eUtility >= neighboursEUtility[move]){
                    tabuList[move] <- 1
                    if(sum(tabuList) > listSize){ #if tabu list is full
                        tabuList[listOrder[1]] <- 0
                        listOrder[1:listSize] <- c(listOrder[2:listSize], 0)
                        } 
                    listOrder[min(which(listOrder == 0))] <- move
                    }
                else if(neighboursEUtility[move] > aspiration) aspiration <- neighboursEUtility[move]
    
                #make the new move
                eUtility <- neighboursEUtility[move]
                config[move] <- abs(config[move]-1)
                configKeep[iter,] <- config
                eUtilityKeep[iter] <- eUtility
                iter <- iter + 1
            }
    
            result = list(aspiration = aspiration,  configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter)
            return(result)
        }
        
    
        #PRELIMINARY SEARCH
        if (verbose) cat("Preliminary search stage...\n")
        result <- preliminarySearch()
        aspiration <- result$aspiration
        configKeep <- result$configKeep
        eUtilityKeep <- result$eUtilityKeep
        iter <- result$iter
    
        #INTENSIFICATION
    	  tempo <- 0
    	  restarts <- 0
        while(tempo < aspiration & restarts < nRestarts){
            if (verbose) cat("Intensification stage...\n")
            eUtility <- max(eUtilityKeep)
            tempo <- aspiration
            config <- configKeep[max(which(eUtilityKeep == max(eUtilityKeep))), ]
            result <- preliminarySearch()  
            aspiration <- result$aspiration
            configKeep <- result$configKeep
            eUtilityKeep <- result$eUtilityKeep
            iter <- result$iter
            restarts <- restarts + 1
        }
    
           	  
        #DIVERSIFICATION
        if (verbose) cat("Diversification stage...\n")
        config <- matrix(0, 1, size) 
        config[sample(1:size, sample(1:size, 1))] <- 1
        eUtility <- objFunc(config)
    
        frequent <- apply(configKeep, 2, function(x)sum(diff(x) != 0))   #create new tabu list from most frequent moves      
        tabuList <- as.numeric(rank(frequent, ties.method =  "random") > (size - listSize))
        listOrder <- sample(which(rank(frequent, ties.method =  "random") > (size - listSize)), listSize)
    
        result <- preliminarySearch()
        iter <- result$iter
        configKeep <- result$configKeep
        eUtilityKeep <- result$eUtilityKeep
    }

    endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], 
    eUtilityKeep = eUtilityKeep[1:(iter - 1)], iters = iters, neigh = neigh, 
    listSize = listSize,  repeatAll = repeatAll)
    class(endResult) = "tabu"  
    return(endResult)
}
