# ----------------------------------------------------
# Optimisation of sampling units multivariate
# allocation with Genetic Algorithm together with strata
# determination 
# Author: Giulio Barcaroli
# Date: 25 September 2015
# ----------------------------------------------------
strataGenalg <- function(errors, strata, cens, strcens, 
    dominio, initialStrata, minnumstr, iter, pops, mut_chance, 
    elitism_rate, addStrataFactor, highvalue, suggestions, realAllocation,
	writeFiles, showPlot) {
    # --------------------------------------------------------------------------
    colnames(strata) <- toupper(colnames(strata))
    colnames(errors) <- toupper(colnames(errors))
    #
    # --------------------------------------------------------------------------
    nvar <- ncol(errors) - 1
    ndom <- nrow(errors)
    nstrcamp <- nrow(strata)
    if (strcens == TRUE) 
        nstrcens <- nrow(cens)
    #
    # --------------------------------------------------------------------------
    # Preparation of solution list
    v <- rep(0,nrow(strata))
	outstrata <- strata
	solution <- list(v,outstrata)
    # --------------------------------------------------------------------------  
	if (writeFiles == TRUE) {
		fileres <- paste("results", dominio, ".txt", sep = "")
		sink(file = fileres)
	}	
    cat("\n---------------------------------------------")
    cat("\nOptimal stratification with Genetic Algorithm")
    cat("\n---------------------------------------------")
    cat("\n *** Parameters ***")
    cat("\n---------------------------")
    cat("\nDomain: ", dominio)
    cat("\nMaximum number of strata: ", initialStrata)
    cat("\nMinimum number of units per stratum: ", minnumstr)
    cat("\nTake-all strata (TRUE/FALSE): ", strcens)
    if (strcens == TRUE) 
        cat("\nnumber of take-all strata : ", nstrcens)
    cat("\nnumber of sampling strata : ", nstrcamp)
    cat("\nNumber of target variables: ", nvar)
    cat("\nNumber of domains: ", ndom)
    cat("\nNumber of GA iterations: ", iter)
    cat("\nDimension of GA population: ", pops)
    cat("\nMutation chance in GA generation: ", mut_chance)
    cat("\nElitism rate in GA generation: ", elitism_rate)
    cat("\nChance to add strata to maximum: ", addStrataFactor)
    cat("\nAllocation with real numbers instead of integers: ", 
        realAllocation)
    if (!is.null(suggestions)) 
        cat("\nSuggestion: ", suggestions[1, ])
    if (writeFiles == TRUE) sink()
    #
    # --------------------------------------------------------------------------
    varloop <- c(1:nvar)
    #
    # --------------------------------------------------------------------------
    # Preparation of take-all strata
    if (strcens == TRUE) {
        vett <- c(rep(1, nrow(cens)))
        censiti <- 1
        cens <- aggrStrata(cens, nvar, vett, censiti, dominio)
        dimcens <- nrow(cens)
    }
    #
    # --------------------------------------------------------------------------
    # 
    # Evaluation function
    evaluate <- function(indices) {
        soluz <- NULL
        v <- NULL
        dimens <- NULL
        censiti <- 0
        strcor <- aggrStrata(strata, nvar, floor(indices), censiti, 
            dominio)
        dimsamp <- nrow(strcor)
        if (strcens == TRUE) 
            strcor <- rbind(strcor, cens)
        dimens <- nrow(strcor)
        # cat('\nCurrent solution:',indices) cat('\nNumber of input
        # strata:',dimens) cat('\n ...sampling strata:',dimsamp)
        # cat('\n ...take-all strata:',dimcens)
        # --------------------------------------------------------------------------
        # Chiamata funzione allocazione multivariata if (dimens <
        # initialStrata)
        soluz <- bethel(strcor, errors, minnumstr, printa = FALSE, 
            realAllocation = realAllocation)
		if (writeFiles == TRUE) {
			sink()
			sink(file = fileres, append = TRUE)
		}
        ntot <- sum(soluz)
        # if (dimens > (initialStrata-1)) ntot <- highvalue
        # cat('\nSolution: ',indices)
#        cat("\nNumber of strata:", nrow(strcor), " Sample cost:", 
#            ntot)
        # cat('\n',floor(indices))
        return(ntot)
        # print(paste('Dimensione: ',round(ntot),' Numero strata:
        # ',dimens))
    }
	evaluateMem <- memoise(evaluate)
    #
    # --------------------------------------------------------------------------
    # Monitoring of processing
    # --------------------------------------------------------------------------

    monitor <- function(obj) {
        # plot the population
        minEval <- min(obj$evaluations)
#		cat("\nSample cost:",minEval)
        if (showPlot == TRUE) plot(obj, type = "trend")
        # plot(dimens,round(rbga.results$best[iter]),type='b',main
        # = '',col='blue') title(main = list('Best sample sizes vs
        # number of strata', cex=1.5, col='red', font=2))
    }
    #
    # --------------------------------------------------------------------------
    # Genetic algorithm execution
    # --------------------------------------------------------------------------
    stringMin <- rep(1, nrow(strata))
    stringMax <- rep(initialStrata, nrow(strata))
	verb = FALSE
	show = FALSE
	if (writeFiles == TRUE) {
		verb = TRUE
		show = TRUE
	}
	    rbga.results <- rbga(stringMin, stringMax, suggestions = suggestions, 
        monitorFunc = monitor, iters = iter, popSize = pops, 
        mutationChance = mut_chance, elitism_rate, addStrataFactor, 
        evalFunc = evaluateMem, verbose = verb, showSettings = show, 
        )
    #
    # --------------------------------------------------------------------------
    # Results
    # --------------------------------------------------------------------------
	if (writeFiles == TRUE) {
#		stmt <- paste("png('plotdom", dominio, ".png',height=5, width=7, units='in', res=144)", sep = "")
		stmt <- paste("pdf('plotdom", dominio, ".pdf',height=5, width=7)", sep = "")
		eval(parse(text = stmt))
	}
    plot(rbga.results)
    title(paste("Domain #", dominio, " - Sample cost", rbga.results$best[iter]), 
        col.main = "red")
    if (writeFiles == TRUE) dev.off()
    plot(rbga.results)
    title(paste("Domain #", dominio, " - Sample cost", rbga.results$best[iter]), 
        col.main = "red")
    # summary(rbga.results, echo = TRUE)
    # print(paste('Sample size:
    # ',round(rbga.results$best[iter]))) cat(' *** Sample size:
    # ',round(rbga.results$best[iter]))
    # --------------------------------------------------------------------------
    # Writing strata corresponding to optimal solution
    # --------------------------------------------------------------------------
    v <- floor(rbga.results$population[rbga.results$evaluations == 
        min(rbga.results$evaluations), ])
    if (class(v) == "matrix") v <- as.vector(v[1, ])
	if (writeFiles == TRUE) {
		stmt <- paste("write.table(v,'solution", dominio, 
					".txt',row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)", 
					sep = "")
		eval(parse(text = stmt))
	}
    censiti <- 0
    strcor <- aggrStrata(strata, nvar, v, censiti, dominio)
 #   if (strcens == TRUE) strcor <- rbind(strcor, cens)
    soluz <- bethel(strcor, errors, minnumstr, printa = FALSE, 
        realAllocation = realAllocation)
	risulta <- cbind(strcor, soluz)
	if (writeFiles == TRUE) {
		sink()
		sink(file = fileres, append = TRUE)
		cat("\n *** Sample cost: ", sum(soluz))
		cat(paste("\n *** Number of strata: ", nrow(strcor)))
		colnames(risulta) <- toupper(colnames(risulta))
		fileout <- paste("outstrata", dominio, ".txt", sep = "")
		write.table(risulta, file = fileout, sep = "\t", row.names = FALSE, 
			col.names = TRUE, quote = FALSE)
		cat("\n...written output to", fileout)
    	sink()
	}
	solution[[1]] <- v
	solution[[2]] <- risulta
	return(solution)
    # End function
}
