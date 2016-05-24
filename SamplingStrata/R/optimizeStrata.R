optimizeStrata <- function (errors, strata, cens = NULL, strcens = FALSE, alldomains = TRUE, 
    dom = NULL, initialStrata = nrow(strata), addStrataFactor = 0.01, 
    minnumstr = 2, iter = 20, pops = 20, mut_chance = 0.05, 
    elitism_rate = 0.2, highvalue = 1e+08, suggestions = NULL, 
    realAllocation = FALSE, writeFiles = FALSE, showPlot = TRUE) 
{
	colnames(errors) <- toupper(colnames(errors))
    colnames(strata) <- toupper(colnames(strata))
    errors$DOMAINVALUE <- as.factor(errors$DOMAINVALUE)
    erro <- split(errors, list(errors$DOMAINVALUE))
    stcamp <- split(strata, list(strata$DOM1))
    if (strcens == TRUE) {
		colnames(cens) <- toupper(colnames(cens))
#       stcens <- split(cens, list(cens$DOM1))
		k <- length(levels(as.factor(strata$DOM1)))
		stcens <- NULL
		for (i in (1:k)) {
			stcens[[i]] <- cens[cens$DOM1 == i,]
 		}
    }

    ndom <- length(levels(as.factor(strata$DOM1)))
    if (alldomains == TRUE) {
        vettsol <- NULL
        outstrata <- NULL
        for (i in 1:ndom) {
            cat("\n *** Domain : ", i, " ",as.character(errors$DOMAINVALUE[i]))
            cat("\n Number of strata : ", nrow(stcamp[[i]]) )
            erro[[i]] <- erro[[i]][, -ncol(errors)]
            cens <- NULL
            flagcens <- strcens
            if (strcens == TRUE) {
                if (nrow(stcens[[i]]) > 0) {
                  cens <- stcens[[i]]
                  flagcens = TRUE
                }
            }
            if (strcens == TRUE) {
                if (nrow(stcens[[i]]) == 0) {
                  cens <- NULL
                  flagcens = FALSE
                }
            }
            if (strcens == FALSE) {
                cens <- NULL
                flagcens = FALSE
            }
            if (nrow(stcamp[[i]]) > 0) {
                solut <- strataGenalg(errors = erro[[i]], 
                  strata = stcamp[[i]], 
                  cens = cens, 
                  strcens = flagcens, 
                  dominio = i, 
                  initialStrata = nrow(stcamp[[i]]), 
                  minnumstr, 
                  iter, 
                  pops, 
                  mut_chance, 
                  elitism_rate, 
                  addStrataFactor, 
                  highvalue, 
                  suggestions, 
                  realAllocation, 
                  writeFiles, 
                  showPlot)
				if (nrow(stcamp[[i]]) == 1) {
	              solut <- list(c(1),stcamp[[i]][c(1:grep("DOM1",colnames(stcamp[[i]])))])
                  error <- data.frame(erro[[i]])
                  strat <- data.frame(solut[[2]])
                  solut[[2]]$SOLUZ <- sum(bethel(strat, error, realAllocation = T))
                  if (solut[[2]]$SOLUZ > solut[[2]]$N) solut[[2]]$SOLUZ <- solut[[2]]$N
                }
                vettsol <- c(vettsol, solut[[1]])
				if (length(outstrata) > 0) colnames(outstrata) <- toupper(colnames(outstrata))
				colnames(solut[[2]]) <- toupper(colnames(solut[[2]]))
                outstrata <- rbind(outstrata, solut[[2]])
            }
 #           if (nrow(stcamp[[i]]) == 0 & flagcens == TRUE) {
 #               sol <- NULL
 #               t <- c(rep(1, nrow(stcens[[i]])))
 #               sol <- aggrStrata(strata = stcens[[i]], nvar = (ncol(errors) - 
 #                 2), vett = t, censiti = 1, dominio = i)
 #               sol$soluz <- sol$N
 #               vettsol <- c(vettsol, t)
 #               outstrata <- rbind(outstrata, sol)
 #           }
        }
    }
    if (alldomains == FALSE) {
        if (dom < 1 | dom > ndom) 
            stop("\nInvalid value of the indicated domain\n")
        i <- dom
        erro[[i]] <- erro[[i]][, -ncol(errors)]
		flagcens <- strcens
		if (strcens == TRUE) {
			flagcens <- TRUE
			colnames(cens) <- toupper(colnames(cens))
			censi <- cens[cens$DOM1 == i,]
			if (nrow(censi) == 0) {
				flagcens <- FALSE
				censi <- NULL
			}
		}
        solut <- strataGenalg(errors = erro[[i]], 
								strata = stcamp[[i]], 
								cens = censi, 
								strcens = flagcens, 
								dominio = i, 
								initialStrata, minnumstr, 
								iter, 
								pops, 
								mut_chance, 
								elitism_rate, 
								addStrataFactor, 
								highvalue, 
								suggestions, 
								realAllocation, 
								writeFiles, 
								showPlot)
        vettsol <- solut[[1]]
        outstrata <- solut[[2]]
    }
    colnames(outstrata) <- toupper(colnames(outstrata))
    dimens <- sum(outstrata$SOLUZ)
    cat("\n *** Sample size : ", dimens)
    cat("\n *** Number of strata : ", nrow(outstrata))
    cat("\n---------------------------")
    if (writeFiles == TRUE) {
        write.table(outstrata, file = "outstrata.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)
        cat("\n...written output to outstrata.txt")
    }
    solution <- list(indices = vettsol, aggr_strata = outstrata)
    return(solution)
}
