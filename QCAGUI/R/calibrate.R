`calibrate` <-
function (x, type="crisp", thresholds = NA, include = TRUE, logistic = TRUE,
          idm = 0.95, ecdf = FALSE, above = 1, below = 1, ...) {
    
    if (!isNamespaceLoaded("QCA")) {
        requireNamespace("QCA", quietly = TRUE)
    }
    
    other.args <- list(...)
    
    if ("q" %in% names(other.args)) {
        above <- other.args$q
    }
    
    if ("p" %in% names(other.args)) {
        below <- other.args$p
    }
    
    if (!is.numeric(x)) {
        cat("\n")
        stop(simpleError("x is not numeric.\n\n"))
    }
    
    if (!(type %in% c("crisp", "fuzzy"))) {
        cat("\n")
        stop(simpleError("Unknown calibration type.\n\n"))
    }
    
    if (all(is.na(thresholds))) {
        cat("\n")
        stop(simpleError("Threshold value(s) not specified.\n\n"))
    }
    
    if (is.character(thresholds) & length(thresholds) == 1) {
        thresholds <- QCA::splitstr(thresholds)
    }
    
    if (type == "crisp") {
        xrange <- range(x, na.rm=TRUE)
        
        if (any(as.numeric(unclass(cut(thresholds, breaks=c(-Inf, xrange, Inf)))) != 2)) {
            cat("\n")
            stop(simpleError("Threshold value(s) outside the range of x.\n\n"))
        }
        
        if (!is.null(names(thresholds))) {
            cat("\n")
            stop(simpleError("Named thresholds require fuzzy type calibration.\n\n"))
        }
        
        return(as.numeric(unclass(cut(x, breaks=c(-Inf, thresholds, Inf), right=!include))) - 1)
        # the built-in findInterval() was interesting, but doesn't cope well with the include argument
    }
    else if (type == "fuzzy") {
        check.equal <- function(x, y) {
            check.vector <- as.logical(unlist(lapply(x, all.equal, y)))
            check.vector[is.na(check.vector)] <- FALSE
            return(check.vector)
        }
        
        lth <- length(thresholds)
        nth <- names(thresholds)
        
        
        if (lth != 3 & lth != 6) {
            cat("\n")
            stop(simpleError("For fuzzy data, there should be either 3 or 6 thresholds\".\n\n"))
        }
        
        if (idm <= 0.5 | idm >= 1) {
            cat("\n")
            stop(simpleError("The inclusion degree of membership has to be bigger than 0.5 and less than 1.\n\n"))
        }
        
        
        if (lth == 3) {
            if (!is.null(names(thresholds))) {
                if (length(unique(nth)) == sum(nth %in% c("e", "c", "i"))) {
                    thresholds <- thresholds[match(c("e", "c", "i"), nth)]
                }
            }
            
            # get rid of the names, if any
            thresholds <- as.vector(thresholds)
            
            thEX <- thresholds[1]
            thCR <- thresholds[2]
            thIN <- thresholds[3]
            
            if (logistic) {
                if (thresholds[1] > thresholds[3]) {
                    thEX <- thresholds[3]
                    thIN <- thresholds[1]
                }
                
                y <- (x < thCR) + 1
                # y is the index of the position in the vector {-1, 1}
                
                result <- 1/(1 + exp(-((x - thCR) * (c(1, -1)[y]*log(idm/(1 - idm))/(c(thIN, thEX)[y] - thCR)))))
                
                if (thresholds[1] > thresholds[3]) {
                    return(1 - result)
                }
                else {
                    return(result)
                }
            }
            else {
                if (any(table(c(thEX, thCR, thIN)) > 1)) {
                    cat("\n")
                    warning(simpleWarning("Some thresholds equal, that should not be equal.\n\n"))
                }
                
                increasing <- TRUE
                
                if (thIN < thCR & thCR < thEX) {
                    increasing <- FALSE
                }      
                
                if (ecdf) {
                    ecdfx <- x[-which(x < min(thresholds))]
                    ecdfx <- ecdfx[-which(ecdfx > max(thresholds))]
                    Fn <- ecdf(ecdfx)
                }
                
                fs <- rep(NA, length(x))    
                for (i in seq(length(x))) {
                    if (increasing) {
                        if (x[i] < thEX | check.equal(x[i], thEX)) {
                            fs[i] <- 0
                        }
                        else if (x[i] < thCR | check.equal(x[i], thCR)) {
                            fs[i] <- (((thEX - x[i])/(thEX - thCR))^below)/2
                            if (ecdf) {
                                fs[i] <- (Fn(x[i])/Fn(thCR))/2
                            }
                        }
                        else if (x[i] < thIN | check.equal(x[i], thIN)) {
                            fs[i] <- 1 - (((thIN - x[i])/(thIN - thCR))^above)/2
                            if (ecdf) {
                                fs[i] <- 1 - ((1 - Fn(x[i]))/(1 - Fn(thCR)))/2
                            }
                        }
                        else {
                            fs[i] <- 1
                        }
                    }
                    else {
                        if (x[i] < thIN | check.equal(x[i], thIN)) {
                            fs[i] <- 1
                        }
                        else if (x[i] < thCR | check.equal(x[i], thCR)) {
                            fs[i] <- 1 - (((thIN - x[i])/(thIN - thCR))^above)/2
                            if (ecdf) {
                                fs[i] <- 1 - (Fn(x[i])/Fn(thCR))/2
                            }
                        }
                        else if (x[i] < thEX | check.equal(x[i], thEX)) {
                            fs[i] <- (((thEX - x[i])/(thEX - thCR))^below)/2
                            if (ecdf) {
                                fs[i] <- ((1 - Fn(x[i]))/(1 - Fn(thCR)))/2
                            }
                        }
                        else {
                            fs[i] <- 0
                        }
                    }
                }
            }
            return(fs)
        }
        else { # 6 thresholds
            if (!is.null(nth)) {
                if (length(unique(nth)) == sum(nth %in% c("e1", "c1", "i1", "i2", "c2", "e2"))) {
                    thresholds <- thresholds[match(c("e1", "c1", "i1", "i2", "c2", "e2"), nth)]
                }
            }
            
            # get rid of the names, if any
            thresholds <- as.vector(thresholds)
            
            thEX1 <- thresholds[1]
            thCR1 <- thresholds[2]
            thIN1 <- thresholds[3]
            thIN2 <- thresholds[4]
            thCR2 <- thresholds[5]
            thEX2 <- thresholds[6]
            
            if (thCR1 < min(thEX1, thIN1) | thCR1 > max(thEX1, thIN1)) {
                cat("\n")
                stop(simpleError("First crossover threshold not between first exclusion and inclusion thresholds.\n\n"))
            }
            
            if (thCR2 < min(thEX2, thIN2) | thCR2 > max(thEX2, thIN2)) {
                cat("\n")
                stop(simpleError("Second crossover threshold not between second exclusion and inclusion thresholds.\n\n"))
            }
            
            
            somequal <- FALSE
            if (any(table(c(thEX1, thCR1, thIN1)) > 1) | any(table(c(thIN2, thCR2, thEX2)) > 1) | thCR1 == thCR2) {
                somequal <- TRUE
            }  
            
            increasing <- TRUE
            if (thIN1 < thCR1 & thCR1 < thEX1 & thEX1 <= thEX2 & thEX2 < thCR2 & thCR2 < thIN2) {
                increasing <- FALSE
            }
            
            if (increasing) {
                if (thEX1 == thEX2) {
                    somequal <- TRUE
                }
            }
            else {
                if (thIN1 == thIN2) {
                    somequal <- TRUE
                }
            }
            
            if (somequal) {
                cat("\n")
                stop(simpleError("Some thresholds equal, that should not be equal.\n\n"))
            }
            
            fs <- rep(NA, length(x))
            for (i in seq(length(x))) {
                if (increasing) {
                    if (x[i] < thEX1 | check.equal(x[i], thEX1)) {
                        fs[i] <- 0
                    }
                    else if (x[i] < thCR1 | check.equal(x[i], thCR1)) {
                        fs[i] <- (((thEX1 - x[i])/(thEX1 - thCR1))^below)/2
                    }
                    else if (x[i] < thIN1) {
                        fs[i] <- 1 - (((thIN1 - x[i])/(thIN1 - thCR1))^above)/2
                    }
                    else if (x[i] < thIN2 | check.equal(x[i], thIN2)) {
                        fs[i] <- 1
                    }
                    else if (x[i] < thCR2 | check.equal(x[i], thCR2)) {
                        fs[i] <- 1 - (((thIN2 - x[i])/(thIN2 - thCR2))^above)/2
                    }
                    else if (x[i] < thEX2 | check.equal(x[i], thEX2)) {
                        fs[i] <- (((thEX2 - x[i])/(thEX2 - thCR2))^below)/2
                    }
                    else {
                        fs[i] <- 0
                    }
                }
                else {
                    if (x[i] < thIN1 | check.equal(x[i], thIN1)) {
                        fs[i] <- 1
                    }
                    else if (x[i] < thCR1 | check.equal(x[i], thCR1)) {
                        fs[i] <- 1 - (((thIN1 - x[i])/(thIN1 - thCR1))^above)/2
                    }
                    else if (x[i] < thEX1) {
                        fs[i] <- (((thEX1 - x[i])/(thEX1 - thCR1))^below)/2
                    }
                    else if (x[i] < thEX2 | check.equal(x[i], thEX2)) {
                        fs[i] <- 0
                    }
                    else if (x[i] < thCR2 | check.equal(x[i], thCR2)) {
                        fs[i] <- (((thEX2 - x[i])/(thEX2 - thCR2))^below)/2
                    }
                    else if (x[i] < thIN2 | check.equal(x[i], thIN2)) {
                        fs[i] <- 1 - (((thIN2 - x[i])/(thIN2 - thCR2))^above)/2
                    }
                    else {
                        fs[i] <- 1
                    }
                }
            }
            return(fs)
        } 
    }
}
