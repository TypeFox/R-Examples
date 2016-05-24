AS.select <-
function (obs, param, sumstats, obspar = NULL, abcmethod = abc, 
    grid = 10, inturn = TRUE, limit = ncol(sumstats), allow.none = FALSE, 
    do.err = FALSE, final.dens = FALSE, errfn = rsse, trace=TRUE, ...) 
{

argl <- list(...)

targind <- match(names(argl), "tol")
targind <- which(!is.na(targind))
margind <- match(names(argl), "method")
margind <- which(!is.na(margind))

if ((length(targind)==0)&identical(abcmethod,abc)){
        eps<-0.01
	argl$tol<-0.01
}
else{
        eps<-eval(argl[[targind]])
}

if ((length(margind)==0)&identical(abcmethod,abc)){
	argl$method<-'rejection'
}

    supp <- range(param)
    if (!is.matrix(obs) | is.data.frame(obs)) {
        obs <- matrix(obs, nrow = 1)
    }
    if (!is.matrix(param) | is.data.frame(param)) {
        param <- as.matrix(param)
    }
    if (!is.matrix(sumstats) | is.data.frame(sumstats)) {
        sumstats <- as.matrix(sumstats)
    }
    if (!is.null(obspar) | is.data.frame(obspar)) {
        if (!is.matrix(obspar)) {
            obspar <- matrix(obspar, byrow = T, ncol = ncol(param))
        }
        if (nrow(obs) != nrow(obspar)) {
            stop("Please supply observed statistics and observed parameter matrices with the same number of rows!\n")
        }
    }

    if (!length(colnames(param))) {
        colnames(param) <- paste("P", 1:ncol(param), sep = "")
    }
    if (!length(colnames(sumstats))) {
        colnames(sumstats) <- paste("C", 1:ncol(sumstats), sep = "")
    }

    nstats <- ncol(sumstats)
    limit <- min(nstats, limit)
    order <- sample(1:nstats, limit, F)
    if (trace) cat("Sumstat order for testing is:", order, "\n")
    I <- numeric(0)
    data <- obs
    argl$param=param

    for (j in 1:length(order)) {
        if ((length(I) == 0)) {
            if (trace) cat("Current subset is: empty Test adding:", order[j], "\n")
            if (allow.none) {
                index1 <- sample(param, ceiling(length(param) * 
                  eps), T)
		argl$target=data[c(I,order[j])]
		argl$sumstat=sumstats[,c(I,order[j])]
                index2 <- do.call(abcmethod,argl)

                if (is.null(index2$adj.values)) {
                  index2 <- index2$unadj.values
                }
                else {
                  index2 <- index2$adj.values
                }
                add <- AS.test(grid, index1, index2, supp)
                if (trace) {
                    if (add) {
                        cat("Significant change to ABC posterior - add\n\n")
                    } else {
                        cat("No significant change to ABC posterior - don't add\n\n")
                    }
                }
            }
            else {
                add <- TRUE
                if (trace) cat("Empty subset not allowed - add \n\n")
            }
        }
        else {
            if (trace) cat("Current subset is:", I, " Test adding:", order[j], "\n")
		argl$target=data[I]
		argl$sumstat=sumstats[,I]
		index1 <-do.call(abcmethod,argl)
            if (is.null(index1$adj.values)) {
                index1 <- index1$unadj.values
            }
            else {
                index1 <- index1$adj.values
            }
		argl$target=data[c(I,order[j])]
		argl$sumstat=sumstats[,c(I,order[j])]
		index2<-do.call(abcmethod,argl)
            if (is.null(index2$adj.values)) {
                index2 <- index2$unadj.values
            }
            else {
                index2 <- index2$adj.values
            }
            add <- AS.test(grid, index1, index2, supp)
            if (trace) {
                if (add) {
                    cat("Significant change to ABC posterior - add \n\n")
                } else {
                    cat("No significant change to ABC posterior - don't add \n\n")
                }
            }
        }
        if (add) {
            I <- c(I, order[j])
            bad <- NULL
            if (length(I) > 1) {
                if (trace) cat("Consider removing previous summaries\n")
                if (trace) cat("Current subset is:", I, " ")
                if (inturn) {
                  for (i in 1:(length(I) - 1)) {
                    if (trace) cat("Test removing:", I[i], "\n")
                    subset2 <- I
                    subset1 <- setdiff(subset2, I[i])
                    if (length(subset1) == 0) {
                      if (allow.none) {
                        i1 <- sample(param, ceiling(length(param) * 
                          eps), T)
                        argl$target=data[subset2]
                        argl$sumstat=sumstats[,subset2]
                        i2<-do.call(abcmethod,argl)
                        if (is.null(i2$adj.values)) {
                          i2 <- i2$unadj.values
                        }
                        else {
                          i2 <- i2$adj.values
                        }
                        add2 <- AS.test(grid, i1, i2, supp)
                        if (trace) {
                            if (add2) {
                                cat("Significant change to ABC posterior - keep\n\n")
                            } else {
                                cat("No significant change to ABC posterior - remove\n\n")
                            }
                        }                        
                      }
                      else {
                        add2 <- TRUE
                        if (trace) cat("Empty subset not allowed - keep\n\n")
                      }
                  }
                    else {
			argl$target=data[subset1]
			argl$sumstat=sumstats[,subset1]
			i1<-do.call(abcmethod,argl)
                      if (is.null(i1$adj.values)) {
                        i1 <- i1$unadj.values
                      }
                      else {
                        i1 <- i1$adj.values
                      }
			argl$target=data[subset2]
			argl$sumstat=sumstats[,subset2]
			i2<-do.call(abcmethod,argl)

                      if (is.null(i2$adj.values)) {
                        i2 <- i2$unadj.values
                      }
                      else {
                        i2 <- i2$adj.values
                      }
                      add2 <- AS.test(grid, i1, i2, supp)
                      if (trace) {
                          if (add2) {
                              cat("Significant change to ABC posterior - keep\n\n")
                          } else {
                              cat("No significant change to ABC posterior - remove\n\n")
                          }
                      } 
                    }
                    if (add2) {
                      I <- subset2
                    }
                    else {
                      I <- subset1
                    }
                  }
                }
                else {
                  for (i in 1:(length(I) - 1)) {
                    if (trace) cat("Test removing:", I[i], "\n")
                    subset2 <- I
                    subset1 <- setdiff(subset2, I[i])
                    if (length(subset1) == 0) {
                      if (allow.none) {
                        i1 <- sample(param, ceiling(length(param) * 
                          eps), T)
			argl$target=data[subset2]
			argl$sumstat=sumstats[,subset2]
			i2<-do.call(abcmethod,argl)
                        if (is.null(i2$adj.values)) {
                          i2 <- i2$unadj.values
                        }
                        else {
                          i2 <- i2$adj.values
                        }
                        add2 <- AS.test(grid, i1, i2, supp)
                        if (add2) {
                            cat("Significant change to ABC posterior - keep\n\n")
                        } else {
                            cat("No significant change to ABC posterior - remove\n\n")
                        }                        
                      }
                      else {
                        if (trace) cat("Empty subset not allowed - keep\n\n")
                        add2 <- TRUE
                      }
                    }
                    else {
			argl$target=data[subset1]
			argl$sumstat=sumstats[,subset1]
			i1<-do.call(abcmethod,argl)
                      if (is.null(i1$adj.values)) {
                        i1 <- i1$unadj.values
                      }
                      else {
                        i1 <- i1$adj.values
                      }
			argl$target=data[subset2]
			argl$sumstat=sumstats[,subset2]
			i2<-do.call(abcmethod,argl)

                      if (is.null(i2$adj.values)) {
                        i2 <- i2$unadj.values
                      }
                      else {
                        i2 <- i2$adj.values
                      }
                      add2 <- AS.test(grid, i1, i2, supp)
                      if (add2) {
                          cat("Significant change to ABC posterior - keep\n\n")
                      } else {
                          cat("No significant change to ABC posterior - remove\n\n")
                      } 
                    }
                    if (!add2) {
                      bad <- c(bad, I[i])
                    }
                  }
                  if (trace) cat("Removing unnecessary summaries:", bad, "\n\n")
                  I <- setdiff(I, bad)
                }
            }
        }
    }
    if (trace) cat("Selected summaries:", I, "\n")
    l <- list()
    best <- I
    if ((length(I) == 0)) {
        vals <- sample(param, ceiling(length(param) * eps), T)
    }
    else {
	argl$target=obs[I]
	argl$sumstat=sumstats[,I]
	vals<-do.call(abcmethod,argl)
        if (is.null(vals$adj.values)) {
            vals <- vals$unadj.values
        }
        else {
            vals <- vals$adj.values
        }
        if (do.err) {
            err <- errfn(vals, obspar, apply(param, 2, var))
        }
    }
    if (do.err) {
        l$err <- err
    }
    if (final.dens) {
        l$post.sample <- matrix(vals, ncol = 1)
    }

    l$best <- matrix(best,nrow=1)
    if(length(I)>0){
    	cm <- combmat(nstats, limit)
    	besti <- rep(0, length = nstats)
    	besti[I] <- 1
    	besti<-which(apply(cm,1,function(v){all(besti==v)}))
	rownames(l$best)<-besti
    }
    else{
	rownames(l$best)<-0
    }

    return(l)
}
