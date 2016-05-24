#' Simulate trees with fixed number of sampled taxa and trait prevalance
#'
#' This function generates a tree that contains a defined number of sampled taxa for each of the two trait states.
#' @param pars a vector of parameters that describe the macroevolutionary processes of the tree.  Parameters are in order of speciation rate for trait state 0, speciation rate for trait state 1, extinction rate for trait state 0, extinction rate for trait state 1, transition rate from state 0 to state 1, transition rate from state 1 to state 0.
#' @param N0 the number of sampled taxa with trait state 0.
#' @param N1 the number of sampled taxa with trait state 1.
#' @param sampling.f a vector of sampling fraction of taxa with trait state 0 and trait state 1.
#' @param max.t the maximum amount of time, above which tree simulation stops and reports the tree as not being able to coalesce.
#' @examples
#' phy <- treesim(pars=c(0.1,0.1,0.05,0.05,0.1,0.1),N0=50,N1=50,sampling.f=c(1,1),max.t=Inf) 
treesim <- function (pars, N0, N1, sampling.f, max.t=Inf) {
    s0 <- pars[1]
    s1 <- pars[2]
    u0 <- pars[3]
    u1 <- pars[4]
    c0 <- pars[5]
    c1 <- pars[6]
    history <- matrix(NA, 1, 3)
    extinction <- numeric()
    n0 <- N0
    n1 <- N1
    n00 <- floor(N0/sampling.f[1])-N0
    n11 <- floor(N1/sampling.f[2])-N1
    idx <- n0 + n1
    lineages0 <- c(1:n0)
    lineages1 <- c((n0 + 1):(n0 + n1))
    tip.label <- c(1:idx)
    tip.state <- c(rep(0, n0), rep(1, n1))
    t <- 0
    while (n1 > 0 || n0 > 0) {
    	if (t >= max.t) {stop("Tree does not coalesce")}
        accept <- 0
        r0 <- ifelse((n0+n00) > 1, s0, 0) + ifelse((n0+n00) > 0, (u0 + c0), 
            0)
        r1 <- ifelse((n1+n11) > 1, s1, 0) + ifelse((n1+n11) > 0, (u1 + c1), 
            0)
        g <- r0 * (n0 + n00 - 1) + r1 * (n1 + n11 - 1)
        c <- r0 * (n0 + n00) + r1 * (n1 + n11)
        if (g == 0) {
            g <- c
        }
        while (accept == 0) {
            dt <- rexp(1, g)
            ps1 <- ifelse((n1 + n11) > 1, yes = (s1 * (n1 + n11 - 1) * exp(-r0 * 
                dt)), no = 0)
            pu1 <- ifelse((n1 + n11) > 0, yes = (u1 * (n1 + n11 + 1) * exp(-(2 * 
                r1 + r0) * dt)), no = 0)
            pc1 <- ifelse((n1 + n11) > 0, yes = (c1 * (n1 + n11 + 1) * exp(-2 * 
                r1 * dt)), no = 0)
            ps0 <- ifelse((n0 + n00) > 1, yes = (s0 * (n0 + n00 - 1) * exp(-r1 * 
                dt)), no = 0)
            pu0 <- ifelse((n0 + n00) > 0, yes = (u0 * (n0 + n00 + 1) * exp(-(2 * 
                r0 + r1) * dt)), no = 0)
            pc0 <- ifelse((n0 + n00) > 0, yes = (c0 * (n0 + n00 + 1) * exp(-2 * 
                r0 * dt)), no = 0)
            p0 <- ps0 + pu0 + pc0
            p1 <- ps1 + pu1 + pc1
            h <- p1 + p0
            if (runif(1) <= h/c) {
                accept <- 1
            }
        }
        t <- t + dt
        state <- as.integer(runif(1) > (p0/h))
        if (!state) {
        	type <- sample(3, 1, FALSE, c(ps0/p0, pu0/p0, pc0/p0))
        	if (type == 2||type == 3) {
        		j <- sample((n0 + n00), 1)
        		if (j<=n0) {
        			a <- 1
        		} else {
        			a <- 0
        		}
        	} else {
        		j <- sample((n0 + n00), 2)
        		if (j[1]<=n0 && j[2]<=n0) {
        			lineage <- lineages0[j]
        			a <- 1
        		} else {
        			a <- 0
        		}
        	}
        } else {
        	type <- sample(3, 1, FALSE, c(ps1/p1, pu1/p1, pc1/p1))
        	if (type == 2||type == 3) {
        		j <- sample((n1 + n11), 1)
        		if (j<=n1) {
        			a <- 1
        		} else {
        			a <- 0
        		}
        	} else {
        		j <- sample((n1 + n11), 2)
        		if (j[1]<=n1 && j[2]<=n1) {
        			lineage <- lineages1[j]
        			a <- 1
        		} else {
        			a <- 0
        		}
        	}
        }
        if (type == 1 && a == 1) {
            if (!state) {
            	if (lineage[1] %in% extinction) {
            		lineages0 <- lineages0[-j[1]]
            		extinction <- extinction[-which(extinction==lineage[1])]
            	} else if (!lineage[1] %in% extinction && lineage[2] %in% extinction) {
            		lineages0 <- lineages0[-j[2]]
            		extinction <- extinction[-which(extinction==lineage[2])]
            	} else {
            		new <- idx + 1
                	history <- rbind(history, cbind(rep(new, 2), lineages0[j], rep(t, 2)))
                	lineages0 <- lineages0[-j]
                	lineages0 <- c(lineages0, new)
                	idx <- idx + 1
                }
                n0 <- n0 - 1
            } else {
            	if (lineage[1] %in% extinction) {
            		lineages1 <- lineages1[-j[1]]
            		extinction <- extinction[-which(extinction==lineage[1])]
            	} else if (!lineage[1] %in% extinction && lineage[2] %in% extinction) {
            		lineages1 <- lineages1[-j[2]]
            		extinction <- extinction[-which(extinction==lineage[2])]
            	} else {
            		new <- idx + 1
                	history <- rbind(history, cbind(rep(new, 2), lineages1[j], rep(t, 2)))
               		lineages1 <- lineages1[-j]
                	lineages1 <- c(lineages1, new)
                	idx <- idx + 1
                }
                n1 <- n1 - 1
            }
        } else if (type == 1 && a == 0) {
        	if (!state) {
        		n00 <- n00 - 1
        	} else {
        		n11 <- n11 - 1
        	}
        } else if (type == 2 && a == 1) {
            new <- idx + 1
            if (!state) {
                lineages0 <- c(lineages0, new)
                extinction <- c(extinction, new)
                n0 <- n0 + 1
            } else {
                lineages1 <- c(lineages1, new)
                extinction <- c(extinction, new)
                n1 <- n1 + 1
            }
            idx <- idx + 1
        } else if (type == 2 && a == 0) {
        	if (!state) {
        		n00 <- n00 + 1
        	} else {
        		n11 <- n11 + 1
        	}
        } else if (type == 3 && a == 1) {
            if (!state) {
                lineages1 <- c(lineages1, lineages0[j])
                lineages0 <- lineages0[-j]
                n0 <- n0 - 1
                n1 <- n1 + 1
            }
            else {
                lineages0 <- c(lineages0, lineages1[j])
                lineages1 <- lineages1[-j]
                n0 <- n0 + 1
                n1 <- n1 - 1
            }
        } else {
        	if (!state) {
                n00 <- n00 - 1
                n11 <- n11 + 1
            }
            else {
                n00 <- n00 + 1
                n11 <- n11 - 1
            }
        }
        if ((n1 == 0 && n0 == 1 && lineages0 > (N0 + N1)) || 
            (n0 == 0 && n1 == 1 && lineages1 > (N0 + N1))) {
            n1 <- 0
            n0 <- 0
        }
    }
    history <- history[-1, ]
    n<-length(history[,1])
    if (floor(n/2)!=(n/2)) {
		history<-history[-n,]
		n<-n-1
		}
		anc<-history[seq(from=1,by=2,length.out=n/2),1]
		des<-matrix(history[,2],ncol=2,byrow=T)
		BL<-numeric(n)
		for (i in 1:(n/2)) {
			off1 <- which(history[,1]==des[i,1])
			off2 <- which(history[,1]==des[i,2])
			BL[2*i-1]<-ifelse(length(off1)==2,history[(2*i-1),3]-history[off1[1],3],history[(2*i-1),3])
			BL[2*i]<-ifelse(length(off2)==2,history[(2*i),3]-history[off2[1],3],history[(2*i),3])
		}
		history1 <- history
		for (i in 1:(n/2)) {
			history1[which(history==anc[i])] <- N0+N1+n/2+1-i
		}
		history<-history1
		BL <- rev(BL)
		history <- cbind(rev(history[,1]),rev(history[,2]))
		tree <- list(edge = history, edge.length=BL, Nnode=n/2, tip.label=as.character(tip.label),tip.state=tip.state)
		class(tree)<-"phylo"
		tree
}
