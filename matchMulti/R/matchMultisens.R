matchMultisens <- function(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, Gamma=1){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     z <- match.data[[treat.name]]
     strata.ranks <- tapply(match.data[[out.name]], match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))
     ws.1 <- 1
     d <- (sum(n.s1 + n.s2))
     ws.2 <- (n.s1 + n.s2) /d
    #Correction for  Multiple Testing
    H <- matrix(NA, n.s, 2)
    colnames(H) <- c("Constant", "Propor")
    H[, 1] <- Q.c <- ((q1- q2))*ws.1
    H[, 2] <- Q.p <- ((q1- q2))*ws.2
    k <- dim(H)[2]
    T.c <- sum(Q.c)/n.s
    T.p <- sum(Q.p)/n.s 
	E.s.c <- ((Gamma-1)/(n.s*(Gamma+1))) *sum(abs(Q.c))
	E.s.p <- ((Gamma-1)/(n.s*(Gamma+1))) *sum(abs(Q.p))
	st <- c((T.c - E.s.c), (T.p - E.s.p))
    cv <- (4*Gamma/(n.s^2*(1+Gamma)^2)) * t(H) %*% H
    dev <- st/sqrt(diag(cv)) 
    bot <- 1/sqrt(outer(diag(cv), diag(cv), "*"))
    cr <- cv * bot
    mx <- max(dev)
    p.sens.1 <- 1 - pmvnorm(lower = rep(-Inf, k), upper = rep(mx, k), corr = cr)
	p.sens.2 <- pmvnorm(lower = rep(-Inf, k), upper = rep(mx, k), corr = cr)
    pval.sens <- min(p.sens.1, p.sens.2)

    cat("Upper bound on one-sided p-value is: ", pval.sens, "\n")
    res <- list(pval=pval.sens)
    }
