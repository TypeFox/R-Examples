maxent.test <- function(model, obs, sub.c, nperm = 99, quick = TRUE, alpha = 0.05, plot = TRUE){
     
    # check input
    if (is.vector(obs)){
        s.names <- names(obs)
        obs <- matrix(obs, 1, length(obs) ) ; dimnames(obs) <- list("set1", s.names)
    }
    if (is.data.frame(obs)) obs <- as.matrix(obs)
    obs.names <- dimnames(obs)
    if (!is.numeric(obs)) stop("obs must only contain numeric values\n")
    if (dim(obs)[2] == 1 && dim(obs)[1] > 1) obs <- t(obs)
    if (any(is.na(obs) )) stop("no NA's allowed\n")
    if (any(obs < 0 | obs > 1)) stop("obs must contain probabilities between 0 and 1\n")
    obs <- t(apply(obs, 1, function(x) x / sum(x))) ; dimnames(obs) <- obs.names
    states <- model$states
    s.names <- dimnames(states)[[2]]
    constr <- model$constr
    prob <- model$prob
    prior <- model$prior
    n.states <- dim(states)[2]
    if (is.vector(constr)){
       means.names <- names(constr)
       constr <- matrix(constr, 1, length(constr) ) ; dimnames(constr) <- list("set1", means.names)
       prior <- matrix(prior, 1, length(prior) ) ; dimnames(prior) <- list("set1", s.names)
       prob <- matrix(prob, 1, length(prob) ) ; dimnames(prob) <- list("set1", s.names)
    }
    n.sets <- dim(constr)[1]
    if (n.sets != dim(obs)[1]) stop("number of rows in obs and constr should be equal\n")
      
    # function to compute statistic
    stat <- function(o, p, q) sum(o * log(p / q))

    # vector to store null values
    values <- rep(NA, nperm)
 
    # test1 - compute p values and stat
    if (missing(sub.c)){
       obs.stat <- stat(obs, prob, prior)
       count <- 0
       outside.ci <- FALSE
       while(!outside.ci && count < nperm){
          count <- count + 1
          prob.temp <- matrix(NA, n.sets, n.states)
             for (j in 1:n.sets){
               shuffled <- sample(1:n.states, n.states)
               states.perm <- states[, shuffled, drop = F]; colnames(states.perm) <- s.names
               constr.perm <- functcomp(t(states.perm), obs[j, , drop = F])
               prob.temp[j, ] <- maxent(constr.perm, states.perm, prior[j, ])$prob
             }
       values[count] <- stat(obs, prob.temp, prior)
       val.temp <- values[1:count]
       p.temp <- (length(val.temp[val.temp >= obs.stat]) + 1) / (length(val.temp) + 1)
       if (quick && p.temp < 1){
          ci.hi <- p.temp + 1.96 * sqrt( (p.temp * (1 - p.temp)) / count )
          ci.lo <- p.temp - 1.96 * sqrt( (p.temp * (1 - p.temp)) / count )
          outside.ci <- ci.hi <= alpha || ci.lo >= alpha
          if (outside.ci) nperm = count
        }
      }
    }

    # test2 - compute p values and stat
    else{
       if (length(sub.c) >= n.states) stop("sub.c contains as many or more elements than there are states\n") 
       if (is.character(sub.c) && !all(sub.c %in% dimnames(states)[[1]]) ) stop("sub.c does not match the names of the state attributes\n")
       if (is.character(sub.c) && !all(sub.c %in% dimnames(constr)[[2]]) ) stop("sub.c does not match the constraint names\n")
       if (is.character(sub.c)) sub.c <- which(dimnames(states)[[1]] %in% sub.c)
       if (!is.vector(sub.c) ) stop("sub.c must be a vector\n")
       prob.a <- maxent(constr[, -sub.c, drop = F], states[-sub.c, , drop = F], prior)$prob
       obs.stat <- stat(obs, prob, prob.a)
       count <- 0
       outside.ci <- FALSE
       while(!outside.ci && count < nperm){
          count <- count + 1
          prob.temp <- matrix(NA, n.sets, n.states)
             for (j in 1:n.sets){
               shuffled <- sample(1:n.states, n.states)
               states.perm <- states
               states.perm[sub.c, ] <- states[sub.c, shuffled, drop = F]; colnames(states.perm) <- s.names
               constr.perm <- functcomp(t(states.perm), obs[j, , drop = F])
               prob.temp[j, ] <- maxent(constr.perm, states.perm, prior[j, ])$prob
             }
       values[count] <- stat(obs, prob.temp, prob.a)
       val.temp <- values[1:count]
       p.temp <- (length(val.temp[val.temp >= obs.stat]) + 1) / (length(val.temp) + 1)
       if (quick && p.temp < 1){
          ci.hi <- p.temp + 1.96 * sqrt( (p.temp * (1 - p.temp)) / count )
          ci.lo <- p.temp - 1.96 * sqrt( (p.temp * (1 - p.temp)) / count )
          outside.ci <- ci.hi <= alpha || ci.lo >= alpha
          if (outside.ci) nperm = count
        }
      }
    }

     # permutational p value and ci
     values <- values[!is.na(values)]
     p.perm <- (length(values[values >= obs.stat]) + 1) / (length(values) + 1)
     p.perm.hi <- p.perm + 1.96 * sqrt( (p.perm * (1 - p.perm)) / nperm )
     p.perm.lo <- p.perm - 1.96 * sqrt( (p.perm * (1 - p.perm)) / nperm )

    # measure of fit beween obs and pred relative to prior
    opqfit <- function(o, p, q) 1 - (sum((o - p)^2) / sum((o - q)^2) )
    r2.op <- cor(as.double(prob), as.double(obs) )^2
    if (length(unique(as.double(prior))) != 1 )   r2.oq <- cor(as.double(prior), as.double(obs) )^2
    else r2.oq = 0
    fit <- opqfit(obs, prob, prior)
    if (!missing(sub.c) ){
         r2.opa <- cor(as.double(prob.a), as.double(obs) )^2
         fit.a <- opqfit(obs, prob.a, prior)
    }
    # plot results
    if (plot){
        if (missing(sub.c) ){
             par(las = 1, mfcol = c(2, 2), oma = c(0, 0, 3, 0))
             fit.text <- bquote(fit[bold(o*","*p)*"|"*bold(q)] ==.(round(fit, 3)))
             r2.op.text <- bquote(italic(r)^2 ==.(round(r2.op, 3)))
             r2.oq.text <- bquote(italic(r)^2 ==.(round(r2.oq, 3)))
             plot(as.double(obs), as.double(prob), xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "predicted probabilities", main = expression(bold(p)(bold(C)*","*~bold(q)) ) )
             abline(0, 1, col = "grey25")
             text(0.1, 0.9, fit.text, cex = 1.2, pos = 4) 
             text(0.1, 0.75, r2.op.text, cex = 1.2, pos = 4) 
             lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", col = "grey50")
             lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", col = "grey50")
             mtext("arithmetic scale", line = 0)
             plot(as.double(obs) + 1e-4, as.double(prob) + 1e-4, xlim = c(1e-4, 1), ylim = c(1e-4, 1), xlab = "observed probabilities", ylab = "predicted probabilities", log = "xy", main = expression(bold(p)(bold(C)*","*~bold(q)) ))
             abline(0, 1, col = "grey25")
             lines(x = c(0.05 + 1e-4, 0.05 + 1e-4), c(0 + 1e-4, 1 + 1e-4), lty = "dashed", col = "grey50")
             lines(x = c(1e-4, 1 + 1e-4), c(0.05+1e-4, 0.05 + 1e-4) , lty = "dashed", col = "grey50")
             mtext("log10 scale, + 1e-4", line = 0)

             plot(as.double(obs), as.double(prior), xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "prior", main = expression(bold(q) ) )
             abline(0, 1, col = "grey25")
             text(0.1, 0.9, r2.oq.text, cex = 1.2, pos = 4)
             lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", col = "grey50")
             lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", col = "grey50")
             mtext("arithmetic scale", line = 0)
             plot(as.double(obs) + 1e-4, as.double(prior) + 1e-4, xlim = c(1e-4, 1), ylim = c(1e-4, 1), xlab = "observed probabilities", ylab = "prior", log = "xy", main = expression(bold(q) ))
             abline(0, 1, col = "grey25")
             lines(x = c(0.05 + 1e-4, 0.05 + 1e-4), c(0 + 1e-4, 1 + 1e-4), lty = "dashed", col = "grey50")
             lines(x = c(1e-4, 1 + 1e-4), c(0.05+1e-4, 0.05 + 1e-4) , lty = "dashed", col = "grey50")
             mtext("log10 scale, + 1e-4", line = 0)
             mtext(expression(H[0] %->% bold(p)(bold(C)*","*~bold(q))==bold(q) ),  cex = 1.2, outer = T, las = 1, line = 1)
             p.text <- bquote(italic(P) ==.(round(p.perm, 3)))
             mtext(p.text,  cex = 1.2, outer = T, las = 1, line = -0.5)
         }



        else{
             par(las = 1, mfcol = c(2, 2), oma = c(0, 0, 3, 0))
             fit.text <- bquote(fit[bold(o*","*p)*"|"*bold(q)] ==.(round(fit, 3)))
             fit.a.text <- bquote(fit[bold(o*","*p)*"|"*bold(q)] ==.(round(fit.a, 3)))
             r2.op.text <- bquote(italic(r)^2 ==.(round(r2.op, 3)))
             r2.opa.text <- bquote(italic(r)^2 ==.(round(r2.opa, 3)))
             plot(as.double(obs), as.double(prob), xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "predicted probabilities", main = expression(bold(p)(bold(A)~union(bold(B)*","*~bold(q))) ) )
             abline(0, 1, col = "grey25")
             text(0.1, 0.9, fit.text, cex = 1.2, pos = 4) 
             text(0.1, 0.75, r2.op.text, cex = 1.2, pos = 4) 
             lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", col = "grey50")
             lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", col = "grey50")
             mtext("arithmetic scale", line = 0)
             plot(as.double(obs) + 1e-4, as.double(prob) + 1e-4, xlim = c(1e-4, 1), ylim = c(1e-4, 1), xlab = "observed probabilities", ylab = "predicted probabilities", log = "xy",  main = expression(bold(p)(bold(A)~union(bold(B)*","*~bold(q))) )  )
             abline(0, 1, col = "grey25")
             lines(x = c(0.05 + 1e-4, 0.05 + 1e-4), c(0 + 1e-4, 1 + 1e-4), lty = "dashed", col = "grey50")
             lines(x = c(1e-4, 1 + 1e-4), c(0.05+1e-4, 0.05 + 1e-4) , lty = "dashed", col = "grey50")
             mtext("log10 scale, + 1e-4", line = 0)

             plot(as.double(obs), as.double(prob.a), xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", main = expression(bold(p)(bold(A)*","*~bold(q)) ) )
             abline(0, 1, col = "grey25")
             text(0.1, 0.9, fit.a.text, cex = 1.2, pos = 4) 
             text(0.1, 0.75, r2.opa.text, cex = 1.2, pos = 4) 
             lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", col = "grey50")
             lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", col = "grey50")
             mtext("arithmetic scale", line = 0)
             plot(as.double(obs) + 1e-4, as.double(prob.a) + 1e-4, xlim = c(1e-4, 1), ylim = c(1e-4, 1), xlab = "observed probabilities", ylab = "", log = "xy", main = expression(bold(p)(bold(A)*","*~bold(q)) ) )
             abline(0, 1, col = "grey25")
             lines(x = c(0.05 + 1e-4, 0.05 + 1e-4), c(0 + 1e-4, 1 + 1e-4), lty = "dashed", col = "grey50")
             lines(x = c(1e-4, 1 + 1e-4), c(0.05+1e-4, 0.05 + 1e-4) , lty = "dashed", col = "grey50")
             mtext("log10 scale, + 1e-4", line = 0)
             mtext(expression(H[0] %->% bold(p)(bold(A)~union(bold(B))*","*~bold(q)) == bold(p)(bold(A)*","*~bold(q)) ),  cex = 1.2, outer = T, las = 1, line = 1)
             p.text <- bquote(italic(P) ==.(round(p.perm, 3)))
             mtext(p.text,  cex = 1.2, outer = T, las = 1, line = -0.5)
        }

    }

    # return results
    res <- list()
    res$fit <- fit
    if (!missing(sub.c) ){
       res$fit.a <- fit.a
       res$r2.a <- r2.opa 
    }
    res$r2 <- r2.op
    if(missing(sub.c)) res$r2.q <- r2.oq
    res$obs.stat <- obs.stat
    res$nperm <- nperm
    res$pval <- p.perm
    res$ci.pval <- c(p.perm.lo, p.perm.hi)
    return(res)
}

