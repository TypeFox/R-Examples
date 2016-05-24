minxent.multiple<-function (q, G, eta,lambda) 
{
fk<-function (lambda, q, G, eta) 
{
     lambda0<-log(sum(q*exp(-lambda%*%G[-1,])))
     (q * exp(-lambda0)*exp(-lambda %*% G[-1,]))%*% t(G) - eta
}

    lambda0 <- log(sum(q * exp(-lambda * G[-1, ])))
    repeat {
        lambda_old <- lambda
        lambda0 <- log(sum(q * exp(-lambda_old %*% G[-1, ])))
        f_old = fk(lambda_old, q = q, G = G, eta = eta)
        dev.ent <- (q * exp(-lambda0) * exp(-lambda_old %*% G[-1, 
            ])) %*% t(G)
        pg <- c(q * exp(-lambda0) * exp(-lambda_old %*% G[-1, 
            ]))
        cov.ent <- cov.wt(t(G), wt = pg, method = "ML", cor = FALSE)
        hess.ent <- cov.ent$cov
        W <- hess.ent[-1, -1]
        Go <- solve(W)
        lambda <- lambda_old + f_old[, -1] %*% Go
        
        if (max(abs(lambda - lambda_old)) < 1e-08) 
            break
    }
    out.lj <- list(estimates = lambda, infmat = Go)
    lambda <- out.lj$estimates
    lambda0 <- log(sum(q * exp(-lambda %*% G[-1, ])))
    pi_solve <- (q * exp(-lambda0) * exp(-lambda %*% G[-1, ]))
    list(Langrangians= c(lambda0,lambda) , Estimates=pi_solve)
}
