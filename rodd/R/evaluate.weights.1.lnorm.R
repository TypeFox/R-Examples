evaluate.weights.1.lnorm <- function(x, w, p, eta, sq.var, theta.fix, theta.var, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon)
{
    w <- rep(1, length(x)) / length(x)
    w_ <- rep(0, length(w))

    Q.res <- matrix(0, length(x), length(x))
    b.res <- numeric(length(x))

    iter <- 0
    while(iter < weights.evaluation.max.iter)
    {
        Q.res <- matrix(0, length(x), length(x))
        b.res <- numeric(length(x))

        iter <- iter + 1
        for(i in 1:length(eta))
            for(j in 1:length(eta))
                if(p[i,j] != 0)
                {
                    ev.eta.i <- eta[[i]](x, theta.fix[[i]]) 
                    sq.sigma.i <- log(1 + sq.var[[i]](x, theta.fix[[i]]) / (ev.eta.i * ev.eta.i))

                    ev.eta.j <- eta[[j]](x, theta.var[[i,j]]) 
                    sq.sigma.j <- log(1 + sq.var[[j]](x, theta.var[[i,j]]) / (ev.eta.j * ev.eta.j))

                    epsilon <- mu(x, eta[[i]], sq.var[[i]], theta.fix[[i]]) - mu(x, eta[[j]], sq.var[[j]], theta.var[[i,j]])
                    jacob <- jacobian(function(theta) mu(x, eta[[j]], sq.var[[j]], theta), theta.var[[i,j]])

                    S <- jacobian(
                        function(theta) 
                        {
                            t <- log(1 + sq.var[[j]](x, theta) / (eta[[j]](x, theta) * eta[[j]](x, theta))) / sq.sigma.i
                            t - log(t)
                        },
				theta.var[[i,j]]
                    )

                    temp <- epsilon * jacob / sq.sigma.i - 0.5 * S

                    Q.res <- Q.res + p[i,j] * temp %*% svd.inverse(t(jacob) %*% diag(w / sq.sigma.i) %*% jacob) %*% t(temp)
                    b.res <- b.res + p[i,j] * (log(sq.sigma.i / sq.sigma.j) + (sq.sigma.j + epsilon * epsilon) / sq.sigma.i)
                }

        w_ <- w

        w <- solve.QP(
            Dmat = nearPD(2 * Q.res)$mat, 
            dvec = b.res, meq = 1, 
            Amat = t(rbind(rep(1, length(x)), diag(rep(1, length(x))))),
            bvec = c(1, rep(0, length(x)))
        )$solution

        #Estimation of the parameters

        for(i in 1:length(eta))
            for(j in 1:length(eta))
                if(p[i,j] != 0)
                    theta.var[[i,j]] <- optim(
                        par = theta.var[[i,j]],
                        function(theta) KLD(x, w, eta[[i]], eta[[j]], sq.var[[i]], sq.var[[j]], theta.fix[[i]], theta)
                    )$par
    }
    x <- x[w > support.epsilon]
    w <- w[w > support.epsilon]
    list(x = x, w = w, theta.var = theta.var)
}