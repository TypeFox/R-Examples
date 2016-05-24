evaluate.weights.1 <- function(x, p, eta, theta.fix, theta.var, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon)
{
    w <- rep(1, length(x)) / length(x)
    w_ <- rep(0, length(w))

    #Q.res <- matrix(0, length(x), length(x))
    #b.res <- numeric(length(x))

    iter <- 0
    while( (Tp(x, w, p, eta, theta.fix, theta.var) / Tp(x, w_, p, eta, theta.fix, theta.var) > 1 + weights.evaluation.epsilon) && (iter < weights.evaluation.max.iter) )
    {
        ###
        Q.res <- matrix(0, length(x), length(x))
        b.res <- numeric(length(x))
        ###
        iter <- iter + 1
        for(i in 1:length(eta))
            for(j in 1:length(eta))
                if(p[i,j] != 0)
                {
                    epsilon <- eta[[i]](x, theta.fix[[i]]) - eta[[j]](x, theta.var[[i,j]])
                    jacob <- jacobian(function(theta) eta[[j]](x, theta), theta.var[[i,j]])
                    Q.res <- Q.res + p[i,j] * (epsilon * jacob) %*% svd.inverse(t(jacob) %*% diag(w) %*% jacob) %*% t(epsilon * jacob)
                    b.res <- b.res + p[i,j] * epsilon * epsilon
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
                        function(theta) Tfs(x, w, eta[[i]], eta[[j]], theta.fix[[i]], theta)
                    )$par
    }
    x <- x[which(w > support.epsilon)]
    w <- w[which(w > support.epsilon)]
    list(x = x, w = w, theta.var = theta.var)
}
