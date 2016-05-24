evaluate.weights.2 <- function(x, p, eta, theta.fix, theta.var, epsilon, jacob, weights.evaluation.epsilon, weights.evaluation.max.iter, support.epsilon)
{
    w <- rep(1, length(x)) / length(x)

    T.curr <- 1
    T.prev <- 0

    iter <- 0
    while( (T.curr / T.prev > 1 + weights.evaluation.epsilon) && (iter < weights.evaluation.max.iter) )
    {
        iter <- iter + 1
        delta_ <- 0
        for(i in 1:length(eta))
            for(j in 1:length(eta))
                if(p[i,j] != 0)
                {
                    theta.var[[i,j]] <- optim(
                        par = theta.var[[i,j]],
                        function(theta) Tfs(x, w, eta[[i]], eta[[j]], theta.fix[[i]], theta)
                    )$par
                    epsilon[[i,j]] <- eta[[i]](x, theta.fix[[i]]) - eta[[j]](x, theta.var[[i,j]])
                    jacob[[i,j]] <- jacobian(function(theta) eta[[j]](x, theta), theta.var[[i,j]])
                    temp <- epsilon[[i,j]] - t(epsilon[[i,j]] * w) %*% jacob[[i,j]] %*% svd.inverse(t(jacob[[i,j]] * w) %*% jacob[[i,j]]) %*% t(jacob[[i,j]])
                    delta_ <- delta_ + p[i,j] * temp * temp
                }

        k.max <- which.max(delta_)
        k.min <- which.min(delta_)

        T.prev <- sum(w * delta_)

        op.res <- optimize(
            f = function(alpha) 
            {
                delta_ <- 0
                w.t <- w_(alpha, w, k.max, k.min)
                for(i in 1:length(eta))
                    for(j in 1:length(eta))
                        if(p[i,j] != 0)
                        {
                            temp <- epsilon[[i,j]] - t(epsilon[[i,j]] * w.t) %*% jacob[[i,j]] %*% svd.inverse(t(jacob[[i,j]] * w.t) %*% jacob[[i,j]]) %*% t(jacob[[i,j]])
                            delta_ <- delta_ + p[i,j] * temp * temp
                        } 
                sum(w.t * delta_)
            },
            interval = c(0, w[k.min]),
            maximum = TRUE)
        w <- w_(op.res$maximum, w, k.max, k.min)
        T.curr <- op.res$objective
        x <- x[which(w > support.epsilon)]
        w <- w[which(w > support.epsilon)]
    }
    list(x = x, w = w, theta.var = theta.var)
}
