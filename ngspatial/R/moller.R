
rand.walk = function(X, A, Z, theta, tol, minit, maxit, s2, eta.max, V, verbose)
{   
    eta.range = c(0, eta.max)
    moller$randWalk(X, A, Z, theta, tol, minit, maxit, s2, eta.range, V, verbose)
}

rand.walk.train = function(X, A, Z, theta, trainit, s2, eta.max, V)
{   
    eta.range = c(0, eta.max)
    moller$randWalkTrain(X, A, Z, theta, trainit, s2, eta.range, V)
}

Moller.train = function(X, A, Z, theta, trainit, s2, eta.max)
{
    p = length(theta)
    V.init = diag(0.01, p)
    V.init[-p, -p] = vcov(glm(Z ~ X - 1, family = binomial))
    sample = rand.walk.train(X, A, Z, theta, trainit, s2, eta.max, V.init)
    V = cov(sample)
    V
}

Moller.run = function(X, A, Z, theta, trainit, tol, minit, maxit, s2, eta.max, verbose)
{
    p = length(theta)
    mcse = numeric(p)
    coefficients = numeric(p)
    result = list()
    V = Moller.train(X, A, Z, theta, trainit, s2, eta.max)
    sample = data.frame(rand.walk(X, A, Z, theta, tol, minit, maxit, s2, eta.max, V, verbose))
    for (i in 1:p)
    {
        temp = bm(sample[, i])
        coefficients[i] = temp$est
        mcse[i] = temp$se
    }
    iter = nrow(sample)
    list(sample = sample, coefficients = coefficients, mcse = mcse, iter = iter, V = V)
}

