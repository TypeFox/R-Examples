predIntNormSimultaneousProb <-
function (n, df = n - 1, n.mean = 1, K, delta.over.sigma, k, 
    m, r, rule, integrate.args.list = NULL) 
{
    prob <- switch(rule, k.of.m = {
        pred.int.norm.k.of.m.on.r.prob(n = n, df = df, n.mean = n.mean, 
            K = K, delta.over.sigma = delta.over.sigma, k = k, 
            m = m, r = r, integrate.args.list = integrate.args.list)
    }, CA = {
        pred.int.norm.CA.on.r.prob(n = n, df = df, n.mean = n.mean, 
            K = K, delta.over.sigma = delta.over.sigma, m = m, 
            r = r)
    }, Modified.CA = {
        pred.int.norm.Modified.CA.on.r.prob(n = n, df = df, n.mean = n.mean, 
            K = K, delta.over.sigma = delta.over.sigma, r = r)
    })
    prob
}
