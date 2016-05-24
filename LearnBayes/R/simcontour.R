simcontour=function (logf, limits, data, m) 
{
LOGF=function(theta, data)
        {
          if(is.matrix(theta)==TRUE){
          val=matrix(0,c(dim(theta)[1],1))
          for (j in 1:dim(theta)[1])
              val[j]=logf(theta[j,],data)
           } 
          else val=logf(theta,data)
          return(val)
        }
    ng = 50
    x0 = seq(limits[1], limits[2], len = ng)
    y0 = seq(limits[3], limits[4], len = ng)
    X = outer(x0, rep(1, ng))
    Y = outer(rep(1, ng), y0)
    n2 = ng^2
    Z = LOGF(cbind(X[1:n2], Y[1:n2]), data)
    Z = Z - max(Z)
    Z = matrix(Z, c(ng, ng))
    d = cbind(X[1:n2], Y[1:n2], Z[1:n2])
    dx = diff(x0[1:2])
    dy = diff(y0[1:2])
    prob = d[, 3]
    prob = exp(prob)
    prob = prob/sum(prob)
    i = sample(2500, m, replace = TRUE, prob = prob)
    return(list(x = d[i, 1] + runif(m) * dx - dx/2, y = d[i, 
        2] + runif(m) * dy - dy/2))
}
