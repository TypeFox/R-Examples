.tau2theta = function(tau, type = "gumbel"){
  if(type == "gumbel")
    1/(1-tau)
  else if(type == "clayton")
    2*tau/ (1-tau)
}
.theta2tau = function(theta, type = "gumbel"){
  if(type == "gumbel")
    1 - 1/theta
  else if(type == "clayton")
    theta / (2 + theta)
}

.Power = function(a, b){a^b}
E = exp(1)

.getSn = function(u){
  n = dim(u)[1];
  dims = dim(u)[2];
  n/(3^dims) - sum(apply(1-u^2, 1, prod)) / (2^(dims-1)) + sum(exp(rowSums(log1p(-pmax(matrix(rep(u, each = n), ncol = dims), apply(u, 2, rep, times = n)))))) / n
}

.ksmooth2 = function(x, v, h){
  sum((1-(x[1] - v[,1])^2 / h[1]^2)^2 * (1-(x[2] - v[,2])^2 / h[2]^2)^2 * ((v[,1]-x[1])^2 <= h[1]^2) * ((v[,2]-x[2])^2 <= h[2]^2)) / (h[1] * h[2]) * 15^2 / 16^2 / dim(v)[1]
}

.integrand = function(x, Lfsample, Lbootsample, h){
  (.ksmooth2(x, Lfsample, h) - .ksmooth2(x, Lbootsample, h))^2
}
