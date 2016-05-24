.V.clay.12 = function(theta, u1, u2){-((theta*(1 + theta)*((u2^theta) + theta*(-(u1^theta) + (u2^theta) + .Power(u1*u2,theta)))*log(u1) + theta*(1 + theta)*((1 + theta)*(u1^theta) + theta*(-1 + (u1^theta))*(u2^theta))*log(u2) + ((u1^theta) + (u2^theta) - .Power(u1*u2,theta))*((theta^2) + (1 + theta)*log(-1 + .Power(u1,-theta) + .Power(u2,-theta))))/((theta^2)*(1 + theta)*(-(u2^theta) + (u1^theta)*(-1 + (u2^theta)))))}
.S.clay.12 = function(theta, u1, u2){((theta^2)*.Power(1 + theta,2)*(1 + 2*theta)*((u1^theta)*.Power(u2,2*theta) - .Power(u1*u2,theta))*.Power(log(u1),2) + 2*theta*.Power(1 + theta,2)*(-.Power(u1*u2,theta) + .Power(u1,2*theta)*(-1 + (u2^theta)))*log(u2) + (theta^2)*.Power(1 + theta,2)*(1 + 2*theta)*(.Power(u1,2*theta)*(u2^theta) - .Power(u1*u2,theta))*.Power(log(u2),2) + 2*theta*.Power(1 + theta,2)*log(u1)*((-1 + (u1^theta))*.Power(u2,2*theta) - .Power(u1*u2,theta) + theta*(1 + 2*theta)*.Power(u1*u2,theta)*log(u2)) + (-.Power(u2,2*theta) + 2*(u1^theta)*.Power(u2,2*theta) + .Power(u1,2*theta)*(-1 + 2*(u2^theta)) - .Power(u1*u2,theta)*(2 + .Power(u1*u2,theta)))*(.Power(theta,3) + 2*.Power(1 + theta,2)*log(-1 + .Power(u1,-theta) + .Power(u2,-theta))))/(.Power(theta,3)*.Power(1 + theta,2)*.Power((u1^theta) + (u2^theta) - .Power(u1*u2,theta),2))}
.V.gumb.12 = function(theta, u1, u2){
  a = .Power(.Power(-log(u1),theta) + .Power(-log(u2),theta),1/theta)
  a2 = a^2
  (-(theta*(.Power(-log(u1),theta)*(1 - 3*a + a2) + .Power(theta,2)*(.Power(-log(u1),theta) - .Power(-log(u2),theta)) + theta*(-1 + a)*(2*.Power(-log(u1),theta) - .Power(-log(u2),theta)))*log(-log(u1))) + (1 + theta*(-1 + a) - 3*a + a2)*(.Power(-log(u1),theta) + .Power(-log(u2),theta))*log(.Power(-log(u1),theta) + .Power(-log(u2),theta)) + theta*(theta*(.Power(-log(u1),theta) + .Power(-log(u2),theta)) - (-(theta*(-1 + a)*(.Power(-log(u1),theta) - 2*.Power(-log(u2),theta))) + .Power(theta,2)*(-.Power(-log(u1),theta) + .Power(-log(u2),theta)) + (1 - 3*a + a2)*.Power(-log(u2),theta))*log(-log(u2))))/(.Power(theta,2)*(-1 + theta + a)*(.Power(-log(u1),theta) + .Power(-log(u2),theta)))
}
.S.gumb.12 = function(theta, u1, u2){
  a = .Power(.Power(-log(u1),theta) + .Power(-log(u2),theta),1/theta)
  a2= a^2
  a3= a^3
  b = .Power(log(u1)*log(u2),theta)
  d = log(.Power(-log(u1),theta) + .Power(-log(u2),theta))
  -((.Power(theta,2)*(.Power(-log(u1),2*theta)*(2 - 2*a + a2)*a+ 2*.Power(theta,4)*b + 5*.Power(theta,3)*(-1 + a)*b + .Power(theta,2)*(.Power(-log(u1),2*theta)*a + 4*b - 9*a*b + 4*a2*b) + theta*(-3*.Power(-log(u1),2*theta)*a + 2*.Power(-log(u1),2*theta)*a2 - b + 4*a*b - 4*a2*b + a3*b))*.Power(log(-log(u1)),2) + (2 + .Power(theta,2) + theta*(-3 + 2*a) - 2*a + a2)*a*(.Power(-log(u1),2*theta) + .Power(-log(u2),2*theta) + 2*b)*.Power(d,2) + 2*theta*d*((-1 + .Power(theta,2)*(-1 + a) + 2*theta*(1 - 3*a + a2) + 4*a - 4*a2 + a3)*(.Power(-log(u1),2*theta) + .Power(-log(u2),2*theta) + 2*b) - (2 + .Power(theta,2) + theta*(-3 + 2*a) - 2*a + a2)*a*(.Power(-log(u2),2*theta) + b)*log(-log(u2))) + .Power(theta,2)*(.Power(theta,2)*(.Power(-log(u1),2*theta) + .Power(-log(u2),2*theta) + 2*b) - 2*(-1 + .Power(theta,2)*(-1 + a) + 2*theta*(1 - 3*a + a2) + 4*a - 4*a2 + a3)*(.Power(-log(u2),2*theta) + b)*log(-log(u2)) + ((2 - 2*a + a2)*a*.Power(-log(u2),2*theta) + 2*.Power(theta,4)*b + 5*.Power(theta,3)*(-1 + a)*b + .Power(theta,2)*(a*.Power(-log(u2),2*theta) + 4*b - 9*a*b + 4*a2*b) + theta*(-3*a*.Power(-log(u2),2*theta) + 2*a2*.Power(-log(u2),2*theta) - b + 4*a*b - 4*a2*b + a3*b))*.Power(log(-log(u2)),2)) - 2*theta*log(-log(u1))*((2 + .Power(theta,2) + theta*(-3 + 2*a) - 2*a + a2)*a*(.Power(-log(u1),2*theta) + b)*d + theta*((-1 + .Power(theta,2)*(-1 + a) + 2*theta*(1 - 3*a + a2) + 4*a - 4*a2 + a3)*(.Power(-log(u1),2*theta) + b)+ (-1 + theta)*(2*.Power(theta,3) + .Power(theta,2)*(-3 + 5*a) + theta*(1 - 5*a + 4*a2) + (2 - 2*a + a2)*a)*b*log(-log(u2)))))/(.Power(theta,4)*.Power(-1 + theta + a,2)*.Power(.Power(-log(u1),theta) + .Power(-log(u2),theta),2)))
}
.V.fran.12 = function(theta, u1, u2){
  (exp(theta*(2 + u2))*(1 + theta*(u1 - u2)) + exp(theta + theta*u1)*(-1 + theta*(1 + u1 - u2)) + exp(theta*(2 + u1))*(1 + theta*(-u1 + u2)) + exp(theta + theta*u2)*(-1 + theta*(1 - u1 + u2)) + exp(theta*(1 + u1 + u2))*(-1 + theta*(-2 + u1 + u2)) + exp(theta*(u1 + u2))*(1 - theta*(-1 + u1 + u2)) + exp(theta)*(1 + theta*(-1 + u1 + u2)) - exp(2*theta)*(1 + theta*(u1 + u2)))/((-1 + exp(theta))*(-exp(theta) + exp(theta + theta*u1) - exp(theta*(u1 + u2)) + exp(theta + theta*u2))*theta)
}

.S.fran.12 = function(theta, u1, u2){
  -((exp(2*theta) + exp(4*theta) + exp(2*theta*(1 + u1)) + exp(2*theta*(2 + u1)) + exp(2*theta*(1 + u2)) + exp(2*theta*(2 + u2)) + exp(2*theta*(u1 + u2)) + exp(2*theta*(1 + u1 + u2)) + exp(3*theta)*(-2 + .Power(theta,2)) + exp(theta*(3 + 2*u1))*(-2 + .Power(theta,2)) + exp(theta*(3 + 2*u2))*(-2 + .Power(theta,2)) + exp(theta + 2*theta*u1 + 2*theta*u2)*(-2 + .Power(theta,2)) - 2*exp(theta*(1 + u1 + 2*u2))*(1 + .Power(theta,2)*.Power(-1 + u1,2)) - 2*exp(theta*(3 + u1 + 2*u2))*(1 + .Power(theta,2)*.Power(-1 + u1,2)) - 2*exp(theta*(2 + u1))*(1 + .Power(theta,2)*.Power(u1,2)) - 2*exp(theta*(4 + u1))*(1 + .Power(theta,2)*.Power(u1,2)) + 2*exp(theta*(3 + u1))*(2 + .Power(theta,2)*(-1 + 2*.Power(u1,2))) + 2*exp(theta*(2 + u1 + 2*u2))*(2 + .Power(theta,2)*(1 - 4*u1 + 2*.Power(u1,2))) + 2*exp(theta*(4 + u1 + u2))*(1 + .Power(theta,2)*.Power(u1 - u2,2)) - 2*exp(theta*(1 + 2*u1 + u2))*(1 + .Power(theta,2)*.Power(-1 + u2,2)) - 2*exp(theta*(3 + 2*u1 + u2))*(1 + .Power(theta,2)*.Power(-1 + u2,2)) - 2*exp(theta*(2 + u2))*(1 + .Power(theta,2)*.Power(u2,2)) - 2*exp(theta*(4 + u2))*(1 + .Power(theta,2)*.Power(u2,2)) + 2*exp(theta*(1 + u1 + u2))*(1 + .Power(theta,2)*.Power(-1 + u1 + u2,2)) - 2*exp(theta*(3 + u1 + u2))*(1 + .Power(theta,2)*(-2 + .Power(u1,2) + u1*(2 - 6*u2) + 2*u2 + .Power(u2,2))) + 2*exp(theta*(3 + u2))*(2 + .Power(theta,2)*(-1 + 2*.Power(u2,2))) + 2*exp(theta*(2 + 2*u1 + u2))*(2 + .Power(theta,2)*(1 - 4*u2 + 2*.Power(u2,2))) - 2*exp(theta*(2 + u1 + u2))*(1 + .Power(theta,2)*(1 + .Power(u1,2) - 4*u2 + .Power(u2,2) + u1*(-4 + 6*u2))))/(.Power(-1 + exp(theta),2)*.Power(exp(theta) - exp(theta + theta*u1) + exp(theta*(u1 + u2)) - exp(theta + theta*u2),2)*.Power(theta,2)))
}

.gumb.12.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  lu1 = -log(u1)
  lu2 = -log(u2)
  a = lu1^x + lu2^x
  xm1 = x - 1
  (lu1^xm1*(xm1 + a^(1/x))*a^(-2 + 1/x)*lu2^xm1)/(exp(a^(1/x))*u1*u2)
}

.fran.12.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  
  au1 = -1 + exp(-x * u1)
  au2 = -1 + exp(-x * u2)
  au12 = au1*au2
  au  = -1 + exp(-x)
  a = au*(1 + au12/au)
  e12 = exp(-x*(u1+u2))
  
  (e12 * au12 * x) / (a^2) - (e12 * x) / a
}

.clay.12.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  u1pt = u1^(-x)
  u2pt = u2^(-x)
  (1+x)*u1pt*u2pt*((u1pt+u2pt-1)^(-1/x-2))/(u1*u2)
}

.l.gumb.12.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  lu1 = -log(u1)
  lu2 = -log(u2)
  a = lu1^x + lu2^x
  xm1 = x - 1
  sum(log((lu1^xm1*(xm1 + a^(1/x))*a^(-2 + 1/x)*lu2^xm1)/(exp(a^(1/x))*u1*u2)))
}

.l.gumb.12.density.t = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  x = 1/(1-x)
  lu1 = -log(u1)
  lu2 = -log(u2)
  a = lu1^x + lu2^x
  xm1 = x - 1
  sum(xm1 * log(lu1) + log(xm1 + a^(1/x)) + (-2 + 1/x) * log(a) + xm1*log(lu2) - a^(1/x) - log(u1*u2))
}

.l.fran.12.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  
  au1 = -1 + exp(-x * u1)
  au2 = -1 + exp(-x * u2)
  au12 = au1*au2
  au  = -1 + exp(-x)
  a = au*(1 + au12/au)
  e12 = exp(-x*(u1+u2))
  
  sum(log((e12 * au12 * x) / (a^2) - (e12 * x) / a))
}

.l.clay.12.density = function(x, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  u1pt = u1^(-x)
  u2pt = u2^(-x)
  sum(log((1+x)*u1pt*u2pt*((u1pt+u2pt-1)^(-1/x-2))/(u1*u2)))
}

.gumb.123.e2 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  lu1 = -log(u1)
  plu12t = lu1^theta + (-log(u2))^theta
  ((lu1^(theta-1))*(plu12t^(1/theta-1)))/(exp((plu12t)^(1/theta))*u1)
}

.fran.123.e2 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  (exp(theta)*(-1 + exp(theta*u2)))/(-exp(theta) + exp(theta*(1 + u1)) - exp(theta*(u1 + u2)) + exp(theta*(1 + u2)))
}

.clay.123.e2 = function(theta, u){
  if(is.null(dim(u))){u1 = u[1];u2 = u[2]}else{u1 = u[,1];u2 = u[,2]}
  u1pt = u1^theta
  u2pt = u2^theta
  -(u2pt/(u1*((-1 + 1/u1pt + 1/u2pt)^(1/theta))*(-u2pt + u1pt*(-1 + u2pt))))
}

.t.12.dens = function(x, u, nu){
  if(is.null(dim(u))){
    u1 = u[1]
    u2 = u[2]
  }else{
    u1 = u[,1]
    u2 = u[,2]
  }
  -((nu*sqrt(1 - (x^2))*(((nu + (u1^2))*(nu + (u2^2)))/(nu^2))^((1 + nu)/2)*gamma(nu/2)*gamma((2 + nu)/2))/((nu*(-1 + (x^2)) - (u1^2) + 2*x*u1*u2 - (u2^2))*((nu - nu*(x^2) + (u1^2) - 2*x*u1*u2 + (u2^2))/(nu - nu*(x^2)))^(nu/2)*(gamma((1 + nu)/2))^2))
}

.lt.12.dens = function(x, u, nu){
  if(is.null(dim(u))){
    u1 = u[1]
    u2 = u[2]
  }else{
    u1 = u[,1]
    u2 = u[,2]
  }
  sum(log(-((nu*sqrt(1 - (x^2))*(((nu + (u1^2))*(nu + (u2^2)))/(nu^2))^((1 + nu)/2)*gamma(nu/2)*gamma((2 + nu)/2))/((nu*(-1 + (x^2)) - (u1^2) + 2*x*u1*u2 - (u2^2))*((nu - nu*(x^2) + (u1^2) - 2*x*u1*u2 + (u2^2))/(nu - nu*(x^2)))^(nu/2)*(gamma((1 + nu)/2))^2))))
}

.grad.lt.12.dens = function(x, u, nu){
  if(is.null(dim(u))){
    u1 = u[1]
    u2 = u[2]
  }else{
    u1 = u[,1]
    u2 = u[,2]
  }
  sum(-(-((1 + nu)*(u1^2)*x) + u1*u2*(2 + nu + nu*(x^2)) - x*((u2^2) + nu*(-1 + u2^2 + x^2)))/((-1 + x^2)*(nu + u1^2 + u2^2 - 2*u1*u2*x - nu*(x^2))))
}

.t.cop.asFunOfX1 = function(x, Lt, x2, dft){
  pCopula(c(x, x2), tCopula(Lt, 2, df = dft, df.fixed = TRUE))
}

.gadens = function(u, sig, dims){exp(-u %*% (solve(sig) - diag(1,dims)) %*% t(t(u)) / 2)/sqrt(det(sig))}

.opt.ga = function(b, psn.sample, m, dims){
  sig = cor(psn.sample[-(((b-1)*m+1):(b*m)),])
  if(m == 1){
    result = .gadens(psn.sample[((b-1)*m+1):(b*m),], sig, dims = dims)
  }else{
    uu = psn.sample[((b-1)*m+1):(b*m),]
    result = apply(uu, 1, FUN = .gadens, sig = sig)
  }
  result 
}

.V.ga = function(u, sig.inv, dims){(sig.inv %*% (u %*% t(u) %*% sig.inv - diag(1, dims)))}
.S.ga.12 = function(u, sig){(-1 + sig[1,2]^4 - 2*u[1]*u[2]*sig[1,2]*(3 + sig[1,2]^2) + (u[1]^2)*(1 + 3*(sig[1,2]^2)) + (u[2]^2)*(1 + 3*(sig[1,2]^2)))/((-1 + sig[1,2]^2)^3)}

.V.t = function(u, rho, nu){(-rho^3 + (2 + nu)*u[1]*u[2] + nu*(rho^2)*u[1]*u[2] - rho*(-1 + (1 + nu)*u[1]^2 + (1 + nu)*u[2]^2))/((-1 + rho^2)*(-1 + rho^2 - u[1]^2 + 2*rho*u[1]*u[2] - u[2]^2))}
.S.t = function(u, rho, nu){(1 + (rho^6) - (u[1]^4) - 2*nu*(rho^5)*u[1]*u[2] + 2*(u[1]^2)*(u[2]^2) - (u[2]^4) - nu*((u[1]^2) + (u[1]^4) + (u[2]^2) + (u[2]^4)) + 4*(rho^3)*u[1]*u[2]*(-2 + (u[1]^2) + (u[2]^2) + nu*(-1 + (u[1]^2) + (u[2]^2))) + 2*rho*u[1]*u[2]*(2*(2 + (u[1]^2) + (u[2]^2)) + nu*(3 + 2*(u[1]^2) + 2*(u[2]^2))) - (rho^2)*(1 + (1 + nu)*(u[1]^4) + 2*(2 + nu)*(u[2]^2) + (1 + nu)*(u[2]^4) + 2*(u[1]^2)*(2 + nu + 7*(u[2]^2) + 5*nu*(u[2]^2))) + (rho^4)*(-1 + (4 + 3*nu)*(u[2]^2) + (u[1]^2)*(4 + nu*(3 - 2*(u[2]^2)))))/((-1 + (rho^2))^2*(1 - (rho^2) + (u[1]^2) - 2*rho*u[1]*u[2] + (u[2]^2))^2)}

