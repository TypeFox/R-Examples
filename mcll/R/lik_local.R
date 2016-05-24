# Standard error estimation
lik_local <- function(a, x, data.t, alp, nodes, weights, low, up) {

    # a : a coefficients for the polynomial function 
    # x : fitting point (length of dimension) 
    # data.t: transformed posterior 
    # alp: alpha smoothing parameter 
    # nodes: quadrature locations 
    # weights: quadrature weights at the nodes 
    # low: lower limit of the Gaussian quadrature integration 
    # up: upper limit of the Gaussian quadrature integration 
 
    ndata <- nrow(data.t)  # data size
    npar <- ncol(data.t)  # dimension

    # bandwidth 
    k = floor(ndata * alp)

    v=as.matrix(data.t)-matrix(rep(as.numeric(x), ndata),nrow=ndata, byrow=TRUE)
    vv = rowSums(v^2)
    V = sqrt(vv)
    
    dk = sort(V)
    h = dk[k]
    
    # weight    
    u = v/h 
    
    find <- function(u0) { as.numeric(all(abs(u0) <1)) }  # u0 is a vector of length npar 
    index0 <- apply(u, 1, find) # apply by row 
    index <- which(index0==1)
    
    uu <- u[index,1:npar]
    w= (1-(abs(uu))^3)^3  
    W = apply(w, 1, prod)

    # polinomial function    
    linear <- apply((t(v[index,]))*a[2:(npar+1)],2,sum)
    quad <- apply((t(v[index,]))^2*a[(1+npar+1): length(a)]*0.5 ,2,sum)
    pol <- a[1] + linear + quad 
    
    # local likelihood 
    term0 <-  W * pol
    fac1 <- max(term0)
    
    # integral 
    uu=(((up-low)/2)*nodes) + ((up+low)/2)
    index <- which(abs(uu/h) < 1) 
    
    linear2 <-  (t(matrix(rep(uu[index], npar),ncol=npar))-x)*a[2:(npar+1)] 
    quad2 <-  (t(matrix(rep(uu[index], npar),ncol=npar))-x)^2*a[(1+npar+1): length(a)]*0.5
    pol2 <- (1-(abs((t(matrix(rep(uu[index], npar),ncol=npar))-x))/h)^3)^3  * exp(linear2 + quad2)
    weight.sum <- apply(weights[index]* t(pol2), 2, sum) 
    int <- ((up-low)/2) * weight.sum

    int <- prod(int)

    term1 = sum(term0) 
    term2 = ndata * exp(a[1]) *int 
    lp = (term1 - term2 ) 
    
    return(-lp)

}
