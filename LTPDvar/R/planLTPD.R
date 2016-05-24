planLTPD <-
function (N, pt, pbar, b = 0.1, cm = 1, method = c("exact", "napprox","ewma2","ewmaSK" ),
intdif = 20,lam=1)
{
intdif = floor(intdif)
method = match.arg(method)
nAPPROXL <- function(N, pt, b, pbar, cm) {
    if (N < 5)
        stop("lot size < 5")
    k0L <- function(n, pt, b) {
        g = function(n, b) 1 - qnorm(b, mean = 0, sd = 1)^2/(2 *
            n - 2)
        h = function(n, pt, b) (g(n, b)/n + qnorm(1 - pt,
            mean = 0, sd = 1)^2/(2 * n - 2))^0.5
        return((qnorm(1 - pt, mean = 0, sd = 1) - qnorm(b,
            mean = 0, sd = 1) * h(n, pt, b))/g(n, b))
    }
    alpha0 <- function(n, pt, b, pbar) pnorm((k0L(n, pt,
        b) - qnorm(1 - pbar))/((1/n) + (k0L(n, pt, b)^2/(2 *
        n - 2)))^0.5)
    fF = function(n, pt, b, pbar, cm) (cm - alpha0(n + 1,
        pt, b, pbar))/(alpha0(n, pt, b, pbar) - alpha0(n +
        1, pt, b, pbar)) + n
    n0 = (ceiling(uniroot(function(n) fF(n, pt, b, pbar,
        cm) - N, c(5, N))$root))
    return(new("ACSPlan", n = n0, k = k0L(n0, pt, b)))
}
plan0 = nAPPROXL(N, pt, b, pbar, cm)
if (method == "napprox")
    return(plan0)
if (method %in% c("exact","ewma2","ewmaSK")) {
    fOptimn <- function(N, pt, b, pbar, cm,type,lam) {
     ImsEW<- function(n, N, pt, pbar, cm = 1, b = 0.1,type,lam)  
      n * cm + (N - n) * (1 - OC(pbar,n,kEWMA(n,pt ,b,type,lam),type,lam));
        Imsf = function(n) ImsEW(n, N, pt, pbar, cm, b,type,lam)
        fMinSearch = function(nl_, nu_) ifelse(nl_ == nu_,
            nl_, ifelse(Imsf(nl_ + floor((nu_ - nl_)/2)) <=
              Imsf(nl_ + floor((nu_ - nl_)/2) + 1), fMinSearch(nl_,
              nl_ + floor((nu_ - nl_)/2)), fMinSearch(nl_ +
              floor((nu_ - nl_)/2) + 1, nu_)))
         init_ = n(nAPPROXL(N, pt, b, pbar, cm))
     nl_init_ = max(c(init_ - intdif, 1))
        #   nl_init_ = 5
        nu_init_ = min(c(N, init_ + intdif))
      #   nu_init_ = 100
        outfOptimn = fMinSearch(nl_init_, nu_init_)
        if (outfOptimn == nu_init_)
            stop(" upper search interval limit reached, increase intdif parameter")
        if (outfOptimn == nl_init_)
            stop(" lower search interval limit reached, increase intdif parameter")
        return(outfOptimn)
    }
 kEWMA <-
function(n_, pt_,b_=0.1,type,lam) {
        k1_ = uniroot(function(k_) OC(pt_,n_,k_,type,lam) - b_, c(1, 8.2), tol = .Machine$double.eps)$root
       return(k1_);
}
    n = fOptimn(N, pt, b, pbar, cm,type=method,lam)
    # return(new("ACSPlan", n = n, k = kL(n, pt, b)))
      return(new("ACSPlan", n = n, k = kEWMA(n, pt,b,type=method,lam)))
}
  }
