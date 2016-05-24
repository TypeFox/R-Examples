dSMR <- function(x, size, df, np = 32, log = FALSE)
{
  if (any(is.nan(x)) == TRUE) warning("The x argument must be numeric!", call. = FALSE)
  if (any(is.nan(size)) == TRUE) warning("The size argument must be numeric!", call. = FALSE)
  if (any(is.nan(df)) == TRUE) warning("The df argument must be numeric!", call. = FALSE)
  if (any(is.nan(size)) == FALSE){
    if(any(size<=1) == TRUE) warning("Sample size must be greater than 1!", call. = FALSE)}
  if (any(is.nan(df)) == FALSE){
    if(any(df <= 0) == TRUE) warning("Degrees of freedom argument (df) must be greater than 0!", call. = FALSE)}
  if (is.numeric(x) == FALSE) stop ("Non-numeric x argument to mathematical function", call. = FALSE)
  if (is.numeric(size) == FALSE) stop ("Non-numeric size argument to mathematical function", call. = FALSE)
  if (is.numeric(df) == FALSE) stop ("Non-numeric df argument to mathematical function", call. = FALSE)
  
  nn <- max(length(x),length(size),length(df),length(np))
  
  
  if (nn == length(x) | length(x) > nn)          xx <- cbind(x) else
    if (length(x) == 1 | length(x) < nn)         xx <- cbind(rep(x, length = nn))
  
  if (nn == length(size)| length(size) > nn)     xx <- cbind(xx, size) else
    if (length(size) == 1 | length(size) < nn)   xx <- cbind(xx,rep(size, length = nn))
  
  if (nn == length(df)| length(df) > nn)         xx <- cbind(xx, df) else
    if (length(df) == 1 | length(df) < nn)       xx <- cbind(xx,rep(df, length = nn))
  
  if (nn == length(np)| length(np) > nn)         xx <- cbind(xx, np) else
    if (length(np) == 1 | length(np) < nn)       xx <- cbind(xx,rep(np, length = nn))
  
  dtched <- function(xx) return(dMR(xx[1], xx[2], xx[3], xx[4]))
  d <- apply(xx, 1, dtched)
  if (log == TRUE) d <- log(d)
  return(d)
}

pSMR <- function(q, size, df, np = 32, lower.tail = TRUE, log.p = FALSE)
{
    if (any(is.nan(q)) == TRUE) warning("The q argument must be numeric!", call. = FALSE)
    if (any(is.nan(size)) == TRUE) warning("The size argument must be numeric!", call. = FALSE)
    if (any(is.nan(df)) == TRUE) warning("The df argument must be numeric!", call. = FALSE)
    if (any(is.nan(size)) == FALSE){
    if(any(size<=1) == TRUE) warning("Sample size must be greater than 1!", call. = FALSE)}
    if (any(is.nan(df)) == FALSE){
      if(any(df <= 0) == TRUE) warning("Degrees of freedom argument df must be greater than 0!", call. = FALSE)}
    if (is.numeric(q) == FALSE) stop ("Non-numeric q argument to mathematical function", call. = FALSE)
    if (is.numeric(size) == FALSE) stop ("Non-numeric size argument to mathematical function", call. = FALSE)
    if (is.numeric(df) == FALSE) stop ("Non-numeric df argument to mathematical function", call. = FALSE)
    
    nn <- max(length(q),length(size),length(df),length(np))
    
    if (nn == length(q) | length(q) > nn)          xx <- cbind(q) else
      if (length(q) == 1 | length(q) < nn)         xx <- cbind(rep(q, length = nn))
    
    if (nn == length(size)| length(size) > nn)     xx <- cbind(xx, size) else
      if (length(size) == 1 | length(size) < nn)   xx <- cbind(xx,rep(size, length = nn))
    
    if (nn == length(df)| length(df) > nn)         xx <- cbind(xx, df) else
      if (length(df) == 1 | length(df) < nn)       xx <- cbind(xx,rep(df, length = nn))
    
    if (nn == length(np)| length(np) > nn)         xx <- cbind(xx, np) else
      if (length(np) == 1 | length(np) < nn)       xx <- cbind(xx,rep(np, length = nn))
    
    dtched <- function(xx) return(pMR(xx[1], xx[2], xx[3], xx[4]))
    p <- apply(xx, 1, dtched)
    if (lower.tail == FALSE) p <- 1 - p
    if (log.p == TRUE) p <- log(p)
    return(p)
}

qSMR <- function(p, size, df, np=32, eps=1e-13, maxit=5000, lower.tail=TRUE, log.p=FALSE)
{
  if (any(is.nan(p)) == TRUE) warning("The p argument must be numeric!", call. = FALSE)
  if (any(is.nan(size)) == TRUE) warning("The size argument must be numeric!", call. = FALSE)
  if (any(is.nan(df)) == TRUE) warning("The df argument must be numeric!", call. = FALSE)
  if (any(is.nan(size)) == FALSE){
    if(any(size<=1) == TRUE) warning("Sample size must be greater than 1!", call. = FALSE)}
  if (any(is.nan(df)) == FALSE){
    if(any(df <= 0) == TRUE) warning("Degrees of freedom argument df must be greater than 0!", call. = FALSE)}
  if (is.numeric(p) == FALSE) stop ("Non-numeric p argument to mathematical function", call. = FALSE)
  if (is.numeric(size) == FALSE) stop ("Non-numeric size argument to mathematical function", call. = FALSE)
  if (is.numeric(df) == FALSE) stop ("Non-numeric df argument to mathematical function", call. = FALSE)  
  if (log.p == TRUE) {    
    if (any(p[!is.nan(p)] > 0)) warning("log(p) > 0 indicates that the probabilities p are not between 0 and 1!", call. = FALSE)
    p <- exp(p)    
  }
  if (lower.tail == FALSE) p <- 1 - p
  if (any(p[!is.nan(p)] > 1) | any(p[!is.nan(p)] < 0)) warning("Probabilities p must be between 0 and 1!", call. = FALSE)
  
  nn <- max(length(p),length(size),length(df),length(np),length(eps),length(maxit))
  
  
  if (nn == length(p) | length(p) > nn)          xx <- cbind(p) else
    if (length(p) == 1 | length(p) < nn)         xx <- cbind(rep(p, length = nn))
    
  if (nn == length(size)| length(size) > nn)     xx <- cbind(xx, size) else
    if (length(size) == 1 | length(size) < nn)   xx <- cbind(xx,rep(size, length = nn))
  
  if (nn == length(df)| length(df) > nn)         xx <- cbind(xx, df) else
    if (length(df) == 1 | length(df) < nn)       xx <- cbind(xx,rep(df, length = nn))
  
  if (nn == length(np)| length(np) > nn)         xx <- cbind(xx, np) else
    if (length(np) == 1 | length(np) < nn)       xx <- cbind(xx,rep(np, length = nn))
    
  if (nn == length(eps)| length(eps) > nn)       xx <- cbind(xx, eps) else
    if (length(eps) == 1 | length(eps) < nn)     xx <- cbind(xx,rep(eps, length = nn))
  
  if (nn == length(maxit)| length(maxit) > nn)   xx <- cbind(xx, maxit) else
    if (length(maxit) == 1 | length(maxit) < nn) xx <- cbind(xx,rep(maxit, length = nn))
  
  dtched <- function(xx) return(qMR(xx[1], xx[2], xx[3], xx[4], xx[5], xx[6]))  
  q <- apply(xx, 1, dtched)  
  return(q)
}

rSMR <- function(n, size,  df = Inf)
{
  if (is.nan(n) == TRUE) warning("The n argument must be numeric!", call. = FALSE)
  if (is.nan(size) == TRUE) warning("The size argument must be numeric!", call. = FALSE)
  if (is.nan(df) == TRUE) warning("The df argument must be numeric!", call. = FALSE)
  if (is.nan(size) == FALSE){
    if(size <= 1) warning("Sample size must be greater than 1!", call. = FALSE)}
  if (is.nan(df) == FALSE){
    if(df <= 0) warning("Degrees of freedom argument df must be greater than 0!", call. = FALSE)}
  if (is.numeric(n) == FALSE) stop ("Non-numeric n argument to mathematical function", call. = FALSE)
  if (is.numeric(size) == FALSE) stop ("Non-numeric size argument to mathematical function", call. = FALSE)
  if (is.numeric(df) == FALSE) stop ("Non-numeric df argument to mathematical function", call. = FALSE)
  if (is.nan(n) == FALSE){
    if (n < 1) warning("The number of simulations must be greater than or equal to 1", call. = FALSE)}  
  x <- rMR(n, size, df)
  return(x)
}