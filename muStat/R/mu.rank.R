`mu.rank` <-
function(x, na.last=TRUE, na.rm=Inf)
{
  if (!is.na(na.last) && na.last == "keep") 
    { na.last <- NA
      na.rm <- FALSE }                   # R compat
  
  if (nas <- !is.na(na.rm))                              # nas allowed?
    if(nas <- length(wna <- which.na(x))) {              # nas present?
        x <- x[-wna]                                     # remove nas
        if (is.inf(na.rm)) na.rm <- is.na(na.last)       # defaults
    }
    if (!nas || na.rm) {
      return(mu.rank.nna(x))
    } else {
    x <- mu.rank.nna(x)
    n <- len(x)
    r <- numeric(n+nas)
    rna <- (nas+1)/2
    r[-wna] <- ifelse1 (is.na(na.last) || (na.last), x,     x+rna)
    r[ wna] <- ifelse1 (is.na(na.last), NA, na.last, n+rna, rna  )
    return(r)
    }
}
