### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### weightings.r
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  

# -  -  -  -  -  -  -  -  -  -  -  -  
# local weightings

# what the hell ;)
lw_tf <- function(m) {
    return(m)
}

# log'ed termfrequency
lw_logtf <- function(m) {
    return( log(m+1) )
}

# binary termfrequency
lw_bintf <- function(m) {
    return( (m>0)*1 )
}

# -  -  -  -  -  -   -  -  -  -  -  -  
# global weightings: 
# Dumais (1992), same in Nakov (2001)

# normalisation
gw_normalisation <- function(m) {
    return ( 1 / sqrt( rowSums((m*m), na.rm = TRUE) ) )
}

# inverse document frequency
# from Dumais (1992), Nakov (2001) uses log not log2
gw_idf <- function(m) {
    df = rowSums(lw_bintf(m), na.rm=TRUE)
    return ( ( log2(ncol(m)/df) + 1 ) )
}

# global frequency * inverse document frequency
# from Nakov (2001)
gw_gfidf <- function(m) {
    gf = rowSums(m, na.rm = TRUE)
    df = rowSums(lw_bintf(m), na.rm=TRUE)
    return ( gf/df )
}

# real entropy from Shannon (1948)
entropy <- function (m) {
    gf = rowSums(m, na.rm = TRUE)
    p = m / gf
    ndocs = ncol(m)
    # shannon resp. turing (there: "weight of evidence")
    # exception:
    #   iff p=0: 0*log(0) = 0
    #   this is solved by rowSums(..., na.rm=TRUE)
    entropy = - rowSums( (p*log(p)) / log(ndocs), na.rm = TRUE )
    return ( entropy )
}

# entropy as in Dumais(1992), Nakov(2001):
# global weighting = 1 + entropy
gw_entropy <- function(m) {
    return ( (1 + entropy(m)) )
}

