WP.check <- function (k, WPfacs, nfac.WP, nfactors, factor.names=Letters[1:nfactors]) 
{
    if (!length(WPfacs)==nfac.WP) stop("WFfacs must be of length ", nfac.WP)
    if (!is.list(WPfacs)) {
        ## invalid blocks
        if (!(is.numeric(WPfacs) | is.character(WPfacs))) 
            stop("WPfacs must be a list of whole plot factor numbers, 
            or a character vector of whole plot factor names or factor letters.")

        ## factor names or character block generators (ABC etc.)
        if (is.character(WPfacs)) {
            ## factor names are made into numeric vector
            hilf <- factor.names
            if (is.list(hilf)) hilf <- names(hilf)
            if (all(WPfacs %in% hilf)) WPfacs <- which(hilf %in% WPfacs)
            else if (all(nchar(WPfacs)==1))
                   WPfacs <- which(Letters %in% WPfacs)   ## not foolproof
            if (is.character(WPfacs)) WPfacs <- as.numeric(chartr("F", " ", WPfacs))
            if (any(is.na(WPfacs))) stop("The character vector WPfacs must contain factor names, 
                        factor letters or factor numbers preceded by capital F.")
            }
    }
    ## generators based on factors
    ## or list of single factor numbers
    if (is.list(WPfacs)){
      if (any(!sapply(WPfacs, length)==1)) 
        stop("All entries in list WPfacs must be single factor numbers.")
      if (any(!sapply(WPfacs, is.numeric))) 
        stop("All entries in list WPfacs must be numbers.")
      if (any(sapply(WPfacs, function(obj) !(obj==floor(obj) & obj>=1 & obj < 2^k)))) 
        stop("All entries in list WPfacs must be integer numbers between 1 and 2^k-1.")
      WPfacs <- as.numeric(WPfacs)
    }
    paste("F",WPfacs,sep="")
}


#WP.check <- function (k, WPs, factor.names) 
#{
#    if (is.numeric(WPs) & length(WPs)==1) {
#            if (!(2^round(log2(WPs)))==WPs) stop ("The number WPs must be an integer power of 2.")
#            return(WPs)
#        }
#    if (!(is.numeric(WPs) | is.character(WPs))) 
#        stop ("WPs must be the number of whole plots or a vector of whole plot factors (position numbers, factor letters or factor names).")
#    if (is.character(WPs)) {
#            hilf <- factor.names
#            if (is.list(hilf)) hilf <- names(hilf)
#            if (all(WPs %in% hilf)) WPs <- which(hilf %in% WPs)
#            else {if (all(nchar(WPs)==1)) WPs <- which(Letters %in% WPs)
#                  else WPs <- as.numeric(chartr("F", " ", WPs))
#                  if (any(is.na(WPs))) stop("The character vector WPs must contain factor names, 
#                        factor letters or factor numbers preceded by capital F.")
#            }
#            }
#     if (is.numeric(WPs)){
#            if (any(!WPs == floor(WPs))) 
#                stop("All entries in WPs must be integer numbers.")
#            if (min(WPs) < 1 | max(WPs) > length(factor.names)) 
#                stop("Factor numbers in WPs must be in the range of 1 to ", length(factor.names), ".")
#         }
#     ## make WPs character
#     paste("F",WPs,sep="")
#}
