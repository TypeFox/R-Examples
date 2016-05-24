
Blue2Red.Col = function(n, alpha = 1){

#===============================================================================
# FUNCTION :  Blue2Red.Col
#-------------------------------------------------------------------------------
# BY: Catherine Timmermans, Institute of Statistics, Biostatistics and Actuarial Sciences, UCLouvain. 
# Last update: 8/19/2010 9:17:13 AM
# Contact: catherine.timmermans@uclouvain.be
#-------------------------------------------------------------------------------
# DESCRIPTION: 
# creates a vector of colors of length n, with color going from blue (negative part) to red  (positive part) with white color corresponding to 0.
#-------------------------------------------------------------------------------
# USAGE:
# Blue2Red.Col(n, alpha = 1)
#-------------------------------------------------------------------------------
# ARGUMENTS:
# n: integer indicating the number of colors to be generated.  
# alpha:  numeric vector of values in the range [0,1] for alpha transparency channel (0 means transparent and 1 means opaque)
#-------------------------------------------------------------------------------
# VALUE:
# a colorscale of length n, with color going from blue (negative part) to red  (positive part) with white color corresponding to 0
#-------------------------------------------------------------------------------
# WARNING/NOTE :
# The use of BUUHWE.plot requires function Blue2Red.Col.r  
#-------------------------------------------------------------------------------
# EXAMPLES:
# > Blue2Red.Col(6, alpha =0.5)
# [1] "#0000FF80" "#6A6AFF80" "#D5D5FF80"
# [4] "#FFAAAA80" "#FF555580" "#FF000080"
# >  Blue2Red.Col(5)
# [1] "#0000FFFF" "#8080FFFF" "#FFFFFFFF"
# [4] "#FF8080FF" "#FF0000FF"
#===============================================================================

    if ((n <- as.integer(n[1])) > 0) { #chek if n =nb color >1
        even.n <- n%%2 == 0           ## 0 if n is even
        k <- n%/%2                    ## floor[ middle of interval]. True value if even.n == 0
        # separate color nb in 2 parts: 
        l1 <- k + 1 - even.n      # midpoint +1  if even // midpoint if odd
        l2 <- n - k + even.n      # midpoint
        c(if (l1 > 0) hsv(h = 4/6, s = seq.int(1, ifelse(even.n, 
            0.5/k, 0), length.out = l1), v = 1, alpha = alpha),                #blue part
            if (l2 > 1) hsv(h = 0, s = seq.int(0,1, length.out = l2)[-1],   #red part
                v = 1, alpha = alpha))
    }
    else character(0)
}