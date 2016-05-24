dgower <- function(x, type=list()){
############ Gower distance for mixed data ###############
# Input:
# x: data matrix. Missing values are indicated as "NA"
# type: list with components,
#      cuant =  position of cuantitative variables
#      bin = position of  binary variables (asymmetric)
#      nom = position of nominal variables
# Output:
# d: distance matrix
###########################################################

dx <- dim(x)
n <- dx[1]
p <- dx[2]


ntype <- names(type)   # given names in type

if (length(type)){
    if ((!is.list(type)) || (is.null(ntype)) || any(ntype == "")){
        stop("invalid ", sQuote("type"), "; must be named list")
    }
   for (nt in ntype){
        cvec <- type[[nt]]
        if (is.numeric(cvec)){
            if (!all(1 <= cvec & cvec <=p)){
                 stop("type$", nt, " must be in 1:ncol(x)")
            }
        } else {
            stop("type$", nt, " must  contain numbers")
        }

    }
    vc <- type$cuant
    vb <- type$bin
    vn <- type$nom
    v <- c(vc, vb, vn)
    if (!all(sort(v) == 1:p)){stop("type must contain speficications for all variables")}
} else {
    stop("type must contain speficications for all variables")
}

# Now we got all variables types in vc, vb and vn

# Are there missing values??
if (any(is.na(x))){
missing <- TRUE
} else{
missing <- FALSE
}

if (missing){
d <- gowerWITHmissing(x, vc, vb, vn)
} else{
d <- gowerNOmissing(x, vc, vb, vn)
}
d <- as.dist(d)
return(d)
}
