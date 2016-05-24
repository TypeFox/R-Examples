### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### lsa.R
### -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
### 
### 2005-11-22:
###    functions renamed to lsa(), as.textmatrix() and fold-in()
###    added warning to lsa(), when 'emtpy' singular values exist
###    chose NOT to replace solve() by 1/...
###    added a routine to as.textmatrix to convert matrices to textmatrices
### 2005-11-08: modified
###    design decision: weighting schemes will not be integrated.
###                     reason: happens before the LSA core operation / is optional
###    design decision: pre-processing also will stay part of textmatrix()
###                     for the same reasons
###    
### 2005-08-29: created.
### 

### lsa (textmatrix, dims) -> LSAspace($u, $v, $d)
### as.textmatrix (LSAspace) -> textmatrix
### fold-in (textmatrix, LSAspace) -> textmatrix

lsa <- function( x, dims=dimcalc_share() ) {
    
    # do the singular value decomposition
    SVD = svd(x)
    
    # if dims is a function, then calculate the number of dims
    if (is.function(dims)) {
        dims = dims(SVD$d)
    }
    if (dims < 2) dims=2
    
    if (any(SVD$d<=sqrt(.Machine$double.eps))) {
        warning("[lsa] - there are singular values which are zero.");
    }
    
    # prepare for returnation
    space = NULL
    space$tk = SVD$u[,1:dims]
    space$dk = SVD$v[,1:dims]
    space$sk = SVD$d[1:dims]
    rownames(space$tk) = rownames(x)
    rownames(space$dk) = colnames(x)
    class(space) = "LSAspace"
    
    # return the LSA space
    return ( space )
    
}

# as.textmatrix: 
# - when given an LSAspace, recalc a textmatrix of 
#   the original format, name it and return it
# - when given a normal matrix, return a textmatrix

as.textmatrix <- function (LSAspace) {
    
    if (inherits(LSAspace,"LSAspace")) {
        
        # convert an lsa-space to a textmatrix
        Y = LSAspace$tk %*% diag(LSAspace$sk) %*% t(LSAspace$dk)
        rownames(Y)=rownames(LSAspace$tk)
        colnames(Y)=rownames(LSAspace$dk)
        class(Y) = "textmatrix"
        environment(Y) = new.env()
        return(Y)
        
    } else if (inherits(LSAspace, "matrix")) {
        
        # convert a matrix to a textmatrix
        class(LSAspace) = "textmatrix"
        environment(LSAspace) = new.env()
        return(LSAspace)
        
    } else {
        stop("[as.textmatrix] - input has to be an LSAspace (or a matrix).")
    }
    
}

fold_in <- function( docvecs, LSAspace ) {
    
    dqs = crossprod( t( crossprod(docvecs,LSAspace$tk) ), solve(diag(LSAspace$sk)) )
    ### alternative: dqs = crossprod( docvecs, crossprod(t(LSAspace$tk), solve(diag(LSAspace$sk))) )
    dtm = crossprod( t( crossprod(t(LSAspace$tk),diag(LSAspace$sk)) ), t(dqs) )
    
    rownames(dtm) = rownames(LSAspace$tk)
    colnames(dtm) = colnames(docvecs)
    
    environment(dtm) = new.env()
    class(dtm) = "textmatrix"
    
    return (dtm)
    
}
