
tabMarg <- function(tab, marg){
    if (!is.array(tab))
        stop("'tab' is not an array")
    if (is.numeric(marg) || is.character(marg)){
        tabMarg__(tab, marg)
    } else {
        if (class(marg)== "formula"){
            marg <- all.vars(marg[[2]])
            tabMarg__(tab, marg)
        } else {
            stop("'marg' must be character or numeric vector or a right hand sided formula")
        }

    }
}

tabPerm <- function(tab, perm){
    if (!is.array(tab))
        stop("'tab' is not an array")
    if (is.numeric(perm) || is.character(perm)){
        aperm__(tab, perm)
    } else {
        if (class(perm)== "formula"){
            perm <- all.vars( perm[[2]] )
            aperm__(tab, perm)
        } else {
            stop("'perm' must be character or numeric vector or right hand sided formula")
        }
    }
}



## Notice: tabSlice does NOT call any external c-code
tabSlice <- function(tab, slice, drop.dimensions=TRUE){
    if (!is.array(tab))
        stop("'tab' must be an array")
    if (is.null(slice))
        return( tab )

    di <- dim(tab)
    dn <- dimnames(tab)
    vn <- names(dn)

    if (length(dn)==0 || length(vn)==0)
        stop("'tabSlice' requires and array with dimnames attribute")

    dn2 <- dn
    di2 <- di

    sn  <- names(slice)
    idx <- vector("list", length(di))
    idx[] <- TRUE
    mar.idx <- match( sn, vn )
    if (any((z <- is.na(mar.idx))))
        stop("Variable(s): ", sn[z], " do not exist in table...")
    for (i in seq_along( sn )){
        k <- mar.idx[i]
        ss <- slice[[i]]
        v <- if (is.character( ss )){
            w <- match( ss, dn[[k]])
            if (any(is.na(w))){
                stop(sprintf("Level(s): %s are not in table\n",
                             toString( ss )))
            }
            w
        } else {
            if (any(is.na( ss )) || max( ss )>di[k] || min( ss )<1)
                stop(sprintf("Invalid slice %s\n",
                             toString( ss )))
            ss
        }
        idx[[k]] <- v
        dn2[[k]] <- dn2[[k]][v]
        di2[k]   <- length(v)
    }
    ##print(idx)
    out <- do.call( "[", c(list(tab), idx) )
    out
    ## ## Dimensions with only one level should not appear as a dimension
    ## j <- which(di2==1)
    ## if (drop.dimensions && length(j)>0){
    ##     dn2 <- dn2[-j]
    ##     di2 <- di2[-j]
    ## }
    ## out <- array(out, dim = di2, dimnames = dn2)

    bbb<- which(unlist(lapply(idx, is.logical)))
    dim(out) <- di[bbb]
    dimnames(out) <- dn[ bbb ]
    if (drop.dimensions)
        out <- .dropRedundantDim( out )
    class(out) <- class(tab)
    out
}

.dropRedundantDim <- function(a){
    dn2 <- attr(a,"dimnames")[dim(a)>1]
    di2 <- unlist(lapply(dn2, length), use.names = FALSE)
    dim(a) <- di2
    attr(a,"dimnames") <- dn2
    a
}


tabMult <- tabMult__
tabDiv  <- tabDiv__
tabDiv0  <- tabDiv0__
tabAdd  <- tabAdd__
tabSubt <- tabSubt__
tabEqual <- tabEqual__
tabAlign <- tabAlign__
tabExpand <- tabExpand__
tabListMult <- tabListMult__
tabListAdd <- tabListAdd__




#' ## Notice: tabSlice does NOT call any external c-code
#' tabSlice2 <- function(x, margin, level, drop.dimensions=TRUE){

#'     slice <- as.list(level)
#'     if (is.character(margin))
#'         names(slice) <- margin
#'     else
#'         names(slice) <- names(dimnames(x))[margin]
#'     tabSlice(x, slice=slice, drop.dimensions=drop.dimensions)
#' }





#' tabMarg <- function(x, margin){
#'     if (!is.array(x))
#'         stop("'x' is not an array")

#'     if (is.numeric(margin)){
#'         margin <- names(dimnames(x))[ margin ]
#'         tabMargc__(x, margin)
#'     } else {
#'         if (is.character(margin)){
#'             tabMargc__(x, margin)
#'         } else {
#'             if (class(margin)== "formula"){
#'                 margin <- all.vars(margin[[2]])
#'                 tabMargc__(x, margin)
#'             } else {
#'                 stop("'margin' must be character or numeric vector")
#'             }
#'         }
#'     }
#' }

#' tabPerm <- function(x, perm){
#'     if (!is.array(x))
#'         stop("'x' is not an array")
#'     if (is.numeric(perm)){
#'         apermi__(x, perm)
#'     } else {
#'         if (is.character(perm)){
#'             apermc__(x, perm)
#'         } else {
#'             if (inherits(perm, "formula")){
#'                 perm <- all.vars( perm[[2]] )
#'                 apermc__(x, perm)
#'             } else {
#'                 stop("'perm' must be character or numeric vector")
#'             }
#'         }
#'     }
#' }
