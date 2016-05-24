################################################
###
### array operations
###
### FIXME: NO check that the dimension match!!
###
################################################

## extendDomain
## Example: a1 is array defined over universe (A,B)
## We create new array a2 defined over (A,B,C,D) such that
## a2[a,b,c,d]=a1[a,b] for all (c,d):
## a1 <- parray(c("A","B"),levels=c(2,2), values=1:4)
## a2 <- extendDomain(a1, list(C=1:2, D=1:2))

arrayExtendDomain <- function(aa, bb){
  ans <- aa
  for (ii in 1:length(bb)){
    cbb <- bb[ii]
    ans <- .arrayExpand(ans, cbb)
  }
  return(ans)
}

arrayCombine <- function(aa.list, aux){

  vn.list<-lapply(aa.list, function(zz) names(dimnames(zz)))
  dn.list<-lapply(aa.list, function(zz) dimnames(zz))

  com.vn  <- unique(unlist(vn.list))
  com.dn  <- unlist(dn.list, recursive=FALSE)[com.vn]
  #str(com.vn)
  aa.list.com <-
    lapply(1:length(aa.list),
           function(ii){
             ll <- com.dn[-match(vn.list[[ii]], com.vn)]
             if (length(ll)>0){
                 ## Udvid de enkelte tabeller til at have samme domæne...
                 zzz <- aa.list[[ii]]
                 ## cat(sprintf("table names: %s extra: %s\n",
                 ##             toString(names(dimnames(zzz))),
                 ##             toString(names(ll))))
                 arrayExtendDomain(zzz, ll)
             } else {
                 ## cat("SOMETHING ELSE...\n")
                 aa.list[[ii]]
             }
           })

  ans <- .arrayListExpand(aa.list.com, aux[1])
  return(ans)
}


.arrayExpand <- function(aa, aux, value=NULL){

    auxn <- names(aux)[1]
    auxl <- aux[[1]]
    ## cat(sprintf("on entry: auxn: %s, auxl: %s\n", toString(auxn), toString(auxl)))

    dn1 <- dimnames(aa)
    vn1 <- names(dn1)
    #cat(sprintf("vn1=%s\n", toString(vn1)))

    aa <- rep(aa, length(auxl))
    ans <- parray(c(vn1, auxn), c(dn1, list(auxn=auxl)), aa)

    if (!is.null(value))
        ans[] <- value

    ## cat(sprintf("on exit: names: %s\n", toString(names(dimnames(ans)))))
    return(ans)
}


## aa.list : list of arrays where all arrays must have the same dimnames
## (and levels) but the permutations of the arrays must not necessarily
## be the same.
.arrayListExpand <- function(aa.list, aux, value=NULL){

    auxn <- names(aux)[1]
    auxl <- aux[[1]]
    ##    cat(sprintf("on entry: auxn: %s, auxl: %s\n", toString(auxn), toString(auxl)))

    if (is.array(aa.list))
        aa.list <- list(aa.list)

    dn1 <- dimnames(aa.list[[1]])
    vn1 <- names(dn1)
    ##cat(sprintf("vn1=%s\n", toString(vn1)))

    if (length(aa.list)>1){
        idx <- rep(1:length(aa.list), length.out=length(auxl)) # Recycle...
        aa.list  <- aa.list[idx]
        ## Check that variables are the same in all arrays.
        ## Permute arrays if necessary

        dn.aa 	<- lapply(aa.list, function(zz) names(dimnames(zz)))
        if (length(dn.aa)>1){
            for (jj in 2:length(dn.aa)){
                bb <- dn.aa[[jj]]
                if (!setequal(vn1,bb)){
                    stop("array names not identical")
                } else {
                    if (!identical(vn1,bb)){
                        aa.list[[jj]] <- tablePerm(aa.list[[jj]], vn1)
                    }
                }
            }
        }
    } else {
        aa.list <- rep(aa.list, length(auxl))
    }
    ans <- parray(c(vn1, auxn), c(dn1, list(auxn=auxl)), unlist(aa.list))

    if (!is.null(value))
        ans[] <- value

    ## cat(sprintf("on exit: names: %s\n", toString(names(dimnames(ans)))))
    return(ans)
}
