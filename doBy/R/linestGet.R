## ######################################################################
##
## Auxillary functions for calculating contrasts, lsmeans etc
##
## Banff, August 2013, Søren Højsgaard
##
## ######################################################################

get_xlevels <- function(obj){
  UseMethod("get_xlevels")
}

get_xlevels.default <- function(obj){
  obj$xlevels
}

get_xlevels.mer <- function(obj){
  ans <- lapply(obj@frame,levels)
  ans <- ans[unlist(lapply(ans, length))>0]

  nn <- names(attr(terms(obj), "dataClasses"))
  ii <- match(names(ans), nn)
  ii <- ii[!is.na(ii)]
  ans <- ans[nn[ii]]
  ans
}

get_xlevels.merMod <- function(obj){
  ##cat("get_xlevels.lmerMod\n")
  ans <- lapply(obj@frame,levels)
  ans <- ans[unlist(lapply(ans, length))>0]

  nn <-attr(terms(obj), "term.labels")
  ii <- match(names(ans), nn)
  ii <- ii[!is.na(ii)]
  ans <- ans[nn[ii]]
  ##print(ans)
  ans
}



get_contrasts <- function(obj){
  UseMethod("get_contrasts")
}

get_contrasts.default <- function(obj){
  obj$contrasts ## FIXME: check RL code
}

get_contrasts.merMod <- function(obj){
  attr(model.matrix(obj), "contrasts")
}

set_xlevels <- function(xlev, at){
  nam     <- names(xlev)
  nn      <- match(nam, names(at))
  nn      <- nn[!is.na(nn)]
  at.fact <- at[nn]
  xlev[names(at[nn])]   <- at.fact
  attr(xlev, "at.fact") <- at.fact
  xlev
}


get_vartypes <- function(object){
    #cat("get_vartypes:\n")
    ## Common!!
    trms <- terms( model.frame( object ))
    trms <- delete.response( trms )
    att  <- attributes( trms )
    #cat("factors:\n"); print(att$factors)
    if (length(att$factors)>0){
        rhs.terms <- rownames(att$factors)[rowSums(att$factors)>0]
        rhs.class <- att$dataClass[match(rhs.terms, names(att$dataClass))]
        nums      <- rhs.terms[rhs.class=="numeric"]
        ##print(list(rhs.terms=rhs.terms, rhs.class=rhs.class, nums=nums))
        fact      <- rhs.terms[rhs.class=="factor"]
    } else {
        nums <- character(0)
        rhs.class <- character(0)
        fact <- character(0)
    }

    ## END

    ## print("here")
    ## print(list(numeric=nums, factor=fact))
    list(numeric=nums, factor=fact)
}

set_covariate_val <- function(xlev, covariateVal){
  nam     <- names(xlev)
  nn      <- match(nam, names(covariateVal))
  nn      <- nn[!is.na(nn)]
  con.con <- covariateVal[nn]
  xlev[names(covariateVal[nn])] <- con.con
  xlev
}


.get_covariate_ave <- function(object, at=NULL, tt=terms(object)){
    ##cat(".get_covariate_ave:\n")
    ## Common!!
    trms <- terms( model.frame( object ))
    trms <- delete.response( trms )
    att  <- attributes( trms )
    ##ff <<- att$factors
    if (length(att$factors)>0){
        rhs.terms <- rownames(att$factors)[rowSums(att$factors)>0]
        rhs.class <- att$dataClass[match(rhs.terms, names(att$dataClass))]
        nums      <- rhs.terms[rhs.class=="numeric"]
    } else {
        nums <- character(0)
        rhs.class <- character(0)
    }
    ##print(nums)
    ## END

    ## calculate average:
    ans  <- lapply(model.frame(object)[,nums, drop=FALSE], mean)
    ##cat("ans:\n"); print(ans)
    nn   <- match(names(ans), names(at))
    nn   <- nn[!is.na(nn)]
    at.num <- at[nn]
    ans[names(at[nn])]  <- at.num
    ##cat("ans (after insertion of 'at' values):\n"); print(ans)
    attr(ans, "at.num") <- at.num
    ans
}



get_X <- function(object, newdata, at=NULL){
  UseMethod("get_X")
}

get_X.default <- function(object, newdata, at=NULL){
    ##cat("get_X.default\n"); print(newdata)
    tt <- terms(object)

    Terms  <- delete.response(tt)
    #' print("HHHHHHHHH")
    #' print(Terms)
    #' print(newdata)

    #' xlev <<- get_xlevels(object)

    offval <- NULL
    #' print(newdata)
    M <- attr(terms(formula(object)),"factors")
    idx<-grep("offset\\(", rownames(M))
    if (length(idx)>0){
        off<- rownames(M)[idx]
        offvar<-gsub("offset\\((.*)\\)", "\\1\\", off)
        if ( !is.null(at) && !is.null(at[offvar]) )
            offval <- at[offvar]
        else
            offval <- 1
        d <- data.frame(rep(offval, nrow(newdata)))
        names(d) <- offvar
        newdata <- cbind(d, newdata)
    }

##    print(newdata)
    mf  <- model.frame(Terms, newdata, xlev = get_xlevels(object))

    X   <- model.matrix(Terms, mf, contrasts.arg = get_contrasts(object))
    attr(X, "assign" )    <- NULL
    attr(X, "contrasts" ) <- NULL
    attr(X, "offset" )    <- offval
    X
}


get_X.merMod <- function(object, newdata, at=NULL){
   tt <- terms(object)
   Terms  <- delete.response(tt)
   mf  <- model.frame(Terms, newdata, xlev = get_xlevels(object))
   X   <- model.matrix(Terms, mf, contrasts.arg = get_contrasts(object))
#  X <- getME(object, "X")
   attr(X,"assign")<-NULL
   attr(X, "contrasts") <- NULL
   X
}







