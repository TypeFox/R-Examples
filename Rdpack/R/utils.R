.tag_in <- function(x, tags){
    a <- attr(x,"Rd_tag")
    is.character(a) && length(a) == 1 && a %in% tags
}

.tag_eq <- function(x, tag){
    a <- attr(x,"Rd_tag")
    is.character(a) && length(a) == 1 && a == tag
}

.ocompare <- function(x,y){  # outer compare
    f <- function(z) sapply(x,  function(w) identical(z,w))
    sapply(y,f)
}

.asym_compare <- function(x,y){
    cmp <- .ocompare(x,y)
    if(length(x) == 1){                      # 2012-10-03 added this if
        cmp <- matrix(cmp, nrow=1)
    }
    if(length(x) == 0){
        inew <- seq_along(y)
        irem <- icom <- integer(0)
    }else if(length(y) == 0){
        irem <- seq_along(x)
        inew <- icom <- integer(0)
    }else{
        inew <- which( ! apply(cmp, 2, any) )               # in y but not in x
        irem <- which( ! apply(cmp, 1, any) )               # in x but not in y
        icom <- which( !(seq_along(x) %in% c(inew, irem)) ) # in both, x and y
    }

    list( i_new     = inew,
          i_removed = irem,
          i_common  = icom )
}


rattr <- function(x,y){                                # restore attributes, see Rdapply below
    if(is.list(x)){
        if(!identical(attributes(x), attributes(y)))
            attributes(x) <- attributes(y)

        for(i in seq_along(x))
            x[[i]] <- Recall(x[[i]], y[[i]])
    }
    x
}

Rdapply <- function(x, ...){                       # needed since rapply loses attributes of x
    wrk <- rapply(x, ...)
    rattr(wrk, x)
}

.pos_in_eqnpos <- function(verb, eqn){  # not a complete solution
    f <- function(x,y) length(x) > length(y)  && all( x[seq_along(y)] == y )

    flag  <- FALSE
    res <- vector(mode="logical", length = length(verb))
    for(i in seq_along(verb)){
        for(j in seq_along(eqn)){
            flag <- f(verb[[i]], eqn[[j]])
            if(flag)
                break
        }
        res[i] <- flag
    }

    res
}
                                                                             # 2012-09-08 new
Rdtagapply <- function(object, FUN, rdtag, classes = "character", how = "replace", ...){
    if(rdtag %in% c("mathVERB", "nonmathVERB", "nonmath")){  # special care for \eqn and \deqn
        feqnpos <- function(x) .tag_in(x, c("\\eqn", "\\deqn"))
        eqnpos <- Rdo_locate(object, feqnpos, lists = TRUE)
        eqnpos <- lapply(eqnpos, function(x) c(x,1))  # 1st arg. of \eqn or \deqn

        if(rdtag == "nonmath"){
            anypos <- Rdo_locate_leaves(object)
            pos <- anypos[ !.pos_in_eqnpos(anypos, eqnpos) ]
        }else{                                   # is it better to use Rdo_locate_leaves here?
            fverbpos <- function(x) .tag_eq(x, "VERB")
            verbpos <- Rdo_locate(object, fverbpos, lists = FALSE)

            pos <- if(rdtag == "mathVERB")                # verbpos[ verbpos %in% eqnpos ]
                       verbpos[ .pos_in_eqnpos(verbpos, eqnpos) ]
                   else # (rdtag == "nonmathVERB")        # verbpos[ !( verbpos %in% eqnpos) ]
                       verbpos[ !.pos_in_eqnpos(verbpos, eqnpos) ]
        }

        for(ind in pos)
            object[[ind]] <- FUN(object[[ind]], ...)

        object
    }else{
        locfun <- function(x,...) if(.tag_eq(x, rdtag))  FUN(x,...) else x
        Rdapply(object,locfun, classes = classes, how = how, ...)
    }
}


.aux_Rdo_locate_leaves <- function(x, i){
    if(is.logical(x)){
        if(isTRUE(x))
            i
        else
            FALSE
    }else{     # here x is a list of integer vectors
        lapply(x, function(x){ c(i,x)})
    }
}

                              # inn1 <- parse_Rd("inner-methods.Rd")
                              # Rdo_locate_leaves(inn1, function(x) grepl("^signature\\(", x))
Rdo_locate_leaves <- function(object, f = function(x) TRUE){
    fxx <- function(x){
        if(is.character(x)){
            return(f(x))         # evaluate `f' for this leaf
        }else if(is.list(x)){
            if(length(x)==0)   # list()
                return(x)
            wrk <- lapply(x, fxx)
            for(i in seq_along(wrk))
                wrk[[i]] <- .aux_Rdo_locate_leaves(wrk[[i]], i)
            indx <- sapply(wrk, function(y) !identical(y,FALSE) )
            wrk <- wrk[indx]
            wrk <- do.call("c",wrk)
            return(wrk)
        }else{
            return(FALSE)     # todo: replace this with a function and additional argument
        }
    }

    fxx(object)
}


.merge_pos0 <- function(a, b){                                                # 2012-09-25 new
    if(is.numeric(b))
        c(a,b)
    else if(is.list(b)){
        res <- vector(mode="list", length=length(b))
        for(i in seq_along(b)){
            res[[i]] <- Recall(a, b[[i]])
        }
        res
    }else if(is.logical(b)){
        if(isTRUE(b))
            a
        else
            as.numeric(NA) # numeric(0)
    }else{
        attr(a, "merge_tag") <- b # obache tozi tag se gubi po-kasno
        a
    }
}

.merge_pos <- function(a, b, ...){                                            # 2012-09-25 new
    if(is.numeric(a) && length(a)==1){
        a <- as.list(1:a)  # todo: process the case a=0
    }

    wrk <- list()
    for(i in a){
        wrk <- c(wrk, .merge_pos0(a[[i]], b[[i]]) )
    }
    res <- wrk
    len <- sapply(res, length)

    res <- res[len>0]   # drop zero length values

    drop <- sapply(res, function(x) any(is.na(x)))
    res[ !drop]
}

.notFALSE <- function(x) !identical(x, FALSE)
.notTRUE  <- function(x) !isTRUE(x)  # !identical(x, TRUE)

                                                         # 2012-10-20 dobavyam argument nested
                                                         # 2012-10-07 dobavyam argument fpos
Rdo_locate <- function(object, f = function(x) TRUE, pos_only = TRUE, lists = FALSE,
                       fpos = NULL, nested = TRUE){
    fxx <- function(x){
        if(is.character(x)){   # a leaf, evaluate `f' for it
            f(x)
        }else if(is.list(x)){
            fli <- if(lists) f(x)
                   else      FALSE

            if(!nested && fli)
                return(TRUE)

            wrk <- if(fli) list(fli)
                   else    list()

            if(length(x)>0)
                wrk <- c(wrk, .merge_pos(length(x), lapply(x, fxx)) )   # recurse
            wrk
        }else{                # this should not happen for Rd objects
            FALSE             # todo: replace this with a function and additional argument?
        }
    }

    pos <- fxx(object)

    f <- if(is.function(fpos))
             function(x) list(pos = x, value = fpos(object, x))
         else if(isTRUE(pos_only))
             NULL
         else if(is.function(pos_only))
             function(x) list(pos = x, value = pos_only(object[[x]]))
         else
             function(x) list(pos = x, value = object[[x]])

    if(!is.null(f))
        lapply(pos, f)
    else
        pos
}

                                                 # insert sep between successive elements of x
.Rdinter <- function(x, sep = Rdo_newline(), before_first = FALSE, after_last = FALSE){
   indx <- c(rbind(seq_along(x), length(x) + 1))
   if(before_first)
       indx <- c( length(x) + 1, indx)
   if(!after_last)
       indx <- head(indx,-1)

   c(x, list(sep))[indx]
}

.Rd_tidy <- function(rdo){                           # todo: optionally remove empty sections?
    tags <- Rdo_tags(rdo)                              # todo: remove multiple empty lines?
    secpos <- which(tags != "TEXT")                   #       remove empty lines at top level?

    secpos <- secpos[secpos > 1]
    if(length(secpos)==0)
        return(rdo)


    nl <- Rdo_newline()
    nlpos <- integer(0)
    for(i in secpos){         # 2012-09-20 smenyam  (!identical(rdo[[i-1]], nl)) s dolnoto.
        if(!Rdo_is_newline(rdo[[i-1]]))   # inache tryabva da se grizha za "srcref" attribute
            nlpos <- c(nlpos, i)
    }

    if(length(nlpos)==0)
        return(rdo)
                                                                          # todo: this is lazy
    for(i in sort(nlpos, TRUE)){   # sort in decreasing order, otherwise indices change
        rdo <- Rdo_insert_element(rdo, nl, i)
    }
    rdo
}

Rdo_remove_srcref <- function(rdo){
    f <- function(x){
        a <- attributes(x)
        if(is.list(x)){
            x <- lapply(x,f)
            if(!is.null(a))        # 2012-10-07 check for NULL
                attributes(x) <- a
        }
        if(!is.null(attr(x,"srcref")))
            attr(x,"srcref") <- NULL
        x
    }
    f(rdo)
}

Rd_combo <- function(rd, f, ..., .MORE){               # todo: allow a directory?
    dots <- list(...)
    flag <- missing(.MORE)
    locfun <- function(x){
                  locwrk <- parse_Rd(x)
                  locargs <- c(list(locwrk), dots)
                  locres <- try( do.call(f, locargs), silent=TRUE)
                  if(flag || inherits(locres,"try-error"))
                      locres
                  else{
                      try( do.call(.MORE, list(locres)), silent=TRUE)
                  }
              }
    lapply(rd, locfun)
}
