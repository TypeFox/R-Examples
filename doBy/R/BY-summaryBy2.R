
summaryBy <-
  function (formula, data=parent.frame(), id=NULL, FUN=mean, keep.names=FALSE,
            p2d=FALSE, order=TRUE, full.dimension=FALSE,
            var.names=NULL, fun.names=NULL,
            ...)
  {
    debug.info <- 0

    zzz <- .get_variables(formula, data, id, debug.info) ## ; str(zzz)
    lhs.num <- zzz$lhs.num
    rhs.grp <- zzz$rhs.grp
    ids.var <- zzz$form.ids.var

    rh.trivial <- length( rhs.grp ) == 0 #; cat(sprintf("rh.trivial=%d\n", rh.trivial))

    rh.string <- .get_rhs_string( data, rhs.grp )
    rh.unique <- unique(rh.string)
    rh.idx    <- match(rh.unique, rh.string)
    rh.string.factor <- factor(rh.string, levels=rh.unique) ## This is important

### Get data for id.vars; use ids.var, data, rh.idx
    ## print(ids.var)
    ## .s <<- ids.var
    ## print(rh.idx)
    ## print(names(data))
    if (length(ids.var)>0){
      id.data <-  data[ rh.idx, ids.var, drop=FALSE ] ##; print(id.data)
    }

### Get lhs data; use lhs.num, data

    lh.data <- do.call(cbind,lapply(paste(lhs.num), function(x)eval(parse(text=x), data)))
    #print(lh.data)
    colnames(lh.data) <- lhs.num


### Function names; use FUN
    funNames <- .get_fun_names( FUN )

### Calculate groupwise statistics
    if (!is.list(FUN))
        FUN <- list(FUN)
    ans <- NULL
    for (ff in 1:length(FUN)) {  ## loop over functions
        ##currFUN <- FUN[[ff]]
        currFUN <- match.fun( FUN[[ff]] )
        for (vv in 1:length(lhs.num)) {  ## loop over variables
            currVAR <- lh.data[,lhs.num[vv]]
            zzz <- tapply(currVAR, rh.string.factor,
                          function(x){ currFUN(x,...) }, simplify=FALSE)
            zzz  <- do.call(rbind, zzz)
            ans  <- cbind(ans, zzz)
        }
    }



    if (!is.null(var.names) && length(var.names)==length(lhs.num))
      lhs.names <- var.names
    else
      lhs.names <- lhs.num

### Set names for columns

    #print(funNames)
    if (!is.null(fun.names) ) ##&& length(fun.names)==length(funNames))
      funNames <- fun.names

    newnames <- .get_col_names(ncol(ans), colnames(ans), funNames, lhs.names, keep.names)
    ##cat(sprintf("newnames    = %s\n", toString( newnames )))
    colnames(ans) <- newnames
    ans <- as.data.frame(ans)


### Pad the rhs data to the result
    if ( !rh.trivial ){
      ans <- cbind(data[rh.idx, rhs.grp, drop=FALSE], ans)
    }

### Pad id.data to result
    ##print(id.data)
    if (length(ids.var)>0){
      ans <- cbind(ans, id.data)
    }

### Must the result have full dimension?
    if (full.dimension){
      rrr <-as.numeric(rh.string.factor)
      ans <- ans[rrr,,drop=FALSE]
    }

### Order the result by the rhs
    if (order==TRUE && !rh.trivial){
      rhs.string  <- paste (rhs.grp, collapse='+')
      ans <- orderBy(as.formula(paste("~", rhs.string)), data=ans)
    }

### Replace '('s and ')'s with '.'s
    if (p2d)
      names(ans) <-  gsub("\\)","\\.", gsub("\\(","\\.",names(ans)))

### Finalize
    rownames(ans) <- 1:nrow(ans)
    if (length(unique(names(ans))) != length(names(ans)))
      warning("dataframe contains replicate names \n", call.=FALSE)

    ans
  }





.get_rhs_string <- function(data, rhs.var, sep.string="@"){

  if (length(rhs.var)==0){
    rep.int("1", nrow(data))
  } else {
    rh.string <- paste(data[,rhs.var[1]])
    if (length( rhs.var ) > 1){
      for (ii in 2:length( rhs.var )){
        rh.string <- paste(rh.string, sep.string, data[, rhs.var[ii]], sep='')
      }
    }
    rh.string
  }
}


.get_fun_names <- function( FUN ){

    if (!is.list(FUN))
      funNames <- paste(deparse(substitute(FUN, env=parent.frame())), collapse = " ")
    else
      funNames <- unlist(lapply(substitute(FUN, env=parent.frame())[-1], function(a) paste(a)))

    ##cat(sprintf("funNames   = %s\n", toString(funNames)))
    funNames
}




.get_col_names <- function(ncol.ans, colNames, funNames, lhs.num, keep.names){
### Names for new variables

  ## Does the columns of ans have names??
  oldnames <- colNames #colnames(ans)
  if (is.null(oldnames))
    hasNames <- 0
  else {
    hasNames <- 1*(prod(nchar(oldnames))>0)
  }
  ##cat(sprintf("hasNames=%i\n", hasNames))

  ## Dim of response (per variable on LHS)
  dimr     <- (ncol.ans)/length(lhs.num)

  only.one.response <- (ncol.ans==length(lhs.num))
  if ( keep.names && only.one.response ){
    newnames <- lhs.num
  } else {
    if (hasNames>0){ ## newnames = var.fun
      funNames <- colNames[1:dimr]
      newnames <- unlist(lapply(lhs.num, function(v){paste(v, funNames, sep='.')}))
    } else {
      if (length(funNames) != dimr){
        funNames <- paste("FUN", 1:dimr, sep='')
        newnames <- unlist(lapply(lhs.num, function(v){paste(v, funNames, sep='.')}))
      } else {
        newnames <- unlist(lapply(funNames, function(x) paste(lhs.num, x, sep='.')))
      }
      if (length(newnames)!=ncol.ans){
        funNames <- paste(funNames, 1:dimr, sep=".")
        newnames  <- unlist(lapply(funNames, function(x) paste(lhs.num, x, sep='.')))
      }
    }
  }
  newnames
}


.get_variables <- function(formula, data, id, debug.info){

  data.var  <- names(data)

  if ( !(class(formula) %in% c("formula", "list")) ){
    stop("'formula' must be a formula or a list")
  }

  if (class(formula) %in% "formula"){
    if (length(formula) != 3) stop("Formula must have a left hand side")
    rhs       <- formula[[3]]
    form.rhs.var   <- all.vars(rhs) ## May contain "." and "1"

    lhs       <- formula[[2]]
    form.lhs.var   <- all.vars(lhs) ## May contain "."
    #print(form.lhs.var)

    zz <- .lhsParse(lhs)
    form.lhs.var <-
      if (length(zz)==1)
        paste(zz)
      else
        paste(unlist(.lhsParse(lhs)))

    ##print(form.lhs.var)

  } else {
    if (length(formula)>=2){
      lhs <- formula[[1]]
      rhs <- formula[[2]]
      form.lhs.var <- lhs
      form.rhs.var <- rhs
    } else {
      stop("Invalid specification of formula")
    }
  }

  if (is.null(id)){
    form.ids.var <- character(0)
  } else {
    if ( !(class(id) %in% c("formula", "character")) ){
      stop("'id' must be a formula or a character vector")
    }
    if (class(id) %in% "formula"){
      form.ids.var   <- all.vars(id)
    } else {
      form.ids.var   <- id
    }
  }

  data.cls  <- lapply(data, class)
  data.num.idx   <- data.cls %in% c("numeric","integer")
  data.num.var   <- data.var[ data.num.idx  ]
  data.fac.var   <- data.var[ !data.num.idx ]

##   print(form.lhs.var)
##   print(data.num.var)
  lhs.num   <- intersect( form.lhs.var, data.num.var )
  rhs.num   <- intersect( form.rhs.var, data.num.var )
  ids.num   <- intersect( form.ids.var, data.num.var )
  lhs.fac   <- intersect( form.lhs.var, data.fac.var )
  rhs.fac   <- intersect( form.rhs.var, data.fac.var )
  ids.fac   <- intersect( form.ids.var, data.fac.var )

  lll <- list(data.var=data.var,
              form.lhs.var=form.lhs.var, form.rhs.var=form.rhs.var, form.ids.var=form.ids.var,
              lhs.num=lhs.num, rhs.num=rhs.num, ids.num=ids.num,
              lhs.fac=lhs.fac, rhs.fac=rhs.fac, ids.fac=ids.fac )

  #if (debug.info>=1)
  ##{ cat("status:\n"); str(lll, vec.len=20) }

  if ( "." %in% form.lhs.var ){ ## need all numeric variables not metioned elswhere on lhs
    form.lhs.var <- setdiff(form.lhs.var, ".")
    lhs.num <- union( form.lhs.var, setdiff(data.num.var, c(rhs.num, ids.num)))
    if ( length( lhs.fac ) > 0 ){
      isSpecial <- rep(NA, length( lhs.fac ))
      for (j in 1:length(lhs.fac)){
        isSpecial[j]<- (class(data[,lhs.fac[j]])[1] %in% c("POSIXt", "Date"))
      }
      lhs.num <- union( lhs.num, lhs.fac[ isSpecial ] )
    }
  } else {
    lhs.num <- form.lhs.var
  }

  ## The grouping variable
  if ("." %in% form.rhs.var){ ## need all factors not mentioned elsewhere as grouping factors
    free.fac <- setdiff(data.fac.var, c(lhs.fac, ids.fac))
    rhs.grp  <- c(setdiff(form.rhs.var, "."), free.fac)
  } else {
    rhs.grp <- form.rhs.var
  }
  rhs.grp <- intersect( rhs.grp, data.var )

  rrr <- list(lhs.num=lhs.num, rhs.fac=rhs.fac, form.ids.var=form.ids.var,
            form.rhs.var=form.rhs.var, rhs.grp=rhs.grp)
  ##str(rrr)
  rrr
 }



.lhsParse <- function(x){
  ##cat(".lhsParse:"); print(x); print(class(x))
  if (class(x)=='name'){
    value <- x
  } else {
    s <- paste(x[[1]])
    value <- switch(s,
                    '+'={  c(.lhsParse(x[[2]]),.lhsParse(x[[3]]))},
                    'I'={  x[[2]]},
                    {  deparse(x)})
  }
  return(value)
}
