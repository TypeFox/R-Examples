#######################################################################
##
## Function: anchors.data() 
## Author  : Jonathan Wand <wand(at)stanford.edu>
##           http://wand.stanford.edu
## Created :  2008-04-20
##
## DESCRIPTION: Function for creating data used by all anchors objects
##              Replaces separate data processing in old functions chopit, anchors, vignette.order, etc
## 
## NOTES:       Primarily used within anchors().
##              Not normally invoked directly by users, but it is available in public scope.
## 
## OUTPUT:      object of class 'anchors.data'
##              containing a list of named vectors/matrices
##
## INPUT:
##   formula  : a list of formula, e.g., list(self=...,vign=...,tau=...,)
##              or a a single formula, e.g., self ~ vign1 + ... + vignJ 
##   data     : matrix or dataframe 
##   method   : type of analysis, choose one
##   subset   : (optional) logical statement based on variables in dataset
##   na.action: usual usage
##   delete   : "minimal" deletes only cases with missing values that affect component of model
##              "listwise" forces listwise deletion on the basis
##              of ALL variables in the formula list EVEN if not vars not used by method
##   debug    : 0=none, 1=extra printing
##
## Modified:  2008-04-29 : JW
##         
#######################################################################

anchors.data <- function(formula,
                         data,
                         method,
                         subset,
                         na.action   = na.omit,
                         na.response = c(NA, 0), 
                         min.response = 1,
                         delete = c("minimal","maximal"), debug=0) {

  delete <- match.arg(delete)

  ## 1. Standardize formula list:
  if (class(formula) == "formula") {
    ## break it apart, LHS is self, RHS are vign
    f1  <- sapply( formula, deparse,  width.cutoff = 500)
    if (length(f1) == 3) {
      f2 <- unlist( strsplit( f1[[3]], "+", fixed=TRUE))
      ## and remake formula in standard form
      formula <- list(self = as.formula(paste( f1[[2]] , "~1", sep="")),
                      vign = as.formula(paste("cbind(",paste(f2, collapse=","),") ~ 1")),
                      tau  = ~ 1)
    } else if (length(f1) == 2) {
      f2 <- unlist( strsplit( f1[[2]], "+", fixed=TRUE))
      ## and remake formula in standard form
      formula <- list(vign = as.formula(paste("cbind(",paste(f2, collapse=","),") ~ 1")),
                      tau  = ~ 1)
    }
  } else if (class(formula) != "list") {
    stop("formula arg must either be a formula or a list of formulas\n")
  }

  ## 2. Extract variable information
  nf <- names(formula)
  no.self <- !("self"   %in% nf)
  no.vign <- !("vign"   %in% nf)
  no.tau  <- !("tau"    %in% nf)
  no.tau1 <- !("tau1"   %in% nf)
  no.cpolr<- !("cpolr" %in% nf)

  if (debug > 0) {
    cat("anchors.data: Missing",
        "self" , no.self,
        "vign" , no.vign,
        "tau"  , no.tau,
        "tau1" , no.tau1,
        "cpolr", no.cpolr,"\n")
  }
  
  cfvignx <- cfself <- cfvign <- cftau <- cftau1 <- cfcpolr <- NULL 
  if      (!no.self) { fself <- as.formula(formula$self ); cfself <- as.character(fself ) } 
  if      (!no.vign) { fvign <- as.formula(formula$vign ); cfvign <- as.character(fvign ) } 
  if      (!no.tau ) { ftau  <- as.formula(formula$tau  ); cftau  <- as.character(ftau  ) } 
  if      (!no.tau1) { ftau1 <- as.formula(formula$tau1 ); cftau1 <- as.character(ftau1 ) } 
  else if (!no.tau ) { ftau1 <- as.formula(formula$tau  ); cftau1 <- as.character(ftau  ) }

  if (no.cpolr  ) { formula$cpolr <- ~ 1; no.cpolr <- FALSE }
  if (!no.cpolr ) { fcpolr<- as.formula(formula$cpolr); cfcpolr<- as.character(fcpolr) } 

  ## coding errors
  if  (!no.self && length(fself) != 3) 
    stop("\n'self=' formula is improperly formed: should be of form ' self = y1 ~ x1 + x2 '\n")
  if (!no.vign && length(fvign) != 3)
    stop("\n'vign=' formula is improperly formed: should be of form ' vign = cbind(z0,z2,z3) ~ 1 '\n")
  if (!no.tau  && length(ftau) != 2)
    stop("\n'tau =' formula is improperly formed: should be of form ' tau  =  ~ x1 + x2 '\n")
  if (!no.tau1 && length(ftau1) != 2)
    stop("\n'tau1=' formula is improperly formed: should be of form ' tau1  =  ~ x1 + x2 '\n")
  if (!no.cpolr  && length(fcpolr) != 2)
    stop("\n'cpolr =' formula is improperly formed: should be of form ' cpolr  =  ~ x1 + x2 '\n")
  ## logic errors
  if (no.self && method %in% c("rank","entropy"))
    stop(paste("must specify a 'self' forumula to use 'anchors(...,method='",method,"')'\n"))
  if (no.cpolr && method == "cpolr")
    stop("must specify a 'cpolr' formula to use 'anchors(...,method='cpolr')'\n")
  if (no.tau   && method == "chopit")
    stop("must specify variables in 'tau' formula to use 'anchors(...,method='chopit')'\n")

  if (debug > 0) cat("anchors.data: Passed all coding tests\n")
  
  if (!no.vign) {
    if (is.na(charmatch("cbind(",cfvign[2]))) {
      cfvignx <- cfvign[2]
    } else {
      cfvignx <- sapply( fvign[[2]], deparse, width.cutoff = 500)[-1]
    }
  }

  ## reformulate formulas to make subsetting easier
  null.omit <- function(x) unique(x[!is.null(x)])

  tmp <- null.omit( c(cfself[3] , cftau1[2], cftau[ 2], "1") )
  fself.chopit <- as.formula(paste(cfself[2] ,"~", paste(tmp , collapse="+" )))

  tmp <- null.omit( c(cftau1[2], cftau[ 2], "1") )
  fvign.chopit <- as.formula(paste("~", paste(tmp , collapse="+" )))

  tmp <- null.omit( c("-1", cfvignx ) )
  fself.rank   <- as.formula(paste(cfself[2] ,"~",paste(tmp,collapse="+")))

  tmp <- null.omit( c(cfvignx, cfcpolr[2] ) )
  fself.cpolr   <- as.formula(paste(cfself[2] ,"~",paste(tmp,collapse="+")))

  fvign.order  <- as.formula(paste(cfvign[2],"~ 1"))

  tmp <- null.omit( c(cfvignx, cfself[3],cftau1[2], cftau[ 2], cfcpolr[2], "1") )
  fall <- as.formula(paste(cfself[2], "~", paste(tmp , collapse="+" )))

  formula2 <- list(fself.chopit=fself.chopit,
                   fvign.chopit=fvign.chopit,
                   fself.rank  =fself.rank  ,
                   fvign.order =fvign.order ,
                   fself.cpolr =fself.cpolr , 
                   fall        =fall        )

  
  if (debug>0) {
    cat("formula\n")
    print(formula2)
  }
  
  ## 3. Build subset data
  ## a) keep maximum data (i.e., use only parts of formula needed for current method)
  ## b) maximal delete across components of chopit (complete rows only) even if method != chopit


  ##
  if (!missing(subset))  {
    e <- substitute(subset)
    r <- eval(e, data, parent.frame())
    if (!is.logical(r)) 
      stop("'subset' must evaluate to logical")
    r <- r & !is.na(r)
    data <- data[r, TRUE, drop = TRUE]
  }

  if (debug > 0) cat("anchors.data: passed subset\n")

  
  if (delete == "maximal")  {

    ## for simplicity, make any missing response value NA
    if (any(!is.na(na.response))) {
      ii <- na.response[!is.na(na.response)]
      for (i in ii)
        data <- replace.value( data, c( cfself[2], cfvignx ), from=i, to=NA,
                              verbose = (debug > 0) )
    }
    ## this is all vars included in model frame
    data <- model.frame( fall, data=data, na.action = na.action)
  }

  if (debug > 0) cat("anchors.data: Passed maximal\n")
  
  y0 <- z0 <- x0 <- NULL
  v0s <- v0s1 <- v0v <- v0v1 <- NULL
  if (method == "chopit") {

    yid <- NULL
    if (!no.self) {

      data.tmp <- model.frame( fself.chopit, data=data, na.action=na.action)
      v0s <- as.matrix(model.matrix(ftau , data.tmp))
      v0s1<- as.matrix(model.matrix(ftau1, data.tmp))
      keep.v0s <- rownames(data.tmp)
      ## in case we have intercept only
      rownames(v0s) <- rownames(v0s1) <- rownames(data.tmp)
      
      ## NOTE: x0 and y0 doesn't get trimmed of cases with missing tau vars:
      data.tmp <- model.frame( fself, data=data, na.action=na.action)
      x0 <- as.matrix(model.matrix(fself, data.tmp))
      y0 <- as.matrix(model.response(data.tmp, type="numeric"))
      ## in case we have intercept only
      rownames(x0) <- rownames(data.tmp)

      ## now subset:
      idx <- match( rownames(y0), keep.v0s )
      idx <- idx[!is.na(idx)]
      x0 <- as.matrix(x0[idx,])
      y0 <- as.matrix(y0[idx,])



      if (NROW(v0s) != NROW(y0))
        stop("ancors.data: mismatch in number of cases in y0 and v0s")
      
      drop.y0 <- which(apply(y0, 1, function(x) { any(x %in% na.response) } ))
      if (length(drop.y0) > 0) {
        v0s <-  v0s[-drop.y0,]
        v0s1<- v0s1[-drop.y0,]
        y0 <- as.matrix(y0[-drop.y0,])
        x0 <- as.matrix(x0[-drop.y0,])
      }
      
      yid <- sort(unique(y0))
      min.yid <- min(yid)
      if (min.yid < min.response)
        stop(paste("There are self response values less than min.response:",
                   paste( yid[yid < min.response ] , collapse=","),
                   "\n"))
      #if (min(yid) < 1)  y0 <- y0 - (min(yid) - 1)
    }

    data.tmp <- model.frame( fvign.chopit, data=data, na.action=na.action)
    v0v <- as.matrix(model.matrix(ftau , data.tmp))
    v0v1<- as.matrix(model.matrix(ftau1, data.tmp))
    ## in case we have intercept only
    rownames(v0v1) <- rownames(v0v) <- rownames(data.tmp)
    keep.v0v <- rownames(data.tmp)
    
    if (debug>2) { cat("anchors.data: v0v(1):\n"); print(keep.v0v) }
#    print(keep.v0v)
    
    ## OL: provided code for treating NA instead of zeros in z0
    ##     but it only worked in limited cases
    ## JW: corrects OL code to handle single vignette case
    ##     and for no covariates (intercept only)
    ##     also generalizes to arbitrary missing values in na.response
    data.tmp <- data[ match( keep.v0v, rownames(data) ) , ]
    z0 <- as.matrix(model.response(model.frame(fvign, data=data.tmp,na.action=NULL),"numeric"))
    keep.z0 <- apply(z0, 1, function(x) { !all(is.na(x) |  x %in% na.response) } )
    
#    cat("xx\n")
#    print( unlist(drop.z0))
#    drop.z0 <- as.matrix(drop.z0[ drop.z0])
#    print(drop.z0)
#    cat("z0\n")
#    print(rownames(drop.z0))
##    print(rownames(v0v))
    
    if ( sum(!keep.z0) > 0) {
      keep.z0 <- names(keep.z0[keep.z0])
      v0v <- as.matrix( v0v[match(  keep.z0 , keep.v0v ),])
      v0v1<- as.matrix(v0v1[match(  keep.z0 , keep.v0v ),])
      z0  <- as.matrix(  z0[match(  keep.z0 , rownames(z0  ) ),])
      if (debug>2) { cat("anchors.data: v0v(2):\n"); print(v0v) }
    }

    zid <- sort(unique(z0))
    zid <- zid[ zid != 0 ]
    min.zid <- min(zid)
    if (min.zid < min.response)
        stop(paste("There are vign response values less than min.response and non-zero:",
                   paste( zid[zid < min.response ] , collapse=","),
                   "\n"))
    
#    if (min(yid) < 1)  z0 <- z0 - (min(yid) - 1)
    
    z0[is.na(z0)] <- 0

    if (debug > 0) cat("anchors.data: finished chopit data\n")
    
  }

  if (method %in% c("rank","entropy")) {

    data.tmp <- model.frame( fself.cpolr, data=data, na.action=na.action)
    y0 <- as.matrix(model.response(data.tmp, "numeric"))
    z0 <- as.matrix(model.matrix( fself.rank, data.tmp))
    x0 <- as.matrix(model.matrix(fcpolr    , data.tmp))

#     data.tmp <- model.frame( fself.rank, data=data, na.action=na.action)
#     y0 <- as.matrix(model.response(data.tmp, "numeric"))
#     z0 <- as.matrix(model.matrix( fself.rank, data.tmp))
    
    drop <- which( apply( cbind(y0,z0), 1, function(x) { any(x %in% na.response) }))
    if (length(drop)>0) {
      y0 <- as.matrix(y0[-drop,])
      z0 <- as.matrix(z0[-drop,])
      x0 <- as.matrix(x0[-drop,])
    }
    if (debug > 0) cat("anchors.data: finished rank data\n")
  }

  
  if (method == "order") {
    data.tmp <- model.frame( fvign.order, data=data, na.action=na.action)
    z0 <- as.matrix(model.response(data.tmp, "numeric"))
    if (NCOL(z0) == 1) stop("Doesn't make sense to look at order of single vignette\n")

    drop <- which( apply( z0, 1, function(x) { any(x %in% na.response) }) )
    if (length(drop)>1) {
      z0 <- as.matrix(z0[-drop,])
#      print(dim(z0))
    }
    if (debug > 0) cat("anchors.data: finished order data\n")
  }

#   cat("test\n")
#   print( cfvign[2] )
#   print( dim(z0) )
#   print( colnames(z0) )
  
  ## JW: does making a matrix matter for any of these?
  if (debug > 0) cat("anchors.data: begin naming\n")
  if(!is.null(x0   ) & NCOL(x0  ) == 1 & is.null(colnames(x0  ))) {  colnames(x0)  <- cfself[3] }
  if(!is.null(y0   ) & NCOL(y0  ) == 1 & is.null(colnames(y0  ))) {  colnames(y0)  <- cfself[2] }
  if(!is.null(v0v  ) & NCOL(v0v ) == 1 & is.null(colnames(v0v ))) {  colnames(v0v) <- cftau[2]  }
  if(!is.null(v0v1 ) & NCOL(v0v1) == 1 & is.null(colnames(v0v1))) {  colnames(v0v1)<- cftau1[2] }
  if(!is.null(z0   ) & NCOL(z0  ) == 1 & is.null(colnames(z0  ))) {  colnames(z0)  <- cfvign[2] }
  if (debug > 0) cat("anchors.data: finished naming\n")
  
  ## quality checks:
#  if (!noself && (count$nobs.y0 != count$nobs.x0 || count$nobs.y0 != count$nobs.v0s)) {
#    stop(paste("Error: Number of cases in components of self-likelihood differ\n",
#               "nobs.y0:",count$nobs.y0,"nobs.x0:",count$nobs.x0,"nobs.v0s:",count$nobs.v0s,"\n"))
#  }
#  if (count$nobs.z0 != count$nobs.v0v ) {
#    stop(paste("Error: Number of cases in components of vign-likelihood differ\n",
#               "nobs.z0:",count$nobs.z0,"nobs.v0v:",count$nobs.v0v,"\n"))
#  }
  

  ff <- function(x) {
    if (is.null(x))
      x
    else
      as.matrix(x)
  }
  
  data <- list(y0  = ff(y0  ),
               z0  = ff(z0  ),
               x0  = ff(x0  ),
               v0v = ff(v0v ),
               v0s = ff(v0s ),
               v0v1= ff(v0v1),
               v0s1= ff(v0s1),
               formula=formula,
               method=method,
               delete=delete
               )
  
  class(data) <- "anchors.data"
  
  return(data)
  
}
