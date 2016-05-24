## Turn a formula into a list of generators (glist)
##
.pFormula2 <- function (formula, varnames, marginal=NULL, interactions=NULL) ##, v.sep = ":", g.sep = "+", ignore.power.value=FALSE) 
{

  used.var <- if (length(marginal) > 0) { marginal } else { varnames }
 
  clformula <- class(formula) ##; cat("class(formula) :", clformula, "\n")
  
  switch(clformula,
         "formula"={
            pow <- .extract.power(formula)
            ##cat(sprintf("A formula is given; power=%d\n", pow))
            if (is.integer(pow)){
              if (identical(pow, -1L)){
                ##cat("The saturated model\n")
                glist <- list(used.var)
              } else {
                if (identical(pow, 1L)){
                  ##cat("The independence model\n")
                  glist <- as.list(used.var)
                } else {
                  pow   <- min(c(pow, length(used.var)))
                  glist <- combn(used.var, pow, simplify=FALSE)
                }               
              }              
            } else {
              ##cat("A proper formula\n")
              glist <- rhsFormula2list(formula)
            }
         },
         "list"={
           glist <- formula
         },
         "graphNEL"={
           glist <- maxCliqueMAT(as.adjMAT(formula))[[1]]
         },
         "matrix"={
           glist <- maxCliqueMAT(formula)[[1]]
         })

  glist <- .check.glist(glist, used.var)  
  if (!is.null(interactions))
    glist <- .set.interactions(glist, interactions)

  value <- list(glist    = glist,
                varNames = uniquePrim(unlist(glist)))
  return(value)
}


.check.glist <- function(glist, used.var){
  if (any(is.na(pmatch(unlistPrim(glist), used.var, 
                       duplicates.ok = TRUE)))) 
    stop("An invalid variable specification has been found\n")
  glist <- lapply(glist, function(x) {
    ii <- pmatch(x, used.var)
    used.var[ii]
  })
  modnames <- uniquePrim(unlist(glist))
  if(any(is.na(match(modnames, used.var))))
    stop("Variables in model not contained in the variable set. Perhaps a problem with 'marginal'?")
  
  glist
}


.set.interactions <- function(glist, interactions){
  zz <- lapply(glist, function(ss){
    if (length(ss)<=interactions){
      list(ss)
    } else {
      combn(ss, interactions, simplify=FALSE)
    }
  })
  zz <- removeRedundant(unlist(zz, recursive=FALSE))
  zz
}


.extract.power <- function (fff) 
{
    mimf <- paste(as.formula(fff))[2]
    mimf.split <- unlist(strsplit(mimf, ""))
    if (length(grep("[:alpha:]", mimf)) > 0) {
        pow <- mimf
    }
    else {
        has.hat <- match("^", mimf.split)
        sub <- unlist(strsplit(mimf, "\\^"))
        if (!is.na(has.hat)) {
            pow <- ifelse(sub[2] == ".", -1, as.numeric(sub[2]))
        }
        else {
            pow <- length(unlist(strsplit(sub, "\\.")))
        }
    }
    return(pow)
}


.extract.power <- function(fff){
  sss <- deparse(fff[[2]])
  has.hat <- length(grep("^\\.\\^", sss)) > 0
  if (!has.hat){
    return (NA)
  } else {
    rest <- gsub("\\.\\^","",sss)
    pp <- strsplit(rest, " ")[[1]][1]
    pow <- suppressWarnings(as.integer(pp))
    if (is.na(pow))
      pow <- -1L
    else
      if (identical(pow, 0L)){
        pow <- 1L
      }
    pow
  }   
}

  
## .glist2str.formula <- function(glist){
##   tmp <- lapply(glist, paste, collapse = v.sep)
##   str.formula <- paste(unlistPrim(tmp), collapse = g.sep, sep = "")
##   str.formula
## }












    ##   if (is it a power formula){
##       ## A power formula
##       formula <- "whatever formula defines on varnames"
##       list.formula <- rhsFormula2list(formula, usedvars)
##     } else {
##       ## Not a power formula
##       formula <- "whatever formula defines on varnames"
##       list.formula <- rhsFormula2list(formula)
##       ## Check if abbreviations are made; leads to
##       list.formula <- updateit (list.formula, usedvars)
##       formula <- list2rhsFormula(list.formula)
##     }














  
##   get.var.of.type <- function(type) {
##         varNames(data)[varTypes(data) == type]
##     }
##     used.var <- get.var.of.type(type)


##     if (!inherits(formula, "formula")) {
##         formula <- list2rhsFormula(formula)
##     }
##     list.formula <- rhsFormula2list(formula)
##     pow <- extract.power(formula)
##     if (!is.numeric(pow)) {
##         if (any(is.na(pmatch(unlist(list.formula), used.var, 
##             duplicates.ok = TRUE)))) 
##             stop("An invalid variable specification has been found\n")
##         list.formula <- lapply(list.formula, function(x) {
##             i <- pmatch(x, used.var)
##             used.var[i]
##         })
##         formula <- list2rhsFormula(list.formula)
##         str.formula <- paste(deparse(formula[[2]]), collapse = "")
##     }
##     else {
##         if (!missing(marginal)) {
##             used.var <- intersect(marginal, used.var)
##         }
##         if (pow == -1) 
##             str.formula <- paste(used.var, collapse = v.sep, 
##                 sep = "")
##         else {
##             pow <- min(c(pow, length(used.var)))
##             tmp <- selectOrder(used.var, pow)
##             str.formula <- paste(unlist(lapply(tmp, paste, collapse = v.sep)), 
##                 collapse = g.sep, sep = "")
##         }
##         formula <- formula(paste("~", str.formula, sep = ""))
##         list.formula <- rhsFormula2list(formula)
##     }
##     num.formula <- lapply(list.formula, function(l) {
##         match(l, used.var)
##     })
##     value <- list(formula = formula, str.formula = str.formula, 
##         num.formula = num.formula, list.formula = list.formula, 
##         gmData = data, varnames = used.var)
##     value



