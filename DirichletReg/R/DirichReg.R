DirichReg <- function(formula,
                      data,
                      model = c("common", "alternative"),
                      subset,
                      sub.comp,                                                 # subcompositions
                      base,
                      weights,
                      control,
                      verbosity = 0
                      ){

  this.call <- match.call()

  if(!(verbosity %in% 0:4)){
    verbosity <- 0L
    warning("invalid value for verbosity.")
  }
  storage.mode(verbosity) <- "integer"

if(verbosity > 0){
  cat("- PREPARING DATA\n")
  if(interactive()) flush.console()
}

  # checks and preliminary work
  if(missing(data)) data <- environment(formula)

  if(missing(formula)){
    stop("specification of \"formula\" is necessary.")
  } else {
    oformula <- formula
  }
  model <- match.arg(model)
  if(missing(control)){
    control <- list(sv = NULL, iterlim = 10000L, tol1 = .Machine$double.eps^(1/2), tol2 = .Machine$double.eps^(3/4))
  } else {
    if(is.null(control$sv))      control$sv       <-   NULL
    if(is.null(control$iterlim)) control$iterlim  <- 10000L
    if(is.null(control$tol1))    control$tol1     <- .Machine$double.eps^(1/2)
    if(is.null(control$tol2))    control$tol2     <- .Machine$double.eps^(3/4)
  }

#>>> get Y #####################################################################
  resp_lang <- oformula[[2L]]
  resp_char <- deparse(resp_lang)

  has_data    <- !missing(data)                                                 # "data" defined?
  Y_in_data   <- ifelse(has_data, resp_char %in% names(data), FALSE)   # response in data
  has_DR_call <- grepl("DR_data", resp_char, fixed = TRUE)             # on the fly transformation?

  if(Y_in_data){
    Y_full <- data[[resp_char]]                                                 # get response from "data"
  } else if(has_DR_call){
    Y_full <- eval(resp_lang)                                                   # get response from an on the fly transformation --- NOT recommended
    warning(paste0(strwrap("The response was transformed by DR_data() on the fly. This is not recommended, consider adapting your code.", width = getOption("width") - 9L, exdent = 9L), collapse = "\n"), call. = FALSE, immediate. = TRUE)
    oformula[[2L]] <- as.symbol("Y_full")
  } else {
    Y_full <- get(resp_char, environment(oformula))                             # get response from the parent frame
  }
  formula <- as.Formula(oformula)

  if(class(Y_full) != "DirichletRegData") stop("the response must be prepared by DR_data")

  if(has_data){
    assign(resp_char, Y_full)
  } else {
    data[[resp_char]] <- Y_full
  }



#  if(ifelse(!missing(data), deparse(oformula[[2]]) %in% names(data), FALSE)){
#    Y_full <- data[[deparse(oformula[[2]])]]                                    # GET THE FULL DRData OBJECT FOR INFORMATIONS (Y_full)
#  } else {                                                                      # handle responses in .GlobalEnv
#    Y_full <- get(deparse(oformula[[2]]), environment(oformula))
#    if(missing(data)){
#      data[[deparse(oformula[[2]])]] <- Y_full
#    } else {
#      assign(deparse(oformula[[2]]), Y_full)
#    }
#  }
#  if(class(Y_full) != "DirichletRegData") stop("the response must be prepared by DR_data")

#<<< get Y #####################################################################

  repar <- ifelse(model == "common", FALSE, TRUE)

  mf <- match.call(expand.dots = FALSE)
# if response was produced by DR_data()
  if(has_DR_call){
    mf[["formula"]][[2L]] <- as.symbol("Y_full")
  }
  mf <- mf[c(1L, match(c("formula", "data", "subset", "weights"), names(mf), 0L))]
  mf[["formula"]] <- as.Formula(mf[["formula"]])
  mf[["drop.unused.levels"]] <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf_formula <- mf
  d <- mf <- eval(mf, parent.frame())

  if("(weights)" %in% names(mf)) weights <- mf[["(weights)"]] else weights <- rep(1, nrow(mf))
  storage.mode(weights) <- "double"

  Y <- model.response(mf, "numeric")



## SUBCOMPOSITIONS
  if(missing(sub.comp)){
    sub.comp <- seq_len(ncol(Y))
  } else {
    if(length(sub.comp) == ncol(Y)) warning("no subcomposition made, because all variables were selected")
    if(length(sub.comp) == (ncol(Y) - 1)) stop("no subcomposition made, because all variables except one were selected")
    if(any((sub.comp < 1) | (sub.comp > ncol(Y)))) stop("subcompositions must contain indices of variables of the Dirichlet data object")
    y_in  <- (seq_len(ncol(Y)))[sub.comp]
    y_out <- (seq_len(ncol(Y)))[-sub.comp]
    y_in_labels <- colnames(Y)[y_in]
    y_out_labels <- paste(colnames(Y)[y_out], sep="", collapse=" + ")
    Y <- cbind(rowSums(Y[,y_out]), Y[,y_in])
    colnames(Y) <- c(y_out_labels, y_in_labels)
  }

  base <- ifelse(missing(base), attr(Y_full, "base"), base)
  if(!(base %in% seq_len(ncol(Y)))) stop("the base variable lies outside the number of variables")

  n.dim <- ncol(Y)

## SANITY CHECKS AND FORMULA EXPANSION
  if(length(formula)[1] != 1) stop("the left hand side of the model must contain one object prepared by DR_data()")

  if(!repar){   # COMMON
    if(length(formula)[2] == 1) for(i in 2:ncol(Y)) attr(formula, "rhs") <- lapply(seq_len(ncol(Y)), function(i) attr(formula, "rhs")[[1]])
    if(length(formula)[2] > ncol(Y)) stop("the right hand side must contain specifications for either one or all variables")
  } else {   # ALTERNATIVE
    if(length(formula)[2] == 1) formula <- as.Formula(formula(formula), ~ 1)
    if(length(formula)[2] > 2) stop("the right hand side can only contain one or two specifications in the alternative parametrization")
  }


  if(!repar){
    X.mats <- lapply(seq_len(ncol(Y)), function(i){ model.matrix(terms(formula, data=data, rhs=i), mf) })
    Z.mat <- NULL
    n.vars <- unlist(lapply(X.mats, ncol))
  } else {
    X.mats <- lapply(seq_len(ncol(Y)), function(i) model.matrix(terms(formula, data=data, rhs=1), mf) )
    Z.mat  <- model.matrix(terms(formula, data=data, rhs=2), mf)
    n.vars <- c(unlist(lapply(X.mats, ncol))[-1], ncol(Z.mat))
  }

################################################################################
### simplify and typecast variables ############################################
################################################################################

  # simplify dependent variable
  Y_fit <- unclass(Y)
  attributes(Y_fit) <- NULL
  dim(Y_fit) <- dim(Y)
  storage.mode(Y_fit) <- "double"

  ### typecasting
  X_fit <- lapply(X.mats, function(this_mat){
    attr(this_mat, "dimnames") <- NULL
    attr(this_mat, "assign") <- NULL
    return(this_mat)
  })
  for(i in seq_along(X_fit)) storage.mode(X_fit[[i]]) <- "double"

  if(!is.null(Z.mat)) storage.mode(Z.mat) <- "double"

  storage.mode(n.dim)  <- "integer"
  storage.mode(n.vars) <- "integer"
  storage.mode(base)   <- "integer"



if(verbosity > 0){
  cat("- COMPUTING STARTING VALUES\n")
  if(interactive()) flush.console()
}

  # compute starting values
  if(is.null(control$sv)){
    starting.vals <- get_starting_values(Y=Y_fit, X.mats=X_fit,
                       Z.mat={if(repar) as.matrix(Z.mat) else Z.mat},
                       repar=repar, base=base, weights=weights) * if(repar){ 1 } else { 1/n.dim }
  } else {
    if(length(control$sv) != n.vars) stop("wrong number of starting values supplied.")
    starting.vals <- control$sv
  }

  parametrization <- ifelse(repar, "alternative", "common")

if(verbosity > 0){
  cat("- ESTIMATING PARAMETERS\n")
  if(interactive()) flush.console()
}

  # fit and store the results
  fit.res <- DirichReg_fit(Y     = Y_fit,
                           X     = X_fit,
                           Z     = as.matrix(Z.mat),
                           sv    = starting.vals,
                           d     = n.dim,
                           k     = n.vars,
                           w     = as.vector(weights),
                           ctls  = control,
                           repar = repar,
                           base  = base,
                           vrb   = verbosity)


  varnames <- colnames(Y)

  coefs <- fit.res$estimate

  if(repar){
    names(coefs) <- unlist(as.vector(c(rep(colnames(X.mats[[1]]),n.dim-1),colnames(Z.mat))))
  } else {
    names(coefs) <- unlist(lapply(X.mats, colnames))
  }


  # FITTED VALUES
  if(repar){

    B <- matrix(0.0, nrow = n.vars[1L], ncol = n.dim)
    B[cbind(rep(seq_len(n.vars[1L]), (n.dim-1L)), rep(seq_len(n.dim)[-base], each=n.vars[1]))] <- coefs[1:((n.dim-1)*n.vars[1])]

    g <- matrix(coefs[((n.dim-1)*n.vars[1]+1):length(coefs)], ncol=1)

    XB <- exp(apply(B, 2L, function(b){ as.matrix(X.mats[[1L]]) %*% b }))
    MU <- apply(XB, 2L, function(x){ x /rowSums(XB) })

    PHI <- exp(as.matrix(Z.mat) %*% g)

    ALPHA <- apply(MU, 2L, "*", PHI)

  } else {

    B <- sapply(seq_len(n.dim), function(i){ coefs[(cumsum(c(0,n.vars))[i]+1) : cumsum(n.vars)[i]] }, simplify = FALSE)

    ALPHA <- sapply(seq_len(n.dim), function(i){ exp(as.matrix(X.mats[[i]]) %*% matrix(B[[i]], ncol=1)) })

    PHI <- rowSums(ALPHA)
    MU  <- apply(ALPHA, 2L, "/", PHI)

  }

  colnames(ALPHA) <- varnames
  colnames(MU) <- varnames

  hessian <- fit.res$hessian

  vcov <- tryCatch(solve(-fit.res$hessian),
                   error=function(x){ return(matrix(NA, nrow=nrow(hessian), ncol=ncol(hessian))) },
                   silent=TRUE)

  if(!repar){   ## COMMON
    coefnames <- apply(cbind(rep(varnames, n.vars), unlist(lapply(X.mats, colnames))), 1, paste, collapse=":")
  } else {   ## ALTERNATIVE
    coefnames <- apply(cbind(rep(c(varnames[-base], "(phi)"), n.vars), c(unlist(lapply(X.mats, colnames)[-base]), colnames(Z.mat))), 1, paste, collapse=":")
  }

  dimnames(hessian) <- list(coefnames, coefnames)
  dimnames(vcov) <- list(coefnames, coefnames)
  shortnames <- names(coefs)
  names(coefs) <- coefnames

  se <- if(!any(is.na(vcov))) sqrt(diag(vcov)) else rep(NA,length(coefs))

  res <- structure(list(
    call            = this.call,
    parametrization = parametrization,
    varnames        = varnames,
    n.vars          = n.vars,
    dims            = length(varnames),
    Y               = Y,
    X               = X.mats,
    Z               = Z.mat,
    sub.comp        = sub.comp,
    base            = base,
    weights         = weights,
    orig.resp       = Y_full,
    data            = data,
    d               = d,
    formula         = formula,
    mf_formula      = mf_formula,
    npar            = length(coefs),
    coefficients    = coefs,
    coefnames       = shortnames,
    fitted.values   = list(mu=MU,phi=PHI,alpha=ALPHA),
    logLik          = fit.res$maximum,
    vcov            = vcov,
    hessian         = hessian,
    se              = se,
    optimization    = list(convergence = fit.res$code,
                          iterations   = fit.res$iterations,
                          bfgs.it      = fit.res$bfgs.it,
                          message      = fit.res$message)
  ),
  class = "DirichletRegModel")

  # remove stuff from maxLik in the parent frame
  for(maxLik_ob in c("lastFuncGrad", "lastFuncParam")){
    if(exists(maxLik_ob, envir = parent.frame(), inherits = FALSE)) rm(list = maxLik_ob, envir = parent.frame(), inherits = FALSE)
  }
  # ... clean up the mess ...
  used_objects <- ls(all.names = TRUE)
  rm(list = c("used_objects", used_objects[used_objects != "res"]))
  # ... and make some space
  on.exit(gc(verbose = FALSE, reset = TRUE), add = TRUE)

  return(res)


}
