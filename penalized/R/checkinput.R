########################################################
# A function that does the common input checking of
# the functions penalized, cvl, profL1, profL2, optL1 and optL2.
########################################################
.checkinput <- function(call, env) {
                      
  # Functions to extract the original input variables
  call <- as.list(call)
  input <- function(str) eval(call[[str]], env)
  input.data <- function(str) eval(call[[str]], data, env)
  missing <- function(str) !str %in% names(call)

  # determine the response
  if (!missing("data")) {
     data <- input("data")
       } else {data <- NULL}
  
  response <- input.data("response")
  if (is(response, "formula")) {
    formula.response <- response
    if (missing("data"))
      response <- eval(attr(terms(response), "variables"), environment(response))[[attr(terms(response), "response")]]
    else
      response <- eval(attr(terms(response, data=data), "variables"), data, environment(response))[[attr(terms(response, data=data), "response")]]
  } else {
    formula.response <- NULL
  }
  
  ##determine whether to perform fused lasso or not
  fusedl <- if (missing("fusedl")){fusedl <- FALSE} else {fusedl <- input("fusedl")}

  # determine the model if missing
  if (missing("model")) {
    if (is(response, "Surv")) model <- "cox"
    else if (is.logical(response) || all(response %in% 0:1) || is.factor(response)) model <- "logistic"
    else if (is.numeric(response)) model <- "linear"
    else stop("Model could not be determined from the input. Please specify the model.")
  } else {
    model <- match.arg(input("model"), c("cox", "logistic", "linear", "poisson"))
  }
  
  # coerce factor response to 0,1
  if (is.factor(response)) {
    if (length(levels(response)) != 2)
      stop("response has ", length(levels(response)), " levels. Only two-level factor response supported.")
    response <- response != levels(response)[1]
  }
    
  # check positivity of survival response
  if (is(response, "Surv") && any(as.vector(response) < 0))
    stop("negative survival times")
  # avoid problems with Design package
  if (is(response, "Surv") && "labelled" %in% class(response))
    class(response) <- "Surv"

  # determine penalized and unpenalized
  if (!missing("penalized")) {
     penalized <- input("penalized")
       } else {
    if (!is.null(formula.response)) {
      penalized <- formula.response
      formula.response <- NULL
    } else
      stop("argument \"penalized\" is missing, with no default")
  }
  if (!missing("unpenalized")) {
    unpenalized <- input.data("unpenalized")
    
  } else {
    if (!is.null(formula.response)) {
      unpenalized <- formula.response
      formula.response <- NULL
    } else {
      unpenalized <- response ~ 1
    }
  }

  # Has the response formula been used?
  if (!is.null(formula.response))
    warning("right hand side of response formula ignored")
  formula.input <- list()     # stores input formula

  # coerce unpenalized into a matrix and find the offset term
  offset <- NULL
  strata <- NULL
  if (is.data.frame(unpenalized) || is.vector(unpenalized)) {
    if (all(sapply(unpenalized, is.numeric))) {
      unpenalized <- as.matrix(unpenalized)
    } else {
      stop("argument \"unpenalized\" could not be coerced into a matrix")
    }
  }
  if (is(unpenalized, "formula")) {
    formula.input$unpenalized <- unpenalized
    if (missing("data")) {
      tup <- terms(unpenalized, specials = "strata") 
      # prevent problems for input ~1 or ~0:
      if ((attr(tup, "response") == 0) && (length(attr(tup, "term.labels")) == 0)) {
        if (attr(tup, "intercept") == 1)
          unpenalized <- terms(response ~ 1)
        else
          unpenalized <- terms(response ~ 0)
      } else {
        offset <- model.offset(model.frame(unpenalized))
        unpenalized <- tup
      }
    } else {
      offset <- model.offset(model.frame(unpenalized, data=data))
      unpenalized <- terms(unpenalized, specials="strata", data=data)
    }
    # Cox: suppress intercept if necessary and extract strata
    if (model == "cox") {
      if (length(attr(unpenalized, "specials")$strata) > 0) {
        strata <- untangle.specials(unpenalized, "strata", 1)
        strata.nrs <- strata$terms                            # indices of the strata variables in the terms object
        strata.nrs2 <- attr(unpenalized, "specials")$strata   # indices of the strata variables in attr(unpenalized, "variables")
        if (missing("data"))
          strata <- strata(eval(attr(unpenalized, "variables"), environment(unpenalized))[strata.nrs2], shortlabel=TRUE)
        else
          strata <- strata(eval(attr(unpenalized, "variables"), data, environment(unpenalized))[strata.nrs2], shortlabel=TRUE)
        unpenalized <- unpenalized[-strata.nrs]
      } 
      attr(unpenalized, "intercept") <- 1
      # prevent problems in case of only unpenalized = ~strata() only
      if ((attr(unpenalized, "response") == 0) && (length(attr(unpenalized, "term.labels")) == 0)) 
        unpenalized <- terms(response ~ 1)
    }
    unpenalized <- model.matrix(unpenalized, data)
    if (model == "cox") 
      unpenalized <- unpenalized[,-1, drop=FALSE]
  }

  # coerce penalized into a matrix
  if (is.data.frame(penalized) || is.vector(penalized))
    if (all(sapply(penalized, is.numeric))) {
      penalized <- as.matrix(penalized)
    } else {
      stop("argument \"penalized\" could not be coerced into a matrix")
    }
  if (is(penalized, "formula")) {
    formula.input$penalized <- penalized
    oldcontrasts <- unlist(options("contrasts"))
    options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    if (missing("data"))
      penalized <- terms(penalized)
    else
      penalized <- terms(penalized, data=data)
    # prevent problems for input ~1 or ~0:
    if (length(attr(penalized, "term.labels")) == 0) 
      penalized <- terms(response ~ 1)
    attr(penalized, "intercept") <- 1
    penalized <- model.matrix(penalized, data)
    penalized <- penalized[,-1,drop=FALSE]
    options(contrasts = oldcontrasts)
  }

  # check missing values in penalized
  if (any(is.na(penalized)))
    stop("missing values in \"penalized\" argument")

  # check dimensions of response, penalized and unpenalized
  if (model == "cox") {
    if (attr(response, "type") == "right")
      n <- length(response)/2 
    else if (attr(response, "type") == "counting")
      n <- length(response)/3
  } else {
    n <- length(response)
  }
  if (nrow(penalized) != n) {
    stop("the length of \"response\" (",n, ") does not match the row count of \"penalized\" (", nrow(penalized), ")")
  }
  if (nrow(unpenalized) != n) {
    stop("the length of \"response\" (",n, ") does not match the row count of \"unpenalized\" (", nrow(unpenalized), ")")
  }

  # expand positive if necessary
  if (!missing("positive"))
    if (length(input("positive")) == 1)
      if (input("positive"))
        positive <- c(logical(ncol(unpenalized)), !logical(ncol(penalized)))
      else
        positive <- c(logical(ncol(unpenalized)), logical(ncol(penalized)))
    else if (length(input("positive")) == ncol(penalized))
      positive <- c(logical(ncol(unpenalized)), input("positive"))
    else
      stop("length of \"positive\" does not match column count of \"penalized\"")
  else
    positive <- logical(ncol(unpenalized) + ncol(penalized))

  # get the value of startbeta
  if (missing("startbeta"))
    startbeta <- numeric(ncol(penalized))
  else {
    startbeta <- input("startbeta")
    if (length(startbeta) != ncol(penalized))
      stop("The length of \"startbeta\" (", length(startbeta), ") does not match the column count of \"penalized\" (", ncol(penalized), ")")
  }
  if (is.null(names(startbeta)))
    names(startbeta) <- colnames(penalized)

  # get the value of startgamma
  if (ncol(unpenalized) > 0) {
    if (is.null(offset)) {
      nullgamma <- switch(model,
        cox = coefficients(coxph(response ~ unpenalized)),
        logistic = coefficients(glm(response ~ 0 + unpenalized, family = binomial)),
        linear = coefficients(lm(response ~ 0 + unpenalized)),
        poisson = coefficients(glm(response ~ 0 + unpenalized, family = poisson))
      )
    } else {
      nullgamma <- switch(model,
        cox = coefficients(coxph(response ~ offset(offset) + unpenalized)),
        logistic = coefficients(glm(response ~ 0 + offset(offset) + unpenalized, family = binomial)),
        linear = coefficients(lm(response ~ 0 + offset(offset) + unpenalized)),
        poisson = coefficients(glm(response ~ 0 + offset(offset) + unpenalized, family = poisson))
      )
    }
    names(nullgamma) <- colnames(unpenalized)
  } else nullgamma <- numeric(0)
  if (missing("startgamma")) {
    startgamma <- nullgamma
  } else {
    startgamma <- input("startgamma")
    if (length(startgamma) != ncol(unpenalized))
      stop("The length of \"startgamma\" (", length(startgamma), ") does not match the column count of \"unpenalized\" (", ncol(unpenalized), ")")
  }
  if (is.null(names(startgamma)))
    names(startgamma) <- colnames(unpenalized)


  # orthogonalize penalized with respect to unpenalized
  if (ncol(unpenalized) > 0) {
    orthogonalizer <- solve(crossprod(unpenalized), crossprod(unpenalized, penalized))
    penalized <- penalized - unpenalized %*% orthogonalizer
  } else {
    orthogonalizer <- matrix(,0,ncol(penalized))
  }
  
  ##getchr
  

    if(is.logical(fusedl)){
      if(fusedl){
     if(ncol(unpenalized)>0){
             chr = c(rep(0,ncol(unpenalized)),rep(1,ncol(penalized)))
             names(chr) = c(colnames(unpenalized),colnames(penalized))
   }else{chr = rep(1,ncol(penalized))
        names(chr) = colnames(penalized)}
    } else {chr = rep(0,ncol(penalized))}
    } else if(!is.logical(fusedl)) {
               if(ncol(unpenalized)>0){
               chr = as.numeric(fusedl)
               chr = c(rep(0,ncol(unpenalized)),chr)
             names(chr) = c(colnames(unpenalized),colnames(penalized))
   }else {names(chr) = colnames(penalized)}
    }

  # Join penalized and unpenalized together
  X <- cbind(unpenalized, penalized)
  beta <- c(startgamma, startbeta)

  # stabilize/standardize
  standardize <- if (missing("standardize")){ FALSE } else {input("standardize")}
  fusedl <- if (is.logical(fusedl)){fusedl <- fusedl} else {fusedl <- TRUE}
  

  
    vars <- apply(X,2,var) * (n-1)/n
    vars[vars == 0] <- 1
    sds <- sqrt(vars)

 if(is.logical(fusedl) && !fusedl){
    beta[beta != 0] <- beta[beta != 0] * sds[beta != 0]
    nullgamma <- nullgamma * sds[1:length(nullgamma)]
    X <- X / matrix(sds, nrow(X), ncol(X), byrow=T)

  # find baselambda1 and baselambda2
  # This lambda1 and lambda2 for unit input lambda1=1 and lambda2=1
  if (standardize) {
    baselambda1 <- c(numeric(ncol(unpenalized)), rep(1, ncol(penalized)))
    baselambda2 <- c(numeric(ncol(unpenalized)), rep(1, ncol(penalized)))
  } else {
    sel <- ncol(unpenalized) + 1:ncol(penalized)
    baselambda1 <- c(numeric(ncol(unpenalized)), 1/sds[sel])
    baselambda2 <- c(numeric(ncol(unpenalized)), 1/vars[sel])
  }
}  else if(is.logical(fusedl) && fusedl){baselambda1 <- c(numeric(ncol(unpenalized)), rep(1, ncol(penalized)))
         baselambda2 <- c(numeric(ncol(unpenalized)), rep(1, ncol(penalized)))
         }
              
  return(list(
    fusedl = fusedl,
    chr = chr, 
    response = response,
    X = X, 
    beta = beta, 
    weights = sds, 
    baselambda1 = baselambda1,
    baselambda2 = baselambda2,
    positive = positive,
    orthogonalizer = orthogonalizer, 
    model = model, 
    nullgamma = nullgamma,
    offset = offset,
    strata = strata,
    formula = formula.input
  ))
}

# Switch functions to choose the appropriate function for a specific model
.modelswitch <- function(model, response, offset, strata) {
  switch(model,
    cox = .coxfit(response, offset, strata),
    logistic = .logitfit(response, offset),
    linear = .lmfit(response, offset),
    poisson = .poissonfit(response, offset)
  )
}

.predictswitch <- function(model, predictions, groups) {
  switch(model, 
    cox = .coxmerge(predictions, groups),
    logistic = .logitmerge(predictions, groups),
    linear = .lmmerge(predictions, groups),
    poisson = .poissonmerge(predictions, groups)
  )
}