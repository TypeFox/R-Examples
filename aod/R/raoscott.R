raoscott <- function(formula = NULL, response = NULL, weights = NULL, group = NULL, data, pooled = FALSE, deff = NULL){
  CALL <- match.call()
### syntax control
  if(missing(formula) & missing(response) | !missing(formula) & !missing(group))
    stop("You must specify either a formula or a response with a grouping variable.")
  if(!missing(formula) & !missing(response))
    stop("You cannot specify a formula AND a response.")
  if(!missing(response) & missing(group))
    stop("You must specify a grouping variable.")

### Build the response and analysed data set

## "formula" syntax
  if(!missing(formula)){
    f <- formula
    if(length(f) != 3)
      stop("The formula ", deparse(f), " is incomplete.")
    explain <- as.character(attr(terms(f), "variables"))[-(1:2)]
    if(length(explain) > 1){
      warning("The formula contains several explanatory variables (", paste(explain, collapse = ", "), ").\n",
              "Only the first one (", explain[1], ") was used.")
      }
    f <- formula(paste(f[2], "~", explain[1]))
    mf <- model.frame(f, data)
    Y <- model.response(mf)
# response is a vector of proportions: y/n ~ group, weights = n
    if(!is.matrix(Y)){
      f <- formula(paste(deparse(f), "+", deparse(substitute(weights))))
      mf <- model.frame(formula(f), data = data)
      names(mf)[ncol(mf)] <- "(weights)"
      w <- model.weights(mf)
      Y <- w * cbind(Y, 1 - Y)
      }
# otherwise: Y is a matrix: cbind(y, n - y) ~ group
    }

## "response" syntax
  else{
    f <- paste(deparse(substitute(response)), "~", deparse(substitute(group)), "+", deparse(substitute(weights)))
    f <- formula(f)
    mf <- model.frame(formula(f), data = data)
    names(mf)[ncol(mf)] <- "(weights)"
    Y <- model.response(mf)
    w <- model.weights(mf)
# response is a matrix: cbind(y, n - y) ~ group
    if(is.matrix(Y)){
      if(!missing(weights))
        warning("The response is a matrix: weights were ignored.")
      }
# response is a vector of proportions: y/n ~ group, weights = n
    else{
      if(any(abs(Y) > 1))
        stop("Proportions cannot be greater than 1.")
      Y <- w * cbind(Y, 1 - Y)
      }
    }

## at the end of this step, response is always a matrix: cbind(y, n - y)

## check response
  n <- rowSums(Y)
  if(any(Y < 0) | any(n <= 0))
    stop("Negative counts and NULL weights are not allowed.")

## Build data set
  datan <- data.frame(n = n, y = Y[, 1], groups = factor(mf[, 2]))

### Computations

# computation of p, vratio and vbin per group
  groups <- sort(unique(datan$groups))
  lev <- levels(groups)
  N <- ntot <- ytot <- p <- vbin <- vratio <- 0
  for(i in 1:length(groups)){
    datatmp <- datan[datan$groups == lev[i], ]
    n <- datatmp$n ; y <- datatmp$y
    N[i] <- nrow(datatmp) ; ntot[i] <- sum(n) ; ytot[i] <- sum(y)
    p[i] <- ytot[i] / ntot[i] ;
    vbin[i] <- p[i] * (1 - p[i]) / ntot[i]
    vratio[i] <- N[i] * (N[i] - 1)^(-1) * ntot[i]^(-2) * sum((y - n * p[i])^2) 
    }

# computation of deff
  if(is.null(deff)){
    deff <- vratio / vbin
    if(pooled){
      df <- length(groups) - 1
      pmean <- sum(ytot) / sum(ntot)
      A <- (1 / df) * (1 / (pmean * (1 - pmean)))
      B <- (1 - ntot / sum(ntot)) * p * (1 - p) * deff
      d <- A * sum(B) ; deff <- rep(d, length(deff))
      }
    }
  else
    deff <- deff

# results
  tab <- data.frame(groups = groups, N = N, n = ntot, y = ytot, p = p, vbin = vbin, vratio = vratio, deff = deff)   
  names(tab)[1] <- as.character(f[3])
  nadj <- tab$n / tab$deff ; yadj <- tab$y / tab$deff
  padj <- sum(yadj) / sum(nadj)
  X2 <- sum((yadj - nadj * padj)^2 / (nadj * padj * (1 - padj)))

# outputs
  new(Class = "drs", CALL = CALL, tab = tab, X2 = X2)
  }

### show method for "donner" and "raoscott"
setMethod("show", signature = "drs",
  function(object){
## Build and print title
    type <- as.character(object@CALL)[1]
    ch_type <- if(type == "raoscott") "(Rao and Scott, 1993)" else "(Donner, 1989)"
    titlet <- paste("\nTest of proportion homogeneity ", ch_type, sep = "")
    undscr <- paste(rep("-", nchar(titlet)-1), collapse = "")
    cat(titlet, "\n")
    cat(undscr, "\n")
## function call
    print(object@CALL)
## computation of test details
    tab <- object@tab
    df <- nrow(tab) - 1     
    P <- pchisq(q = object@X2, df = df, lower.tail = FALSE)
    cat("N = ", sum(tab$N), " clusters, n = ", sum(tab$n), " subjects, y = ", sum(tab$y), " cases, I = ", 
        nrow(tab), " groups.\n", sep = "")
    if(type == "raoscott")
      cat("\nData and design effects:\n")
    if(type == "donner")
      cat("\nData and correction factors:\n")
    print(format(tab, digits = 4))
    if(type == "donner")
      cat("\nIntra-cluster correlation (anova estimate):", round(object@rho, digits = 4), "\n")
    cat("\nAdjusted chi-squared test:\n")
    cat("X2 = ", round(object@X2, 1), ", df = ", df, ", P(> X2) = ", round(P, digits = 4), "\n", sep = "")    
    })
