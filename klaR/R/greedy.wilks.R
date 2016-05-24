greedy.wilks <- function(X, ...) 
{
  UseMethod("greedy.wilks")
}

greedy.wilks.formula <- function(formula, data = NULL, ...) 
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- greedy.wilks(x, grouping, ...)
    res$formula <- as.formula(paste(as.character(Terms[[2]]), "~", 
                                    paste(res$results$vars, collapse = "+")))
    return(res) 
}

greedy.wilks.default <- function(X, grouping, niveau = 0.2, ...){
 
 namesset <- FALSE
 if (is.null(colnames(X))) {colnames(X) <- 1:ncol(X); namesset <- TRUE}
 
 ### Help Function   
 Lambda <- function(X, grouping){       
    # function to calculate the eigenvalues of model X with its vector of groupings grouping
  n     <- nrow(X)                      # number of observations
  n.vek <- table(grouping)              # numbers of observations in each group
  C1 <- H <- matrix(0, nrow=n, ncol=n)

  for(i in seq(along=levels(grouping))){
    einsi        <- as.numeric(grouping == levels(grouping)[i])
    C1           <- C1 + (diag(einsi) - (einsi %*% t(einsi)) / n.vek[i])
  }

  W <- t(X) %*% C1 %*% X       # "within-groups" SSP matrix
 
  H[] <- -1/n
  diag(H) <- 1 - 1/n
  TS <- t(X) %*% H %*%X        # "total" SSP matrix 
  B  <- TS-W                   # "between-groups" SSP matrix
  A  <- solve(W) %*% B
  evek <- matrix(eigen(A, symmetric = FALSE)$values, ncol = 1)  # evaluating for the eigenvalues
  evek <- Re(evek[abs(Im(evek)) < .Machine$double.eps^0.5])  # erasing complex eigenvalues (in case of existence)
  evek <- sort(evek[evek > 0], decreasing = TRUE)     # sorting the positive eigenvalues, biggest at first
  lambda <- prod(1/(1 + evek))              # calcutating Wilks' lambda
  return(lambda)                            # Output: Wilks' lambda
 }  ### END Help Function


 gpVar <- deparse(substitute(grouping))
 if(!is.null(ncol(grouping)) && ncol(grouping) > 1) 
    stop("only one grouping variable supported")
 
 X        <- as.matrix(X)               # data matrix
 grouping <- factor(grouping)           # vector of groupings
 k        <- length(levels(grouping))   # number of existing groups/populations
 n        <- nrow(X)                    # number of observations
 if(k < 2) 
    stop("at least two levels ", "required to do anything meaningful")
 if(n < 2)
    stop("n > 1 observations ", "required to do anything meaningful")
 
 # finding the first variable to accept in the model:
 p           <- ncol(X)                 # number of samples
 first.pwert <- first.fstat <- numeric(p)
 for(j in 1:p){
  e1       <- aov(X[,j] ~ grouping)
  e1sum    <- summary(e1)[[1]]
  first.pwert[j] <- e1sum[["Pr(>F)"]][1]
  first.fstat[j] <- e1sum[["F value"]][1]
 }
 min.pwert <- min(first.pwert)
 if(min.pwert < niveau) {                       # condition for stopping the forward-selection
    a               <- which.min(first.pwert)
        # pick variable which seperates the groups most, this is variable a, the one with the smallest p-value
    auswahl         <- colnames(X)[a]
    weiter          <- TRUE                     # parameter to report the necessity of another variable selection
    X.mod           <- as.matrix(X[,auswahl])   # model of the registered variables
    colnames(X.mod) <- auswahl
    X               <- X[, colnames(X) != auswahl]  # data-matrix of the selection of variables
    p               <- ncol(X)
    wilks           <- Lambda(matrix(X.mod), grouping)  # Wilks' lambda of the selected model
    Fstat <- Fstat2 <- first.fstat[a]           # approx. F-statistic of the Wilks lambda of the selected model / partial Wilks' lambda
    pwert <- pwert2 <- min.pwert                # appropriate p-value of Fstat and Fstat2
 } else {weiter <- FALSE}   
 
 # finding the next variable to accept in the model:
 while(weiter && !is.null(p) && p>0){
   auswahl     <- NULL          # vector for the names of the significant variables of the selection
   im          <- length(wilks) # dimension of selected variables in the model
   ausw.wilks  <- ausw.Fstat  <- ausw.pwert  <- ausw.Fstat2 <- ausw.pwert2 <- numeric(p)    
   ## temporary:
   # - Wilk's lambdas for the current existing X.mod + additive variables
   # - approx. F-statistics for ausw.wilks
   # - p-values for aus.Fstat
   # - approx. F-statistics for the partial Wilks' lambdas (see wpart)
   # - p-values for ausw.Fstat2
   
   for(j in 1:p){
        Mat           <- as.matrix(cbind(X.mod, X[,j]))     # the model to be tested
        colnames(Mat) <- c(colnames(X.mod), colnames(X)[j])
        e2            <- manova(Mat ~ grouping)
        e2sum         <- summary(e2, test = "Wilks")[[4]][1,]
        ausw.wilks[j] <- e2sum[2]
        ausw.Fstat[j] <- e2sum[3] 
        ausw.pwert[j] <- e2sum[6] 
        wpart         <- ausw.wilks[j]/wilks[im]    # partial Wilks' lambda: (lambda of the tested model) / (lambda of the curent exit. model X.mod)
        ausw.Fstat2[j]<- (1-wpart) / wpart * (n-k-im) / (k-1)
        ausw.pwert2[j]<- 1 - pf(ausw.Fstat2[j], k-1, n-k-im+1)
   }
   a <- which.min(ausw.wilks)[1]                # most significant variable a (with the smalles wilks' lambda)
   if(ausw.pwert2[a] < niveau){                 # condition for stopping the forward-selection
      namen2          <- colnames(X)  
      auswahl         <- namen2[a]              # see a
      namen           <- colnames(X.mod)
      X.mod           <- as.matrix(cbind(X.mod, X[,auswahl]))
      colnames(X.mod) <- c(namen, auswahl)
      X               <- as.matrix(X[, namen2 != auswahl])
      p               <- ncol (X)
      if(p == 1) colnames(X) <- namen2[namen2 != auswahl]
      if(p == 0) weiter      <- FALSE               # the case of selecting all variables: non is left to be selected 
      wilks           <- c(wilks, ausw.wilks[a])    # Wilks' lambdas of the selected models
      Fstat           <- c(Fstat, ausw.Fstat[a])    # approx. F-statistics of the Wilks' lambdas
      Fstat2          <- c(Fstat2, ausw.Fstat2[a])  # approx. F-statistics of the partial Wilks' lambdas
      pwert           <- c(pwert, ausw.pwert[a])    # appropriate p-values of Fstat
      pwert2          <- c(pwert2, ausw.pwert2[a])  # appropriate p-values of Fstat2
   } else {weiter <- FALSE}   
 }

        # Output: The names of the variables in the order of its selection, the appropriate Wilks' lambdas,
        #         the approximated F-statistics and its p-values, the F-statistics of the partial Wilks's lambdas 
        #         and its p-values.
 if(!exists("X.mod"))
    stop("unable to perform required calculations, perhaps not enough observations?")
 vars <- colnames(X.mod)
 if (namesset) vars<-as.numeric(vars)
 resDat <- data.frame(vars = vars, Wilks.lambda = wilks, F.statistics.overall = Fstat, 
                      p.value.overall = pwert, F.statistics.diff = Fstat2, p.value.diff = pwert2)
 resForm <- as.formula(paste(gpVar, "~", paste(vars, collapse = "+")))
 res <- list(results = resDat, formula = resForm)
 class(res) <- "greedy.wilks"
 return(res) 
}


print.greedy.wilks <- function(x, ...)
{
  cat("Formula containing included variables:", "\n\n")
  print(x$formula, ...)
  cat("\n\nValues calculated in each step of the selection procedure:", "\n\n")
  print(x$results, ...)
  invisible(x)
}
