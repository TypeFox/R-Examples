################################################################################
################################################################################
################################## multiPIM.R ##################################
################################################################################
################################################################################

## Author: Stephan Ritter

## File containg definition of main function for multiPIM package.
## In addition, 4 character vectors are made available to the user.

## The character vectors are:

## all.bin.cands
## all.cont.cands
## default.bin.cands
## default.cont.cands

## The function is called multiPIM

################################################################################
############################## define candidates ###############################
################################################################################

## these objects make the lists of valid candidates available to the user

## IF YOU CHANGE THESE, remember to change the versions inside the multiPIM
## function as well. Also, there should not be any overlap between elements of
## all.bin.cands and all.cont.cands

all.bin.cands <- c("polyclass", "penalized.bin", "main.terms.logistic",
                   "rpart.bin")
all.cont.cands <- c("polymars", "lars", "main.terms.linear", "penalized.cont",
                    "rpart.cont")
default.bin.cands <- all.bin.cands[1:3]
default.cont.cands <- all.cont.cands[1:3]


################################################################################
################################# main function ################################
################################################################################

multiPIM <- function(Y, A, W = NULL, ## data frames
                     estimator = c("TMLE", "DR-IPCW", "IPCW", "G-COMP"),
                     g.method = "main.terms.logistic", g.sl.cands = NULL,
                     g.num.folds = NULL, g.num.splits = NULL,
                     Q.method = "sl", Q.sl.cands = "default", ## ignored if
                                                       ## estimator is "IPCW"
                     Q.num.folds = 5, Q.num.splits = 1, ## ignored if
                                                       ## estimator is "IPCW"
                     Q.type = NULL, ## override the default mechanism for
                                    ## choosing between logistic- vs.
                                    ## linear-type regression
                                    ## This may not be NULL if check.input is
                                    ## FALSE and estimator is "DR-IPCW". Will be
                                    ## ignored if estimator is "IPCW"
                     adjust.for.other.As = TRUE, # ignored if ncol(A) == 1
                     truncate = 0.05, ## can also set this to FALSE
                     return.final.models = TRUE,
                     na.action, ## not yet implemented
                     check.input = TRUE, ## may set to false when when running
                                         ## the bootstrap
                     verbose = FALSE,
                     extra.cands = NULL,
                     standardize = TRUE, ## applies to lars and penalized cands
                     ...) {

  ## Save start time

  start.time <- Sys.time()
  
################################ check input ###################################

  ## check check.input

  if( !(identical(check.input, TRUE) || identical(check.input, FALSE)) )
    stop("the check.input argument may only be either TRUE or FALSE")

  ## check Q.type

  if(identical(Q.type, "continuous.outcome")
     || identical(Q.type, "binary.outcome")) {
    Q.type.override <- TRUE
  } else Q.type.override <- FALSE

  ## check estimator

  estimator = match.arg(estimator)

  if(estimator == "IPCW") {
    do.Q <- FALSE
  } else do.Q <- TRUE

  if(estimator == "G-COMP") {
    do.g <- FALSE
  } else do.g <- TRUE
     
  ## names of the built-in candidates

  all.bin.cands <- c("polyclass", "penalized.bin", "main.terms.logistic",
                     "rpart.bin")
  all.cont.cands <- c("polymars", "lars", "main.terms.linear", "penalized.cont",
                      "rpart.cont") 
  default.bin.cands <- all.bin.cands[1:3]
  default.cont.cands <- all.cont.cands[1:3]

  if(!check.input) {

    ## make sure either Q.type was specified or estimator is "IPCW"

    if(!Q.type.override && do.Q)
      stop('If check.input is FALSE and estimator is not "IPCW"\n',
           'then Q.type must be specified as either of\n',
           '"continuous.outcome" or "binary.outcome"')

  } else { ## check the rest of the input

    if(!(identical(verbose, TRUE) || identical(verbose, FALSE)))
      stop("the argument verbose must be given as either TRUE or FALSE")

    if(verbose) cat("checking input\n")

    ## find out whether Y or A is missing

    if(missing(Y)) stop("must supply a data frame Y containing outcomes")
    if(missing(A)) stop("must supply a data frame A containing exposures ",
                        "or treatments")

    ## check Y ##

    ## is it a data frame?

    if(!identical(class(Y), "data.frame")) stop("Y must be a data frame")

    ## check dimensions of Y

    if(ncol(Y) < 1) stop("Y must have at least one column")
    if(nrow(Y) < 2) stop("Y must have at least two rows")

    ## test for missing values

    if(any(is.na(Y))) stop("Y has missing values. This is not yet supported.")

    ## check class/mode of the columns make sure mode is numeric and class
    ## isn't factor

    for(j in 1:ncol(Y)) {

      if(identical(class(Y[[j]]), "factor"))
        stop('Y may not contain columns of class "factor"')
      if(!identical(mode(Y[[j]]), "numeric"))
        stop('all columns of Y must have mode "numeric"')

    }

    ## check each column has a name

    if(!identical(length(names(Y)), ncol(Y)) || any(is.na(names(Y))))
      stop("each column of the data frame Y must have a valid name")

    ## check extra.cands

    if(!is.null(extra.cands)) {

      if(!is.list(extra.cands) || is.null(names(extra.cands)) ||
         any(sapply(extra.cands, class) != "function"))
        stop("if extra.cands is non-null, it must be a named list of functions")

      all.cands <- c(names(extra.cands), all.bin.cands, all.cont.cands)

      if(!identical(all.cands, make.names(all.cands, unique = TRUE)))
        stop("the names of extra.cands must be syntactically valid\n",
             "(e.g. they must start with letters), they must be unique,\n",
             "(i.e. no repetition), and they must be different from\n",
             "the elements of all.bin.cands and all.cont.cands")
    }

    ## go back and check on Q.type (but only if do.Q is TRUE)
    if(do.Q) {

      if(!Q.type.override) {

        if(!is.null(Q.type))
          stop('if estimator is not "IPCW", then Q.type must be:\n',
               'either NULL or "continuous.outcome" or "binary.outcome"')

        ## Find out if the outcomes are binary (i.e. just have 0 or 1) or
        ## continuous.
        ## If some are binary and others are not throw a warning, but continue
        ## as if all were continuous. Set Q.type accordingly.

        col.is.binary <- rep(FALSE, ncol(Y))

        for(j in 1:ncol(Y))
          if( all( Y[[j]]==0 | Y[[j]]==1 ) )
            col.is.binary[j] <- TRUE

        if(all(col.is.binary)) {
          Q.type <- "binary.outcome"
        } else {
          Q.type <- "continuous.outcome"
          if(any(col.is.binary))
            warning("Some columns of Y are binary (0/1) and others are not.\n",
                    "You may want to run multiPIM twice,\ni.e. separately ",
                    "on the binary and on the continuous outcomes.")
        }

        rm(col.is.binary)

      } ## end if not Q.type.override


      ## check Q.method ##

      ## first check that length is 1

      if(length(Q.method) != 1)
        stop("Q.method must be a length 1 character vector giving\n",
             "a regression method to use for modeling Q")

      if(Q.method != "sl") {

        if(identical(Q.type, "binary.outcome")) {

          if(!(Q.method %in% c(all.bin.cands, names(extra.cands))))
            stop("The type of regression for Q has been selected as\n",
                 "binary outcome. Thus Q.method must be a length one\n",
                 'character vector whose value is either "sl", or an element\n',
                 "of the vector all.bin.cands, or one of the names of the\n",
                 "extra.cands list of functions.")

        } else if(identical(Q.type, "continuous.outcome")) {

          if(!(Q.method %in% c(all.cont.cands, names(extra.cands))))
            stop("The type of regression for Q has been selected as\n",
                 "continuous outcome. Thus Q.method must be a length one\n",
                 'character vector whose value is either "sl", or an element\n',
                 "of the vector all.cont.cands, or one of the names of the\n",
                 "extra.cands list of functions.")

        } else stop("Internal Error #1: Q.type has some weird value")

      } else { ## else Q.method should be "sl"
               ## so need to check Q.sl.cands, Q.num.folds and Q.num.splits

        if(length(Q.sl.cands) == 1) {

          if(!(identical(Q.sl.cands, "all")
               || identical(Q.sl.cands, "default")))
            stop('if Q.sl.cands has length 1, its value must be\n',
                 'either "all" or "default".\n',
                 'Otherwise, at least 2 candidates must be specified.\n',
                 'Type ?multiPIM for more details.')

        } else if(length(Q.sl.cands) > 1) {

          if(!identical(unique(Q.sl.cands), Q.sl.cands))
            stop("there must be no repetition of candidates ",
                 "in the Q.sl.cands argument")

          if(identical(Q.type, "binary.outcome")) {

            if(!all(Q.sl.cands %in% c(all.bin.cands, names(extra.cands))))
              stop('the type of regression for Q has been selected as\n',
                   'binary outcome type (either because this was specified\n',
                   'by the Q.type argument, or because of the contents of',
                   'Y).\n',
                   'Thus, the elements of Q.sl.cands must match with those of',
                   '\nthe vector all.bin.cands, or with the names of the\n',
                   'extra.cands list of functions, if it is non-null.')

          } else if(identical(Q.type, "continuous.outcome")) {

            if(!all(Q.sl.cands %in% c(all.cont.cands, names(extra.cands))))
              stop('the type of regression for Q has been selected as\n',
                   'continuous outcome type (either because this was\n',
                   'specified by the Q.type argument, or because of the\n',
                   'contents of Y). Thus, the elements of Q.sl.cands\n',
                   'must match with those of the vector all.cont.cands,\n',
                   'or with the names of the extra.cands list of functions,\n',
                   'if it is non-null.')

          } else stop("Internal Error #2: Q.type has some weird value")

        } else stop('if Q.method is "sl" then the candidates to use\n',
                    'must be specified in the Q.sl.cands argument,\n',
                    'which should be a character vector specifying\n',
                    'at least two candidates. The length one character',
                    ' vectors\n "all" and "default" are also allowed.')

        ## check the Q.num.folds and Q.num.splits

        if(!identical(mode(Q.num.folds), "numeric")
           || length(Q.num.folds) != 1 || is.na(Q.num.folds)
           || Q.num.folds < 2 || Q.num.folds > nrow(Y)
           || (Q.num.folds %% 1) != 0 )
          stop('if Q.method is "sl", then Q.num.folds\n',
               "must be a length 1 numeric vector with an integer value\n",
               "greater than or equal to 2 and less than or equal to nrow(Y)")

        if(!identical(mode(Q.num.splits), "numeric")
           || length(Q.num.splits) != 1 || is.na(Q.num.splits)
           || Q.num.splits < 1 || (Q.num.splits %% 1) != 0 )
          stop("Q.num.splits must be a single integer value that is >= 1")

      } ## end else for which Q.method should be "sl"
      
    } ## end if do.Q

    ## now check A ##

    ## is it a data frame?

    if(!identical(class(A), "data.frame")) stop("A must be a data frame")

    ## check dimensions of A

    if(ncol(A) < 1) stop("A must have at least one column")

    if(nrow(Y) != nrow(A)) stop("Y and A must have the same number of rows")

    ## test A for missing values

    if(any(is.na(A))) stop("A has missing values. This is not yet supported.")

    ## check class/mode of the columns make sure mode is numeric and class
    ## isn't factor

    for(j in 1:ncol(A)) {

      if(identical(class(A[[j]]), "factor"))
        stop('A may not contain columns of class "factor"')
      if(!identical(mode(A[[j]]), "numeric"))
        stop('all columns of A must have mode "numeric"')
    }

    ## check to make sure that all values in A are either 0 or 1

    if( !all( A==0 | A==1 ) )
      stop("All values in A must be 0 or 1")

    ## check whether each column has a name

    if(!identical(length(names(A)), ncol(A)) || any(is.na(names(A))))
      stop("each column of the data frame A must have a non-NA name")

    ## check the names of Y and A (test with make.names)

    if(!identical(c(names(Y), names(A)),
                  make.names(c(names(Y), names(A)), unique = TRUE)))
      stop("The names of Y and A must be valid variable names.\n",
           "They must all start with letters and not numbers,\n",
           "and they must all be unique.")

    ## check adjust.for.other.As if ncol(A) > 1

    if(ncol(A) > 1) {

      if(!(identical(adjust.for.other.As, TRUE) ||
           identical(adjust.for.other.As, FALSE)))
        stop('if ncol(A) is greater than 1, then adjust.for.other.As\nmust be ',
             'a single logical value, either TRUE or FALSE')
    }

    ## check stuff for g modeling if do.g

    if(do.g) {

      ## check g.method

      if(length(g.method) != 1 ||
         !(g.method %in% c(all.bin.cands, "sl", names(extra.cands))))
        stop('If estimator is not "G-comp", then g.method must be \n',
             'a length 1 character vector which is either\n',
             'the name of a binary outcome candidate,\n',
             'or "sl" for super learning')
      
      ## check g.sl.cands, g.num.folds and g.num.splits if g.method is "sl"

      if(g.method == "sl") {

        if(!all(g.sl.cands %in% c(all.bin.cands, names(extra.cands))))
          stop('if g.method is "sl" then g.sl.cands must be specified.\n',
               'See help file for possible values')

        if(!identical(unique(g.sl.cands), g.sl.cands) ||
           length(g.sl.cands) < 2)
          stop('if g.method is "sl", then g.sl.cands must contain at least\n',
               "2 different valid candidate names,\nand no repetitions are ",
               "allowed")

        if( !identical(mode(g.num.folds), "numeric") || length(g.num.folds) != 1
           || is.na(g.num.folds) || g.num.folds < 2 || g.num.folds > nrow(Y)
           || (g.num.folds %% 1) != 0 )
          stop('if g.method is "sl", then g.num.folds must be\n',
               "an integer value greater than or equal to 2\n",
               "and less than or equal to nrow(Y)")

        if(!identical(mode(g.num.splits), "numeric")
           || length(g.num.splits) != 1
           || is.na(g.num.splits) || g.num.splits < 1
           || (g.num.splits %% 1) != 0 )
          stop("g.num.splits must be a single integer value that is >= 1")
      }
    } ## end if do.g
    
    ## check W ##

    if(!is.null(W)) {

      ## if we are in this if block, that means a non-null W was supplied

      ## make sure W is a data frame

      if(!identical(class(W), "data.frame"))
        stop("if a W is supplied, it must be a data frame")

      ## check dimensions

      if(ncol(W) < 1)
        stop("if a W is supplied, it must have at least one column")

      if(!identical(nrow(W), nrow(Y)))
        stop("if a W is supplied, it must have the same number of rows\n",
             "as Y and A")

      ## test W for missing values

      if(any(is.na(W))) stop("W has missing values. This is not yet supported.")

      ## check class/mode of the columns of W -- may be numeric or factor

      for(j in 1:ncol(W)) {

        if(!identical(mode(W[[j]]), "numeric"))
          stop('The columns of W must have mode "numeric"\n',
               '(class may be "numeric", "integer" or "factor")')

        if(!(identical(class(W[[j]]), "numeric")
             || identical(class(W[[j]]), "integer")
             || ("factor" %in% class(W[[j]])))) ## allow ordered factors
          stop('Columns of W must have class "numeric", "integer" or "factor"')

      }

      ## check names of W using make.names and also compare to A and Y
      ## throw error in either case

      if(!identical(names(W), make.names(names(W), unique = TRUE)))
        stop("the names of W must be valid variable names:\n",
             "they must start with letters, not numbers, ",
             "and they must be unique")

      if(any( names(W) %in% c(names(Y), names(A)) ))
        stop("W has one or more names which are the same as those of A ",
             "and/or Y")

    } else { ## else W is NULL

      if(ncol(A) == 1 || !adjust.for.other.As)
        stop("if A has only one column, or if adjust.for.other.As is FALSE,\n",
             "then you must supply a data frame W with at least one column.")
    }

    ## check truncate (but ignore if !do.g)

    if(do.g) {
      if(!identical(truncate, FALSE)) {

        if(length(truncate) != 1 || mode(truncate) != "numeric" || truncate <= 0
           || truncate > 0.5)
          stop("the truncate argument must be either FALSE or\n",
               "a number greater than 0 and less than 0.5")
      }
    }

    ## check return.final.models

    if(!(identical(return.final.models, TRUE)
         || identical(return.final.models, FALSE)))
      stop("return.final.models must be either TRUE or FALSE")

    ## check standardize
    
    if(!(identical(standardize, TRUE)
         || identical(standardize, FALSE)))
      stop("the standardize argument must be either TRUE or FALSE")
        
  } ## end else check rest of input

################################################################################
############################# end input checking ###############################
################################################################################



################################################################################
################################# functions ####################################
################################################################################


############## split function which won't warn when nobs %% v != 0 #############

  randomized.split <- function(nobs, v) {

    return(split(sample(nobs), rep(1:v, length = nobs)))

    ## got the idea for using rep with a length argument from susan
    ## think she got it from eric

  }

########## function for evaluating candidates in binary outcome case ###########

  get.log.likelihood <- function(predicted.probs, response) {

    ## assumes "response" is a 0/1 variable

    return(sum(log(ifelse(response == 1, predicted.probs,
                          1 - predicted.probs))))

  }

######## function for evaluating candidates in continuous outcome case #########

  get.MSE <- function(estimates, actual.outcome) {

    return( mean( (estimates - actual.outcome) ^2 ) )

  }

################################################################################
############################# start actual stuff ###############################
################################################################################

## convert Q.sl.cands == "all" or "default" to actual candidates if applicable

  if(do.Q && Q.method == "sl") {

    if(identical(Q.sl.cands, "all")) {

      if(identical(Q.type, "binary.outcome")) {
        Q.sl.cands <- all.bin.cands
      } else Q.sl.cands <- all.cont.cands

    } else if(identical(Q.sl.cands, "default")) {

      if(identical(Q.type, "binary.outcome")) {
        Q.sl.cands <- default.bin.cands
      } else Q.sl.cands <- default.cont.cands
    }
  }

####################### function to expand factors #############################

  ## input is a data frame, output will be a matrix

  expand.factors <- function(X) {

    nobs <- nrow(X)
    
    ## make a dummy y

    y <- numeric(nobs)

    data1 <- cbind(y = y, X)

    f1 <- formula(data1)

    ## factors will be expanded according to getOption("contrasts")
    
    mm1 <- model.matrix(f1, data = data1)[, -1, drop = FALSE]
    
 } 
############### Generate list of functions for the candidates ##################

  if(!is.null(extra.cands)) {
    candidate.functions <- extra.cands
  } else candidate.functions <- vector("list", length = 0)

  rm(extra.cands)

  ## Add binary outcome candidates to candidate function list

  if(g.method == "polyclass"
     || (g.method == "sl" && "polyclass" %in% g.sl.cands)
     || (do.Q && identical(Q.type, "binary.outcome")
         && (Q.method == "polyclass"
             || (Q.method == "sl" && "polyclass" %in% Q.sl.cands)))) {

    candidate.functions$polyclass <- function(X, Y, newX, force, standardize) {

      result <- vector("list", length = 4)
      names(result) <- c("preds", "fitted", "main.model", "force.model")
      class(result) <- "polyclass.result"

      expanded.X <- expand.factors(X)
      expanded.newX <- expand.factors(newX)

      result$main.model <- polyclass(as.vector(Y), expanded.X)

      if(missing(force) || force %in% result$main.model$fcts[, 1]) {
        result$fitted <- ppolyclass(cov = expanded.X,
                                   fit = result$main.model)[,2]
        result$preds <- ppolyclass(cov = expanded.newX,
                                   fit = result$main.model)[,2]
        return(result)
      }
      
      if(nrow(result$main.model$fcts) == 1
         && all(is.na(result$main.model$fcts[1, 1:4]))) {

        glm.data <- cbind(Y, X[, force, drop = FALSE])
        glm.formula <- paste(names(Y), "~", names(X)[force], sep = "")
        result$force.model <- glm(glm.formula, data = glm.data,
                                  family = binomial, model = FALSE,
                                  y = FALSE)
        result$fitted <- predict(result$force.model, type = "response")
        result$preds <- predict(result$force.model, type = "response",
                                newdata = newX[, force, drop = FALSE])
        return(result)
      }

      ## Truncate so that you don't get Inf's that crash glm

      pred.probs <- ppolyclass(cov = expanded.X, fit = result$main.model)[,2]
      pred.probs[pred.probs < 0.001] <- 0.001
      pred.probs[pred.probs > 0.999] <- 0.999

      glm.data <- cbind(Y, X[, force, drop = FALSE],
                        preds = qlogis(pred.probs))

      glm.formula <- paste(names(Y), "~", names(X)[force], "*preds", sep = "")
      result$force.model <- glm(glm.formula, data = glm.data, family = binomial,
                                model = FALSE, y = FALSE)

      result$fitted <- predict(result$force.model, type = "response")

      ## Don't truncate the predictions on newX
      ## (predict.glm can handle Inf's without crashing)

      result$preds <- predict(result$force.model, type = "response",
                              newdata = cbind(newX[, force, drop = FALSE],
                                preds = qlogis(ppolyclass(cov = expanded.newX,
                                               fit = result$main.model)[,2])))
      return(result)
    }
  }

  if(g.method == "penalized.bin"
     || (g.method == "sl" && "penalized.bin" %in% g.sl.cands)
     || (do.Q && identical(Q.type, "binary.outcome")
         && (Q.method == "penalized.bin"
             || (Q.method == "sl" && "penalized.bin" %in% Q.sl.cands)))) {

    candidate.functions$penalized.bin <- function(X, Y, newX, force,
                                                  standardize) {

      result <- vector("list", length = 4)
      names(result) <- c("preds", "fitted", "main.model", "force.model")
      class(result) <- "penalized.bin.result"

      expanded.X <- data.frame(expand.factors(X))
      expanded.newX <- data.frame(expand.factors(newX))

      if(missing(force)) {

        pen.formula <- as.formula(paste("~", paste(names(expanded.X),
                                                   collapse = "+"),
                                        sep = ""))
      
        lambdas.and.cvls <- profL1(Y[, 1], pen.formula, data = expanded.X,
                                   lambda2 = 0, model = "logistic", fold = 5,
                                   standardize = standardize, steps = 10,
                                   trace = FALSE)[c("lambda", "cvl")]
        lambda1.to.use <-
          lambdas.and.cvls$lambda[which.max(lambdas.and.cvls$cvl)]

        result$main.model <- penalized(Y[, 1], pen.formula,
                                       data = expanded.X, trace = FALSE,
                                       lambda1 = lambda1.to.use,
                                        standardize = standardize,
                                       model = "logistic")
        result$preds <- predict(result$main.model, penalized = pen.formula,
                                data = expanded.newX)
        result$fitted <- predict(result$main.model, penalized = pen.formula,
                                 data = expanded.X)
        
        
        return(result)
      }

      pen.formula <- as.formula(paste("~", paste(names(expanded.X)[-force],
                                                 collapse = "+"),
                                      sep = ""))
      unpen.formula <- as.formula(paste("~", names(expanded.X)[force],
                                        sep = ""))

      lambdas.and.cvls <- profL1(Y[, 1], pen.formula, unpen.formula,
                                 data = expanded.X,
                                 steps = 10, lambda2 = 0, model = "logistic",
                                 trace = FALSE, fold = 5,
                                 standardize = standardize)[c("lambda", "cvl")]
      lambda1.to.use <- lambdas.and.cvls$lambda[which.max(lambdas.and.cvls$cvl)]

      result$main.model <- penalized(Y[, 1], pen.formula, unpen.formula,
                                     data = expanded.X,
                                     standardize = standardize,
                                     lambda1 = lambda1.to.use,
                                     model = "logistic", trace = FALSE)
      result$preds <- predict(result$main.model,
                              penalized = pen.formula,
                              unpenalized = unpen.formula,
                              data = expanded.newX)
      result$fitted <- predict(result$main.model,
                              penalized = pen.formula,
                              unpenalized = unpen.formula,
                              data = expanded.X)
      return(result)
    }
  }

  ## if(g.method == "randomForest.bin"
  ##    || (g.method == "sl" && "randomForest.bin" %in% g.sl.cands)
  ##    || (do.Q && identical(Q.type, "binary.outcome")
  ##        && (Q.method == "randomForest.bin"
  ##            || (Q.method == "sl" && "randomForest.bin" %in% Q.sl.cands)))) {

  ##   candidate.functions$randomForest.bin <- function(X, Y, newX, force,
  ##                                                    standardize) {

  ##     result <- vector("list", length = 3)
  ##     names(result) <- c("preds", "main.model", "force.model")
  ##     class(result) <- "randomForest.bin.result"
      
  ##     result$main.model <- randomForest(X, factor(Y[, 1]))

  ##     if(missing(force)) {
        
  ##       result$preds <- predict(result$main.model, newX, type = "prob")[, "1"]
  ##       return(result)
  ##     }

  ##     ## Truncate so that you don't get Inf's that crash glm

  ##     pred.probs <- predict(result$main.model, type = "prob")[, "1"]
  ##     pred.probs[pred.probs < 0.001] <- 0.001
  ##     pred.probs[pred.probs > 0.999] <- 0.999

  ##     glm.data <- cbind(Y, X[, force, drop = FALSE],
  ##                       preds = qlogis(pred.probs))

  ##     glm.formula <- paste(names(Y), "~", names(X)[force], "*preds", sep = "")
  ##     result$force.model <- glm(glm.formula, data = glm.data, family = binomial,
  ##                               model = FALSE, y = FALSE)

  ##     ## Don't truncate the predictions on newX
  ##     ## (predict.glm can handle Inf's without crashing)

  ##     result$preds <- predict(result$force.model, type = "response",
  ##                             newdata = cbind(newX[, force, drop = FALSE],
  ##                               preds = qlogis(predict(result$main.model, newX,
  ##                                 type = "prob")[, "1"])))
  ##     return(result)
  ##   }
  ## }

  if(g.method == "main.terms.logistic"
     || (g.method == "sl" && "main.terms.logistic" %in% g.sl.cands)
     || (do.Q && identical(Q.type, "binary.outcome")
         && (Q.method == "main.terms.logistic"
             || (Q.method == "sl" && "main.terms.logistic" %in% Q.sl.cands)))) {

    candidate.functions$main.terms.logistic <- function(X, Y, newX, force,
                                                        standardize) {

      result <- vector("list", length = 4)
      names(result) <- c("preds", "fitted", "main.model", "force.model")
      class(result) <- "main.terms.logistic.result"

      formula.string <- paste(names(Y), "~", paste(names(X), collapse = "+"),
                              sep = "")

      result$main.model <- glm(formula.string, data = cbind(Y, X),
                               family = binomial, model = FALSE, y = FALSE)

      result$preds <- predict(result$main.model, newdata = newX,
                              type = "response")
      result$fitted <- predict(result$main.model, type = "response")
      return(result)
    }
  }

  if(g.method == "rpart.bin"
     || (g.method == "sl" && "rpart.bin" %in% g.sl.cands)
     || (do.Q && identical(Q.type, "binary.outcome")
         && (Q.method == "rpart.bin"
             || (Q.method == "sl" && "rpart.bin" %in% Q.sl.cands)))) {

    candidate.functions$rpart.bin <- function(X, Y, newX, force, standardize) {

      result <- vector("list", length = 4)
      names(result) <- c("preds", "fitted", "main.model", "force.model")
      class(result) <- "rpart.result"

      formula.string <- paste(names(Y), "~", paste(names(X), collapse = "+"),
                              sep = "")
      result$main.model <- rpart(formula.string, method = "class",
                                 data = cbind(Y, X))

      if(missing(force) || result$main.model$frame$var[1] == names(X)[force]) {

        result$preds <- predict(result$main.model, newdata = newX,
                                type = "prob")[, 2]
        result$fitted <- predict(result$main.model, type = "prob")[, 2]
        return(result)
      }
      
      if(nrow(result$main.model$frame) == 1) {

        glm.data <- cbind(Y, X[, force, drop = FALSE])
        glm.formula <- paste(names(Y), "~", names(X)[force], sep = "")
        result$force.model <- glm(glm.formula, data = glm.data,
                                  family = binomial,
                                  model = FALSE, y = FALSE)

        result$preds <- predict(result$force.model, type = "response",
                                newdata = newX[, force, drop = FALSE])
        result$fitted <- predict(result$force.model, type = "response")
        return(result)
      }

      ## Truncate so that you don't get Inf's that crash glm

      pred.probs <- predict(result$main.model, type = "prob")[, 2]
      pred.probs[pred.probs < 0.001] <- 0.001
      pred.probs[pred.probs > 0.999] <- 0.999

      glm.data <- cbind(Y, X[, force, drop = FALSE], preds = qlogis(pred.probs))
      glm.formula <- paste(names(Y), "~", names(X)[force], "*preds", sep = "")
      result$force.model <- glm(glm.formula, data = glm.data, family = binomial,
                                model = FALSE, y = FALSE)

      ## Don't truncate the predictions on newX

      result$preds <- predict(result$force.model, type = "response",
                              newdata = cbind(newX[, force, drop = FALSE],
                                preds = qlogis(predict(result$main.model,
                                  newdata = newX, type = "prob")[, 2])))
      result$fitted <- predict(result$force.model, type = "response")
      return(result)
    }
  }

  ## Add continuous outcome candidate functions

  if(do.Q && identical(Q.type, "continuous.outcome")) {

    if(Q.method == "polymars"
       || (Q.method == "sl" && "polymars" %in% Q.sl.cands)) {

      candidate.functions$polymars <- function(X, Y, newX, force, standardize) {

        result <- vector("list", length = 4)
        names(result) <- c("preds", "fitted", "main.model", "force.model")
        class(result) <- "polymars.result"

        expanded.X <- expand.factors(X)
        expanded.newX <- expand.factors(newX)
        
        startmodel <- c(force, NA, NA, NA, 1)
        dim(startmodel) <- c(1, 5)

        result$main.model <- polymars(Y[, 1], expanded.X,
                                      startmodel = startmodel)
        result$preds <- as.vector(predict(result$main.model, x = expanded.newX))
        result$fitted <- as.vector(predict(result$main.model, x = expanded.X))
        
        return(result)
      }
    }

    if(Q.method == "lars"
       || (Q.method == "sl" && "lars" %in% Q.sl.cands)) {

      candidate.functions$lars <- function(X, Y, newX, force, standardize) {

        result <- vector("list", length = 4)
        names(result) <- c("preds", "fitted", "main.model", "force.model")
        class(result) <- "lars.result"

        cv.grid <- 0:50/50

        expanded.X <- expand.factors(X)
        expanded.newX <- expand.factors(newX)
        
        cv.output <- cv.lars(expanded.X, Y[, 1], type = "lar",
                             normalize = standardize, mode = "fraction",
                             index = cv.grid, plot.it = FALSE)$cv

        fraction.to.use <- cv.grid[which.min(cv.output)]

        result$main.model <- lars(expanded.X, Y[, 1], type = "lar",
                                  normalize = standardize)

        coefs <- coef(result$main.model, mode = "fraction", s = fraction.to.use)

        new.preds <- predict(result$main.model, expanded.newX,
                             mode = "fraction",
                             s = fraction.to.use, type = "fit")$fit
        lars.fitted <- predict(result$main.model, expanded.X,
                                   mode = "fraction", s = fraction.to.use,
                                   type = "fit")$fit
        
        if(coefs[force] != 0) {

          result$fitted <- lars.fitted
          result$preds <- new.preds
          return(result)
        }
        
        if(all(coefs == 0)) {

          lm.data <- cbind(Y, X[, force, drop = FALSE])
          lm.formula <- paste(names(Y), "~", names(X)[force], sep = "")
          result$force.model <- lm(lm.formula, data = lm.data)

          result$fitted <- predict(result$force.model, type = "response")
          result$preds <- predict(result$force.model, type = "response",
                                  newdata = newX[, force, drop = FALSE])

          return(result)
        }

        lm.data <- cbind(Y, X[, force, drop = FALSE],
                         preds = lars.fitted)
        lm.formula <- paste(names(Y), "~", names(X)[force], "*preds", sep = "")
        result$force.model <- lm(lm.formula, data = lm.data)

        result$fitted <- predict(result$force.model, type = "response")
        result$preds <- predict(result$force.model, type = "response",
                                newdata = cbind(newX[, force, drop = FALSE],
                                                preds = new.preds))

        return(result)
      }
    }

    ## if(Q.method == "randomForest.cont"
    ##    || (Q.method == "sl" && "randomForest.cont" %in% Q.sl.cands)) {

    ##   candidate.functions$randomForest.cont <- function(X, Y, newX, force,
    ##                                                     standardize) {

    ##     result <- vector("list", length = 3)
    ##     names(result) <- c("preds", "main.model", "force.model")
    ##     class(result) <- "randomForest.cont.result"

    ##     result$main.model <- randomForest(X, Y[, 1])

    ##     lm.data <- cbind(Y, X[, force, drop = FALSE],
    ##                      preds = predict(result$main.model))
    ##     lm.formula <- paste(names(Y), "~", names(X)[force], "*preds", sep = "")
    ##     result$force.model <- lm(lm.formula, data = lm.data,
    ##                              model = FALSE, y = FALSE)
    ##     result$preds <- predict(result$force.model, type = "response",
    ##                             newdata = cbind(newX[, force, drop = FALSE],
    ##                               preds = predict(result$main.model, newX)))
    ##     return(result)
    ##   }
    ## }

    if(Q.method == "main.terms.linear"
       || (Q.method == "sl" && "main.terms.linear" %in% Q.sl.cands)) {

      candidate.functions$main.terms.linear <- function(X, Y, newX, force,
                                                        standardize) {

        result <- vector("list", length = 4)
        names(result) <- c("preds", "fitted", "main.model", "force.model")
        class(result) <- "main.terms.linear.result"

        formula.string <- paste(names(Y), "~", paste(names(X), collapse = "+"),
                                sep = "")
        result$main.model <- lm(formula.string, data = cbind(Y, X),
                                model = TRUE, y = FALSE)
        result$preds <- predict(result$main.model, newdata = newX,
                                type = "response")
        result$fitted <- predict(result$main.model, type = "response")

        return(result)
      }
    }

    if(Q.method == "penalized.cont"
       || (Q.method == "sl" && "penalized.cont" %in% Q.sl.cands)) {

      candidate.functions$penalized.cont <- function(X, Y, newX, force,
                                                     standardize) {

        result <- vector("list", length = 4)
        names(result) <- c("preds", "fitted", "main.model", "force.model")
        class(result) <- "penalized.cont.result"

        expanded.X <- data.frame(expand.factors(X))
        expanded.newX <- data.frame(expand.factors(newX))

        pen.formula <- as.formula(paste("~", paste(names(expanded.X)[-force],
                                                   collapse = "+"),
                                        sep = ""))
        unpen.formula <- as.formula(paste("~", names(expanded.X)[force],
                                          sep = ""))

        lambdas.and.cvls <- profL1(Y[, 1], pen.formula, unpen.formula,
                                   data = expanded.X,
                                   steps = 10, lambda2 = 0, model = "linear",
                                   trace = FALSE, fold = 5,
                                   standardize = standardize)[c("lambda",
                                                                "cvl")]
        lambda1.to.use <-
          lambdas.and.cvls$lambda[which.max(lambdas.and.cvls$cvl)]

        result$main.model <- penalized(Y[, 1], pen.formula, unpen.formula,
                                       data = expanded.X,
                                       lambda1 = lambda1.to.use,
                                       standardize = standardize,
                                       model = "linear", trace = FALSE)
        result$preds <- predict(result$main.model, penalized = pen.formula,
                                unpenalized = unpen.formula,
                                data = expanded.newX)[, 1]
        result$fitted <- predict(result$main.model, penalized = pen.formula,
                                unpenalized = unpen.formula,
                                data = expanded.X)[, 1]
        return(result)
      }
    }

    if(Q.method == "rpart.cont"
       || (Q.method == "sl" && "rpart.cont" %in% Q.sl.cands)) {

      candidate.functions$rpart.cont <- function(X, Y, newX, force,
                                                 standardize) {

        result <- vector("list", length = 4)
        names(result) <- c("preds", "fitted", "main.model", "force.model")
        class(result) <- "rpart.cont.result"

        formula.string <- paste(names(Y), "~", paste(names(X), collapse = "+"),
                                sep = "")
        result$main.model <- rpart(formula.string, method = "anova",
                                   data = cbind(Y, X))

        if(result$main.model$frame$var[1] == names(X)[force]) {
          result$preds <- predict(result$main.model, newdata = newX,
                                  type = "vector")
          result$fitted <- predict(result$main.model, type = "vector")
          return(result)
        }
        
        if(nrow(result$main.model$frame) == 1) {

          lm.data <- cbind(Y, X[, force, drop = FALSE])
          lm.formula <- paste(names(Y), "~", names(X)[force], sep = "")
          result$force.model <- lm(lm.formula, data = lm.data)
          result$preds <- predict(result$force.model, type = "response",
                                  newdata = newX[, force, drop = FALSE])
          result$fitted <- predict(result$force.model, type = "response")
          return(result)
        }

        lm.data <- cbind(Y, X[, force, drop = FALSE],
                         preds = predict(result$main.model, type = "vector"))
        lm.formula <- paste(names(Y), "~", names(X)[force], "*preds", sep = "")
        result$force.model <- lm(lm.formula, data = lm.data)
        result$preds <-
          predict(result$force.model, type = "response",
                  newdata = cbind(newX[, force, drop = FALSE],
                                  preds = predict(result$main.model,
                                                  newdata = newX,
                                                  type = "vector")))
        result$fitted <- predict(result$force.model, type = "response")
        return(result)
      }
    }
  }

######################### Get Ready for Main Loop ##############################

  ## store the actual parameter estimates in this matrix:

  param.estimates <- matrix(0, nrow = ncol(A), ncol = ncol(Y))
  dimnames(param.estimates) <- list(exposure = names(A), outcome = names(Y))

  ## make a copy of that to store plug-in standard errors

  if(estimator != "G-COMP") {
    stand.errs <- param.estimates
  } else stand.errs <- NA
  
  ## lists for returning the final g and Q models
  ## and arrays to store cross-validated risks for returning to user

  g.final.models <- NULL
  Q.final.models <- NULL
  g.cv.risk.array <- NA
  Q.cv.risk.array <- NA

  if(return.final.models) {
  
    ## create list to store final g models

    if(do.g)
      g.final.models <- vector("list", length = ncol(A))

    if(do.Q) {

      ## create list of lists to store final Q models
  
      Q.final.models <- vector("list", length = ncol(A))
      for(i in 1:ncol(A))
        Q.final.models[[i]] <- vector("list", length = ncol(Y))
    }
  }

  ## prepare for modeling the g's

  if(do.g) {

    ## keep track of time taken for g modeling
    
    g.time <- 0

    if(g.method == "sl") {

      ## keep track of time for doing x-val part of sl for g

      g.sl.time <- 0

      ## and time by candidate

      g.sl.cand.times <- rep(0, times = length(g.sl.cands))
      names(g.sl.cand.times) <- g.sl.cands
      
      ## generate matrix to store x-validation predictions

      g.Z <- matrix(0, nrow = nrow(A), ncol = length(g.sl.cands))
      colnames(g.Z) <- g.sl.cands

      ## array to store x-val results

      g.cv.risk.array <- array(0, dim = c(ncol(A), g.num.splits,
                                          length(g.sl.cands)),
                               dimnames = list(exposure = names(A),
                                               split.num = 1:g.num.splits,
                                               candidate = g.sl.cands))

      ## character vector to store which candidate won for each exposure

      g.winning.cands <- vector("character", length = ncol(A))
      names(g.winning.cands) <- names(A)
    }
  }

  ## prepare for modeling the Q's

  if(do.Q) {

    ## keep track of time taken for Q modeling

    Q.time <- 0
    
    if(Q.method == "sl") {

      ## keep track of time for doing x-val part of sl for Q

      Q.sl.time <- 0

      ## and time by candidate

      Q.sl.cand.times <- rep(0, times = length(Q.sl.cands))
      names(Q.sl.cand.times) <- Q.sl.cands

      ## generate matrix to store x-validation predictions
 
      Q.Z <- matrix(0, nrow = nrow(Y), ncol = length(Q.sl.cands))
      colnames(Q.Z) <- Q.sl.cands
      
      ## array to store x-val results

      Q.cv.risk.array <- array(0, dim = c(ncol(A), ncol(Y), Q.num.splits,
                                          length(Q.sl.cands)),
                               dimnames = list(exposure = names(A),
                                               outcome = names(Y),
                                               split.num = 1:Q.num.splits,
                                               candidate = Q.sl.cands))

      ## character matrix to store winning Q candidates

      Q.winning.cands <- matrix("", nrow = ncol(A), ncol = ncol(Y))
      dimnames(Q.winning.cands) <- dimnames(param.estimates)
    }
  }
  
  if(verbose) cat("Starting main loop over exposures\n")

  for(i in 1:ncol(A)) { ## main exposure loop

    if(verbose) cat("starting on exposure number", i, "of", ncol(A),
                    "in main loop over the exposures\n")

    if(do.g) {

      g.start.time <- Sys.time()
      
      if(g.method == "sl") {

        for(split.num in 1:g.num.splits) {

          g.split.index.list <- randomized.split(nrow(A), g.num.folds)

          for(index in g.split.index.list) {

            if(is.null(W)) {

              X.training <- A[-index, -i, drop = FALSE]
              X.validation <- A[index, -i, drop = FALSE]

            } else if(ncol(A) == 1 || !adjust.for.other.As) {

              X.training <- W[-index, , drop = FALSE]
              X.validation <- W[index, , drop = FALSE]

            } else {

              X.training <- cbind(A[-index, -i, drop = FALSE],
                                  W[-index, , drop = FALSE])
              X.validation <- cbind(A[index, -i, drop = FALSE],
                                    W[index, , drop = FALSE])

            }

            Y.training <- A[-index, i, drop = FALSE]

            g.sl.start.time <- Sys.time()
            
            for(candidate in g.sl.cands) {

              g.cand.start.time <- Sys.time()
              
              g.Z[index, candidate] <-
                tryCatch(candidate.functions[[candidate]](
                                    X.training, Y.training, X.validation,
                                    standardize = standardize)$preds,
                         error = function(e) {

                           warning("the candidate ", candidate,
                                   " threw an error during cross-validation.\n",
                                   "The call was: ", e$call, "\n",
                                   "The error message was: ",
                                   e[["message"]], "\n",
                                   "This candidate will be ignored for this ",
                                   "g model,\nwhich is for exposure ", i, " (",
                                   names(A)[i], ")")
                           as.double(NA)
                         })

              g.cand.end.time <- Sys.time()

              g.sl.cand.times[candidate] <-
                (g.sl.cand.times[candidate]
                 + as.double(difftime(g.cand.end.time, g.cand.start.time,
                                      units = "secs")))
            }

            g.sl.end.time <- Sys.time()

            g.sl.time <- g.sl.time + as.double(difftime(g.sl.end.time,
                                                        g.sl.start.time,
                                                        units = "secs"))

          } ## end fold loop

          g.cv.risk.array[i, split.num, ] <- apply(g.Z, 2, get.log.likelihood,
                                                   A[,i])

        } ## end split loop

        ## select best candidate

        if(g.num.splits == 1) {

          ## subsetting g.cv.risk.array will cause dim to be dropped
        
          max.index <- which.max(g.cv.risk.array[i, , ])

        } else {

          ## g.num.splits is > 1, so average over the splits

          max.index <- which.max(apply(g.cv.risk.array[i, , ], 2, mean))

        }

        if(length(max.index) == 0)
          stop("there were no non-NA results in the ",
               "cross-validation for the g model for exposure ", i, " (",
               names(A)[i], ")")

        g.winning.cands[i] <- g.sl.cands[max.index]

        if(verbose) cat("Winning candidate for exposure", i, "is",
                        g.winning.cands[i], "\n")

##################### calculate g(0,W) using the correct model #################

        ## which model?

        fun.to.use <- g.winning.cands[i]

      } else { ## g.method is not "sl"

        fun.to.use <- g.method
      }

      ## use correct predictors (depending on W and adjust.for.other.As)

      if(is.null(W)) {

        g.predictors <- A[, -i, drop = FALSE]

      } else if(ncol(A) == 1 || !adjust.for.other.As) {

        g.predictors <- W

      } else {

        g.predictors <- cbind(A[ , -i, drop = FALSE], W)

      }

      g.model <- candidate.functions[[fun.to.use]](g.predictors,
                                                   A[, i, drop = FALSE],
                                                   g.predictors,
                                                   standardize = standardize)
      g0W <- 1 - g.model$preds

      if(return.final.models) g.final.models[[i]] <- g.model
    
      rm(g.predictors)

      ## truncate unless truncate is FALSE

      truncation.occurred <- FALSE

      if(!identical(truncate, FALSE)) {

        if(any(g0W < truncate)) {
          truncation.occurred <- TRUE
          g0W[g0W < truncate] <- truncate
        }
      }

      g.end.time <- Sys.time()

      g.time <- g.time + as.double(difftime(g.end.time, g.start.time,
                                            units = "secs"))
      
    } ## end if do.g
      
    ## prepare for Q modeling loop

    if(do.Q) {

      if(is.null(W)) {
        
        Q.predictors <- A
        force <- i

      } else if(ncol(A) == 1 || !adjust.for.other.As) {

        Q.predictors <- cbind(A[, i, drop = FALSE], W)
        force <- 1

      } else {

        Q.predictors <- cbind(A, W)
        force <- i
        
      }

      Q.predictors.0 <- Q.predictors
      Q.predictors.0[, force] <- 0
    }

    if(verbose)
        cat("Starting inner loop over outcomes\n")

    for(j in 1:ncol(Y)) { ## inner loop over outcomes

      if(verbose) cat("Starting on outcome number", j, "of", ncol(Y),
                      "in inner loop\n")

      if(do.Q) {

        Q.start.time <- Sys.time()
      
        if(Q.method == "sl") {

          for(split.num in 1:Q.num.splits) {

            Q.split.index.list <- randomized.split(nrow(Y), Q.num.folds)

            for(index in Q.split.index.list) {

              Y.training <- Y[-index, j, drop = FALSE]

              if(is.null(W)) {

                X.training <- A[-index, , drop = FALSE]
                X.validation <- A[index, , drop = FALSE]

              } else if(ncol(A) == 1 || !adjust.for.other.As) {

                X.training <- cbind(A[-index, i, drop = FALSE],
                                    W[-index, , drop = FALSE])
                X.validation <- cbind(A[index, i, drop = FALSE],
                                      W[index, , drop = FALSE])

              } else {

                X.training <- cbind(A[-index, , drop = FALSE],
                                    W[-index, , drop = FALSE])
                X.validation <- cbind(A[index, , drop = FALSE],
                                      W[index, , drop = FALSE])
              }

              Q.sl.start.time <- Sys.time()
              
              for(candidate in Q.sl.cands) {

                Q.cand.start.time <- Sys.time()
                
                tryCatch(Q.Z[index, candidate] <-
                         candidate.functions[[candidate]](
                                         X.training, Y.training, X.validation,
                                         force, standardize)$preds,
                         error = function(e){

                           warning("the candidate ", candidate,
                                   " threw an error during cross-validation.\n",
                                   "The call was: ", e$call, "\n",
                                   "The error message was: ",
                                   e[["message"]], "\n",
                                   "This candidate will be ignored for this ",
                                   "Q model,\nwhich is for exposure ", i, " (",
                                   names(A)[i], ") and outcome ", j, "(",
                                   names(Y)[j], ")")
                           as.double(NA)
                         })

                Q.cand.end.time <- Sys.time()

                Q.sl.cand.times[candidate] <-
                  (Q.sl.cand.times[candidate]
                   + as.double(difftime(Q.cand.end.time, Q.cand.start.time,
                                        units = "secs")))
              }

              Q.sl.end.time <- Sys.time()

              Q.sl.time <- Q.sl.time + as.double(difftime(Q.sl.end.time,
                                                          Q.sl.start.time,
                                                          units = "secs"))

            } ## end fold loop

            if(identical(Q.type, "binary.outcome")) {

              Q.cv.risk.array[i, j, split.num, ] <- apply(Q.Z, 2,
                                                          get.log.likelihood,
                                                          Y[, j])
            } else {

              Q.cv.risk.array[i, j, split.num, ] <- apply(Q.Z, 2, get.MSE,
                                                          Y[, j])

            }

          } ## end split loop

          ## select best candidate

          if(Q.num.splits == 1) {

            ## subsetting will cause dim to be dropped

            final.cv.results <- Q.cv.risk.array[i, j, , ] 

          } else{

            ## average over the splits

            final.cv.results <- apply(Q.cv.risk.array[i, j, , ], 2, mean)

          }

          if(all(is.na(final.cv.results)))
            stop("there were no non-NA results for the cross-validation\n",
                 "for the Q model for exposure ", i, " and outcome ", j)

          if(identical(Q.type, "binary.outcome")) {

            Q.winning.cands[i, j] <- Q.sl.cands[which.max(final.cv.results)]

          } else {

            Q.winning.cands[i, j] <- Q.sl.cands[which.min(final.cv.results)]

          }

          if(verbose) cat("Winning Q candidate for exposure", i, "and outcome",
                          j, "is", Q.winning.cands[i, j], "\n")

          ## do actual Q models depending on method or winning candidate

          ## which model?

          fun.to.use <- Q.winning.cands[i, j]

        } else { ## Q.method is not "sl"

          fun.to.use <- Q.method
        }
        
        ## correct predictors were already set above

        Q.model <-
          candidate.functions[[fun.to.use]](Q.predictors, Y[, j, drop = FALSE],
                                            Q.predictors.0, force, standardize)
        Q0W <- Q.model$preds

        if(return.final.models)
          Q.final.models[[i]][[j]] <- Q.model

        ## Done with Q modeling, record the time
        
        Q.end.time <- Sys.time()

        Q.time <- Q.time + as.double(difftime(Q.end.time, Q.start.time,
                                              units = "secs"))
        
        Ai.logical <- !as.logical(A[,i])

        indicator.Ai.is.zero <- as.double(Ai.logical)

        if(estimator == "TMLE") {

          ZAW <- indicator.Ai.is.zero / ifelse(Ai.logical, 1 - g.model$preds,
                                               g.model$preds)

          Ybar <- mean(Y[, j])
            
          if(identical(Q.type, "binary.outcome")) {

            eps.hat <- glm(Y[, j] ~ 0 + ZAW,
                           offset = qlogis(Q.model$fitted),
                           family = binomial)$coeff

            QhatEps0W <- plogis(qlogis(Q0W) + eps.hat / g0W)

            EhatY0 <- mean(QhatEps0W)

            param.estimates[i, j] <- EhatY0 - Ybar

            IC <- ((ZAW * (Y[, j] - QhatEps0W) + QhatEps0W - EhatY0)
                   - (Y[, j] - Ybar))
            stand.errs[i, j] <- sd(IC) / sqrt(nrow(Y))

          } else { ## continuous outcome

            eps.hat <- lm(Y[, j] ~ 0 + ZAW, offset = Q.model$fitted)$coeff

            QhatEps0W <- Q0W + eps.hat / g0W

            EhatY0 <- mean(QhatEps0W)

            param.estimates[i, j] <- EhatY0 - Ybar
            IC <- ((ZAW * (Y[, j] - QhatEps0W) + QhatEps0W - EhatY0)
                   - (Y[, j] - Ybar))
            stand.errs[i, j] <- sd(IC) / sqrt(nrow(Y))

          }

        } else if(estimator == "DR-IPCW") {

          ## Calculate the DR-IPCW estimate for exposure i and outcome j

          weights <- indicator.Ai.is.zero / g0W

          final.vector <- ( weights * Y[, j] -
                           (indicator.Ai.is.zero - g0W) * Q0W / g0W
                           - Y[, j] )
          param.estimates[i, j] <- mean(final.vector)
          stand.errs[i, j] <- sd(final.vector) / sqrt(nrow(Y))

        } else if(estimator == "G-COMP") {

          param.estimates[i, j] <- mean(Q0W) - mean(Y[, j])

        } else stop("Internal Error 55")

      } else { ## do.Q is FALSE, use regular IPCW

        final.vector <- (as.double(!as.logical(A[,i])) * Y[,j] / g0W
                         - Y[, j] )
        param.estimates[i, j] <- mean(final.vector)
        stand.errs[i, j] <- sd(final.vector) / sqrt(nrow(Y))

      } ## end else for which do.Q is FALSE
    } ## end inner loop over outcomes
  } ## end main loop over exposures

  if(!do.Q || Q.method != "sl") {
    Q.sl.cands <- NA
    Q.num.folds <- NA
    Q.num.splits <- NA
    Q.winning.cands <- NA
    Q.sl.time <- NA
    Q.sl.cand.times <- NA
  }

  if(!do.Q) {
    Q.method <- NA
    Q.type <- NA
    Q.time <- NA
  }

  if(!do.g || g.method != "sl") {

    g.sl.cands <- NA
    g.num.folds <- NA
    g.num.splits <- NA
    g.winning.cands <- NA
    g.sl.time <- NA
    g.sl.cand.times <- NA
  }

  if(!do.g) {
    g.method <- NA
    truncation.occurred <- NA
    truncate <- NA
    g.time <- NA
  }
  
  if(is.null(W)) {
    W.names <- NA
  } else W.names <- names(W)

  if(ncol(A) == 1)
    adjust.for.other.As <- NA

  ## get end time
  
  end.time <- Sys.time()
  
  return.val <- list(param.estimates = param.estimates,
                     plug.in.stand.errs = stand.errs,
                     call = match.call(),
                     num.exposures = ncol(A),
                     num.outcomes = ncol(Y),
                     W.names = W.names,
                     estimator = estimator,
                     g.method = g.method,
                     g.sl.cands = g.sl.cands,
                     g.winning.cands = g.winning.cands,
                     g.cv.risk.array = g.cv.risk.array,
                     g.final.models = g.final.models,
                     g.num.folds = g.num.folds,
                     g.num.splits = g.num.splits,
                     Q.method = Q.method,
                     Q.sl.cands = Q.sl.cands,
                     Q.winning.cands = Q.winning.cands,
                     Q.cv.risk.array = Q.cv.risk.array,
                     Q.final.models = Q.final.models,
                     Q.num.folds = Q.num.folds,
                     Q.num.splits = Q.num.splits,
                     Q.type = Q.type,
                     adjust.for.other.As = adjust.for.other.As,
                     truncate = truncate,
                     truncation.occurred = truncation.occurred,
                     standardize = standardize,
                     boot.param.array = NULL,
                     main.time = as.double(difftime(end.time, start.time,
                                                    units = "secs")),
                     g.time = g.time,
                     Q.time = Q.time,
                     g.sl.time = g.sl.time,
                     Q.sl.time = Q.sl.time,
                     g.sl.cand.times = g.sl.cand.times,
                     Q.sl.cand.times = Q.sl.cand.times)

  class(return.val) <- "multiPIM"

  return(return.val)

} ## end multiPIM function
