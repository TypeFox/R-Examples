###############################################################################
##
## $Id: matchit.R 389 2007-01-10 04:28:44Z enos $
##
## Functions have been prepended with "." and need to be called with
## namespace prefix since they're not exported.
##
################################################################################

.matchit <- function(f,
                     data,
                     treat.var = NULL,
                     exact     = character(),
                     method    = "greedy",
                     n.matches = 1,
                     verbose   = FALSE
                     ){ 

  ## "f is a formula such as y ~ x + z. "data" is a data frame
  ## containing values for the formula. "verbose" is a logical value
  ## specifying whether extra information should be output.

  stopifnot(
            is.logical(verbose),
            is.data.frame(data),
            is.character(exact),
            all(exact %in% names(data))
            )

  stopifnot(method %in% c("greedy","sample","random"))

  if(method %in% "greedy" && n.matches > 1){
    stop("Greedy matching cannot provide more than 1 match. Set n.matches = 1.")
  }
  
  ## Compute treat.var and covarites from the formula 'f', if one is
  ## supplied.  At the very least, we must be armed with a treat.var.
  
  covariates <- NULL
  if(!missing(f) && !is.null(f)){

    stopifnot(is(f, "formula"),
              all(all.vars(f) %in% names(data)))

    if(length(f) != 3){
      stop("f must be a two-sided formula")
    }
    
    if(!is.null(treat.var)){
      warning("treat.var parameter overridden by formula f")
    }
    
    treat.var  <- all.vars(getResponseFormula(f))
    covariates <- all.vars(getCovariateFormula(f))

    if(isTRUE(all.equal(method, "random"))){
      warning("Covariates in formula ignored when using random matching")
    }
  }
  else if(method %in% c("greedy", "sample")){
    stop("Formula required for methods 'greedy' and 'sample'")
  }
  else{
    stopifnot(!is.null(treat.var))
  }

  ## No missing data in any of the required columns is allowed.
  ##
  ## Note that we used to allow missing data here, but we were only
  ## reporting to use via cat/warning that a certain number were
  ## removed.  A result object is required for keeping track of
  ## removed NA's in a systematic way, so we will add back NA handling
  ## when we include one, possibly with an 'na.rm' parameter.

  if(any(is.na(data[c(treat.var, covariates, exact)]))){
    stop("NA's not allowed in formula term, 'treat.var', and 'exact' data columns")
  }

  ## Prep the data by converting character vectors to factors, etc.,
  ## and checking some properties.  May stop.
  
  x <- .data.prep(treat.var, covariates, data)

  ## Calulate propensity score if necessary.
  
  if(method %in% c("greedy", "sample")){
    x$ps <- .ps(f, x)
  }
  
  if(length(exact) > 0){

    res <- do.call(rbind, lapply(split(x, do.call(paste, as.list(x[exact]))),
                                 function(x){
                                   .do.matching(method, x, treat.var, n.matches)
                                 }
                                 )
                   )
  }
  else{
    res <- .do.matching(method, x, treat.var, n.matches)
  }
  
  ## Returns a matrix of matches.  The row names are the names of the
  ## treated observations, and each column contains a set of
  ## matches. 

  res
}

## .do.matching provides some basic checks on the input to the various
## matching sub-functions.

.do.matching <- function(method, x, treat.var, n.matches){

  if(sum(x[[treat.var]]) == 0){
    res <- NULL
  }
  else if(sum(!x[[treat.var]]) == 0){
    stop("No controls found for treatment")
  }
  else{
    
    method.fun <- paste(".", method, ".matching", sep = "")
    res <- do.call(method.fun, list(x = x,
                                    treat.var = treat.var,
                                    n.matches = n.matches))
  }

  res
}

.greedy.matching <- function(x, treat.var, n.matches){

  ## "x" is a data frame containing a vector of propensity scores,
  ## "ps", and a logical vector, "treat.var", specifying whether or
  ## not a unit received treatment.  Returns "x" with an additional
  ## column appended, "matches", which specifies the row name that an
  ## observation has been matched to.

  stopifnot("ps" %in% names(x))

  ## creates a vector of NAs to keep track of all the matches made

  x$matches <- NA

  ## loop continues until all treatments have been matched.  When loop
  ## exits, the "matches" column of "x" will be be populated for
  ## units that have been matched.

  for(i in 1:sum(x[[treat.var]])){

    ## The order in which the ".next.treatment" function returns
    ## treatment units affects the result!  The current implementation
    ## is to choose the treatment unit with the largest propensity
    ## score that has not yet been matched.  Because the algorithm is
    ## "greedy", assigning matches in a different order may result in a
    ## different result.
    
    a.treatment     <- .next.treatment(x, treat.var)

    ## determines the differences in propensity scores between the
    ## treated unit and all other units

    distances       <- abs(x[a.treatment, "ps"] - x$ps)

    a.best.match    <- .best.match(distances, x, treat.var)
    
    ## Don't use NA as an index vector
    
    if(!is.na(a.best.match)){
      x$matches[a.treatment]  <- row.names(x[a.best.match,])
      x$matches[a.best.match] <- row.names(x[a.treatment,])
    }
  }

  matrix(x[x[[treat.var]], "matches"],
         nrow = sum(x[[treat.var]]),
         ncol = 1,
         dimnames = list(row.names(x[x[[treat.var]],]), "1")
         )

}

.sample.matching <- function(x, treat.var, n.matches){

  ## "x" is a data frame containing a column of propensity scores,
  ## "ps", and a column of logical values named "treat.var" that
  ## indicates whether an observation has reserved treatment.
  ## n.matches is the number of sets of matches to be created.
  ## Returns a i x j matrix where i equals the number of observations
  ## that receive treatment and j equals the number of observations
  ## that don't recieve treatment.  The row names are the row names of
  ## the units that receive treatment, and the cell values are the row
  ## names of the observations to which they have been matched.

  stopifnot("ps" %in% names(x))

  ## creates an sum(x[[treat.var]]) x sum(!x[[treat.var]]) matrix of
  ## differences in propensity scores.

  distances <- abs(outer(x[x[[treat.var]], "ps"], x[!x[[treat.var]], "ps"], "-"))
  dimnames(distances) <- list(row.names(x[x[[treat.var]],]),
                              row.names(x[!x[[treat.var]],]))

  ## calculates the probability matrix

  probs <- .prob.matrix(distances)

  ## creates a matrix to store all the matches

  .sample.matrix(probs, n.matches)

}

.random.matching <- function(x, treat.var, n.matches){

  treated <- row.names(x)[x[[treat.var]]]
  
  res <- matrix(nrow = length(treated),
                ncol = n.matches,
                dimnames = list(treated, 1:n.matches)
                )

  for(i in 1:n.matches){

    ## Sampling here is equal weighted, with replacement.  Note that
    ## treatment obs are not allowed in the match.

    ## sample() doesn't work on a charactor vector of length 1, so
    ## handle separately.

    if(sum(!x[[treat.var]]) == 1){
      res[,i] <- row.names(x)[!x[[treat.var]]]
    }
    else{
      res[,i] <- sample(row.names(x)[!x[[treat.var]]],
                        size = nrow(res),
                        replace = TRUE)
    }
  }

  res

}

.prob.matrix <- function(distances){

  ## distances is a matrix of differences in propensity scores.  The
  ## row names of distance are the IDs of the treated units and the
  ## column names are the IDs of the control units.  Returns a matrix
  ## where the value of each cell is the probability that the match
  ## specified by the row column intersection will be made

  res <- matrix(0,
                nrow = nrow(distances),
                ncol = ncol(distances),
                dimnames = dimnames(distances)
                )

  ## Given a distance between two observations, calculates the
  ## probability that these 2 units would be matched.  Smaller
  ## distances are better so we subtract the distances from 1.  Having
  ## taken this difference, raising to the 20th power gives more
  ## preference for the smaller distances

  res <- apply(distances, 2, function(x){(1 - x)^100})
  res

}

.sample.matrix <- function(probs, n.matches){

  ## "probs" is a matrix of probabilities.  The probabilities must be
  ## non-negative.  "n.matches" is how many samples will be drawn from
  ## the matrix.  Returns a nrow(probs) x n.matches matrix of the
  ## samples drawn.

  ## there must be at least one value with a probability > 0 or
  ## sample() throws an error.  All probabilities must be non-negative

  if(!isTRUE(all(apply(probs, 1, function(x){any(x > 0)})))){
    stop("No control has a probability > 0.", call. = FALSE)
  }
  
  if(!isTRUE(all(sign(probs) >= 0))){
    stop("All probabilities must be non-negative.", call. = FALSE)
  }

  res <- matrix(nrow = nrow(probs),
                         ncol = n.matches,
                         dimnames = list(dimnames(probs)[[1]], 1:n.matches)
                         )

  ## fills in res by sampling the "probs" matrix by row.  
  
  indices <- 1:ncol(probs)

  for(i in 1:nrow(probs)){

    res[i,] <- dimnames(probs)[[2]][sample(indices,
                                           size = n.matches,
                                           replace = TRUE,
                                           prob = probs[i,])]
    
  }

  res

}

.ps <- function(formula, data){

  ## Calculates propensity scores using the glm

  ## method excludes NAs from the calculations, but does not
  ## excise them from the returned, numeric vector

  fitted(glm(formula, data, family = binomial("logit")))

}

.next.treatment <- function(x, treat.var){
  
  ## x is a data frame containing columns "treat.var", "ps", and
  ## "matches".  The next unit to be matched has the greatest value in
  ## "ps" of all observations with a a value of NA in the "matches"
  ## column.  Returns the row index, 1, 2, 3, ... (not row name) of
  ## the treatment to be matched next.  The "matches" column in x has
  ## a value of NA if a unit has not been matched or a value equal to
  ## the row name of the unit that this unit will be matched with.
  
  stopifnot(
            "matches" %in% names(x),
            "ps" %in% names(x)
            )

  x[!is.na(x[["matches"]]) | !x[[treat.var]], "ps"] <- NA
  
  ## if there are observations with the same value, return the first

  which.max(x$ps)
}

.best.match <- function(distances, x, treat.var){
  
  ## "distances" is a numeric vector of differences in propensity
  ## scores. "x" is a data frame containing at least three columns,
  ## "treat.var", "ps", and "matches".  "treat.var" is the column name
  ## of the logical vector indicating whether or not a unit receives
  ## treatment.  Returns the row index of the treatment with the least
  ## distance.

  ## sets distance to NA for observations which have already been
  ## matched and are not treatments

  distances[which(x[[treat.var]] | 
                  !is.na(x[["matches"]]))] <- NA

  ## returns NA if there are no eligible controls to match

  if(isTRUE(all(!is.na(x[!x[[treat.var]], "matches"])))){
    return(NA)
  }

  which.min(distances)
}

.data.prep <- function(treat.var, covariates = NULL, data){

  ## Checks that user-supplied data is in proper form, corrects it if
  ## necessary, then returns it "treat.var" is a column of logical or
  ## binary values in "data" that specifies which units received
  ## treatment.  Returns data with the values in "treat.var" converted
  ## to logical values if they were binary, otherwise returns "data"
  ## as it was.  Stops if the treatment column is neither logical or
  ## binary.

  if(!is.logical(data[[treat.var]])){
    stop("Response variable must be logical.", call. = FALSE)
  }

  if(sum(data[[treat.var]]) == 0){
    stop("There must be at least one treated value")
  }

  if(sum(!data[[treat.var]]) == 0){
    stop("There must be at least one non-treated value")
  }
    
  
  ## Any character covariates must be stored as factors for "glm" to
  ## work properly.  If covariates = NULL, nothing is done.
  
  for(i in covariates){
    if(identical(class(data[[i]]), "character")){
      data[[i]] <- as.factor(data[[i]])
    }
  }

  data

}

.remove.nas <- function(data, treat.var, covariates = NULL, verbose = FALSE){
  
  ## "data" is a data frame containing columns for "treat.var" and the
  ## "covariates". "verbose" is a logical scalar and tells the user
  ## where the NA data is in "data".  Returns a data frame with rows
  ## removed that have an NA in the "treat.var" or "covariates"
  ## columns

  omit <- FALSE

  cols <- treat.var
  if(!is.null(covariates)){
    cols <- c(treat.var, covariates)
  }
  
  for(i in cols){

    x <- data[[i]]
    x <- is.na(x)
    omit <- omit | x
  }
  
  res <- data[!omit, ] 
      
  if(verbose){

    ## sum.nas works within a closure to calculate the number of NAs in
    ## each column of a data frame

    sum.nas <- function(x){sum(is.na(x))}

    ## nas is a named, numeric vector where the names are the names of
    ## the columns in "data" and the values are the number of NAs in
    ## each column

    nas <- sapply(data, sum.nas)
    
    ## loops through the different columns of "data" and prints if the
    ## column is referenced by a.formula, "f", and has at least 1 NA

    for(i in names(nas)){

      if(isTRUE(i %in% cols)){
        cat(i, "contains", nas[i], "NA(s)\n")
      }
    }

    cat("\nOriginal data has", nrow(data), "observation(s).",
        nrow(data[data[[treat.var]],]), "receive(s) treatment.\n")

    cat("Treated data has", nrow(res), "observation(s).",
        nrow(res[res[[treat.var]],]), "receive(s) treatment.\n")
  }

  res

}
