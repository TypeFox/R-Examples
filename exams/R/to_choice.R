num_to_schoice <- function(
  correct,                       ## correct numeric solution
  wrong = NULL,                  ## optional wrong numerical solutions
  range = c(0.5, 1.5) * correct, ## range of random wrong solutions
  delta = 1,                     ## minimal distance between solutions
  digits = 2,                    ## digits that should be displayed
  method = c("runif", "delta"),  ## method for generating random results
  sign = FALSE,                  ## randomly change sign?
  verbose = getOption("num_to_choice_warnings") ## display warnings?
)
{
  ## verbosity
  verbose <- if(is.null(verbose)) TRUE else isTRUE(verbose)

  ## make sure that correct solution is in range
  range <- range(c(correct, range))
  
  ## make sure range encompasses wrong solutions as well (if any)
  ## and use a maximum of two wrong solutions
  if(!is.null(wrong)) {
    wrong <- as.numeric(na.omit(wrong))
    range <- range(c(wrong, range))
    if(length(wrong) > 2) wrong <- sample(wrong, 2)
  }

  ## sanity checks of arguments
  if((length(wrong) > 0 && any(abs(round2(wrong, digits = digits) - round2(correct, digits = digits)) < delta)) |
     (length(wrong) > 1 && min(abs(dist(round2(wrong, digits = digits)))) < delta)) {
    if(verbose) warning("specified 'wrong' is too small for 'delta'")
    return(NULL)
  }
  if(diff(range) < delta | max(abs(range - correct)) < delta) {
    if(verbose) warning("specified 'range' is too small for 'delta'")
    return(NULL)
  }
  
  ## include sign changes among wrong solutions?
  if(isTRUE(sign)) sign <- sample(0:(4 - length(wrong)), 1)
  sign <- as.integer(sign)

  ## draw random wrong solutions
  ok <- FALSE
  nle <- sample(0:4, 1)  
  if(abs(correct - range[1]) < (nle * delta) | abs(correct - range[2]) < ((4 - nle) * delta)){
    if(verbose) warning("specified 'range' is too small for 'delta'")
    return(NULL)  
  }
  while(!ok) {
    rand <- switch(match.arg(method),  
      "runif" = c(runif(nle, range[1], correct), runif(4 - nle, correct, range[2])),
      "delta" = c(sample(seq(round2(correct, digits = digits) - delta, range[1], by = -delta), nle),
        sample(seq(round2(correct, digits = digits) + delta, range[2], by = delta), 4 - nle))
    )
    if(sign == 0L) {
      solution <- c(correct, wrong, sample(rand, 4 - length(wrong)))
    } else {
      solution <- c(correct, wrong, sample(rand, 4 - length(wrong) - sign), sample(-c(correct, wrong, rand), sign))
    }
    solution <- round2(solution, digits = digits)   
    ok <- length(unique(solution)) == 5 & min(abs(dist(solution))) >= delta
  }

  ## return shuffled solutions
  o <- sample(1:5)
  list(
    solutions = c(TRUE, rep(FALSE, 4))[o],
    questions = paste("$", format(solution, nsmall = digits, trim = TRUE)[o], "$", sep = "")
  )
}

matrix_to_mchoice <- function(
  x,              ## correct result matrix
  y = NULL,       ## (optional) vector with (potentially) wrong comparions
  lower = FALSE,  ## only elements on lower triangle?
  name = "a",     ## base name for the matrix elements
  comparisons = c("==", "<", ">", "<=", ">=")
                  ## possible (in)equality symbols
)
{
  ## input matrix
  x <- as.matrix(x)
  d <- dim(x)
  
  ## potentially wrong comparions
  if(is.null(y)) y <- sample(-round(max(abs(x))):round(max(abs(x))), 5, replace = TRUE)
  stopifnot(length(y) == 5)

  ## set up potential index pairs and draw random sample of 5 pairs
  ix <- as.matrix(expand.grid(row = 1:d[1], col = 1:d[2]))
  if(lower) ix <- ix[ix[,1] >= ix[,2], , drop = FALSE]
  ix <- ix[rep(sample(1:nrow(ix)), length.out = 5), , drop = FALSE]
    
  ## randomly use value of matrix or y
  prob <- runif(1, 0.3, 0.8)
  y <- ifelse(sample(c(TRUE, FALSE), 5, replace = TRUE, prob = c(prob, 1 - prob)), x[ix], y)
  
  ## randomly choose (in)equality type
  comp <- comp_latex <- sample(comparisons, 5, replace = TRUE)
  comp_latex[comp == "=="] <- "="
  comp_latex[comp == "<="] <- "\\le"
  comp_latex[comp == ">="] <- "\\ge"
  
  ## questions/solution/explanation generation
  questions <- character(5)
  solutions <- logical(5)
  explanations <- character(5)
  for(i in 1:5) {
    solutions[i] <- eval(parse(text = paste(x[ix][i], comp[i], y[i])))
    questions[i] <- paste("$", name, "_{", ix[i,1], ix[i,2], "} ", comp_latex[i], " ", y[i], "$", sep = "")
    explanations[i] <- paste("$", name, "_{", ix[i,1], ix[i,2], "} = ", x[ix][i], 
      if(solutions[i]) "$" else paste(" \\not", comp_latex[i], " ", y[i], "$", sep = ""), sep = "")
  }
  return(list(
    questions = questions,
    solutions = solutions,
    explanations = explanations))
}

matrix_to_schoice <- function(
  x,              ## correct result matrix
  y = NULL,       ## (optional) vector with (potentially) wrong comparions
  lower = FALSE,  ## only elements on lower triangle?
  name = "a",     ## base name for the matrix elements
  delta = 0.5,    ## minimal distance between solutions
  digits = 0      ## digits that should be displayed
)
{
  ## input matrix
  x <- as.matrix(x)
  d <- dim(x)
  
  ## potentially wrong comparions
  if(is.null(y)) y <- -round(max(abs(x))):round(max(abs(x)))
  stopifnot(length(y) >= 5)

  ## set up potential index pairs and draw random sample of 5 pairs
  ix <- as.matrix(expand.grid(row = 1:d[1], col = 1:d[2]))
  if(lower) ix <- ix[ix[,1] >= ix[,2], , drop = FALSE]
  ix <- ix[sample(1:nrow(ix)), , drop = FALSE]
  ix0 <- ix[1, ]
  ix <- ix[-1, , drop = FALSE]

  ## correct/wrong/random solutions
  correct <- x[ix0[1], ix0[2]]
  wrong <- x[ix]
  wrong <- unique(wrong[wrong != correct])
  random <- unique(y[!(y %in% c(correct, wrong))])
  if(length(c(wrong, random)) < 4) {
    warning("'y' contains too few potentially wrong comparisons")
    return(NULL)
  }
  
  ## set up solution vector with correct and 4 wrong/random solutions
  solution <- correct
  while(length(unique(solution)) != 5 || min(abs(dist(solution))) < delta) {
    if(length(random) < 1) {
      solution <- c(correct, sample(wrong, 4))
    } else if(length(random) == 1) {
      solution <- c(correct, sample(wrong, 3), random)
    } else {
      nw <- sample(max(0, 4 - length(random)):min(3, length(wrong)), 1)
      solution <- c(correct, sample(wrong, nw), sample(random, 4 - nw))
    }
    solution <- round2(solution, digits = digits)
  }
    
  ## return shuffled solutions
  o <- sample(1:5)
  list(
    index = ix0,
    name = paste("$", name, "_{", ix0[1], ix0[2], "}", c("", paste("=", correct)), "$", sep = ""),
    solutions = c(TRUE, rep(FALSE, 4))[o],
    questions = paste("$", format(solution, nsmall = digits)[o], "$", sep = "")
  )
}

det_to_schoice <- function(
  x,              ## correct 2x2 result matrix
  y = NULL,       ## wrong 2x2 result matrix
  range = NULL,   ## (optional) range of determinant
  delta = 0.5,    ## minimal distance between solutions
  digits = 0      ## digits that should be displayed
)
{
  ## input matrix
  x <- as.matrix(x)
  d <- dim(x)
  stopifnot(d[1] == 2 & d[2] == 2)
  
  ## range of determinants
  if(is.null(range)) {
    m <- max(2, max(abs(x)))
    range <- c(-1.2 * m^2, 1.2 * m^2)
  } else {
    range <- range(range)
  }

  ## correct result and typical mistake
  correct <- x[1,1] * x[2,2] - x[2,1] * x[1,2]
  wrong <- x[1,1] * x[2,2] + x[2,1] * x[1,2]

  ## correct determinant of wrong matrix y (if any)
  if(!is.null(y)) wrong <- c(wrong, y[1,1] * y[2,2] - y[2,1] * y[1,2])
  wrong <- wrong[wrong != correct]
  solution <- c(correct, wrong)

  ## random results
  if(isTRUE(all.equal(as.vector(x), round(as.vector(x))))) {
    range <- round(range[1]):round(range[2])
    range <- range[!(range %in% solution)]
    while(length(unique(solution)) != 5 || min(abs(dist(solution))) < delta) {
      solution <- c(correct, wrong, sample(range, 4 - length(wrong)))
      solution <- round(solution, digits = digits)
    }
  } else {
    while(length(unique(solution)) != 5 || min(abs(dist(solution))) < delta) {
      solution <- c(correct, wrong, runif(4 - length(wrong), range[1], range[2]))
      solution <- round(solution, digits = digits)
    }
  }
    
  ## return shuffled solutions
  o <- sample(1:5)
  list(
    solutions = c(TRUE, rep(FALSE, 4))[o],
    questions = paste("$", format(solution, nsmall = digits, trim = TRUE)[o], "$", sep = "")
  )
}
