################################################################################
##
## $Id: weight.R 386 2007-01-10 04:12:20Z enos $
##
## Form weights from a data.frame and an in.var.
##
################################################################################

weight <- function(x, in.var, type, size, sides, weight.var = NULL, verbose = FALSE){

  ## Does x have the in.var in question?

  stopifnot(in.var %in% names(x))
  stopifnot(is.numeric(x[[in.var]]))

  x <- x[c(in.var, weight.var)]
  
  ## Figure out what type we're working with
  
  type.choices <- c("relative","equal","linear","sigmoid","centroid","complex")
  type         <- type.choices[pmatch(type, type.choices)]
  if(is.na(type)){
    stop("Invalid type")
  }

  sides <- unique(sides)
  if(!(all(sides %in% c("long", "short")) &&
       length(sides > 0))){
    stop("Invalid sides")
  }
  
  if(is.character(size)){

    size.choices <- c("all","decile","quintile","quartile","tercile","demile")
    size         <- size.choices[pmatch(size, size.choices)]
    if(is.na(type)){
      stop("Invalid size")
    }
    if(type == "relative" & !is.numeric(size) & size != "all"){
      stop("size must be set to all for type relative")
    }
    
    ## We base size on non-na in.vars.  Always use the floor for
    ## sizes that are not integers.
    
    num.non.na <- sum(!is.na(x[[in.var]]))
    size       <- floor(switch(size,
                               all      = num.non.na,
                               decile   = num.non.na / 10,
                               quintile = num.non.na / 5,
                               quartile = num.non.na / 4,
                               tercile  = num.non.na / 3,
                               demile   = num.non.na / 2,
                               stop("Invalid size")))
  }

  ## At this point a character value for size will have been converted
  ## to a number.
  
  if(is.numeric(size)){

    ## Is size too small?

    if(type != "relative" && size < 1){
      stop("Size is < 1 per side.  Increase size parameter or number non-na in.var values.")
    }
    
    ## Are there enough values in in.var to go around?

    side.denom <- ifelse(all(c("long","short") %in% sides), 2, 1)
    if(type != "relative" && size > floor(sum(!is.na(x[[in.var]])) / side.denom)){
      stop("Size too large for the number of non-na ranks")
    }
    
  }
  else{
    stop("Invalid size")
  }

  if(verbose){
    cat("Creating", type, "weights using size", size,
        "from", sum(!is.na(x[[in.var]])), "candidates\n")
  }

  if(type != "relative"){
  
    ## NA out those outside of the size range.
    
    x$r <- rank(x[[in.var]], na.last = "keep", ties.method = "first")
    is.na(x[[in.var]]) <- ! (x$r <= size | max(x$r, na.rm = TRUE) - size < x$r)


    ## The strategy for the rest of this function is to form weights as
    ## if we were dealing with both the long and short sides.  A problem
    ## arises, however, if there are fewer than twice the number of
    ## non-na ranks and we want ranks for only one side.  To deal with
    ## this, we add filler entries according to which single side we
    ## want in the end.

    x$orig <- TRUE
    num.to.fill <- size * 2 - sum(!is.na(x$r))

    ## It should never happen that we have to perform a fill and we're
    ## dealing with both long and short.  The num-non-na checks above
    ## will catch this, but let's just be sure.

    if(num.to.fill > 0 && length(sides) == 2){
      stop("Something is wrong: shouldn't be filling if working on long and short sides")
    }
    
    if(num.to.fill > 0){
      fill.val            <- ifelse(sides == "long", -Inf, Inf)
      fill.rows           <- x[0,][1:num.to.fill,]
      fill.rows[[in.var]] <- fill.val
      fill.rows$orig      <- FALSE

      x <- rbind(x, fill.rows)
    }
    
    x$r <- rank(x[[in.var]], na.last = "keep", ties.method = "first")
  }
    
  ## Construct relative weights for each formation type
  

  if(type == "equal"){
    x$rel.weight <- ifelse(x$r <= size, -1, 1)
  }
  else if(type == "linear"){

    ## Subtracting the median rank from the rank gives weights spaced
    ## by 1.  If there are an even number of weights, we must add +/-
    ## 0.5 to each side to make the smallest weight 1 (or, in absolute
    ## weight space, equal to the amount between successive linear
    ## weights).

    x$r <- x$r - median(x$r, na.rm = TRUE)
    x$r <- x$r + sign(x$r) * x$r %% 1

    x$rel.weight <- x$r
  }
  else if(type == "sigmoid"){

    ## The sigmoid calculation below requires linear-spaced weights on
    ## [-1,1].
    
    x$r <- x$r - median(x$r, na.rm = TRUE)
    x$r <- x$r + sign(x$r) * x$r %% 1
    
    x$r <- x$r / max(x$r, na.rm = TRUE)
    x$rel.weight <- sign(x$r) * 1 / (1 + exp( -(abs(x$r) - 0.3) / 0.1))

  }
  else if(type == "centroid"){
    
    n <- sum(!is.na(x[[in.var]]))
    a <- 0.375
    x$rel.weight <- - qnorm((n + 1 - x$r - a) / (n - 2 * a + 1))
    
  }
  else if(type == "complex"){

    ## This piecewise function requires linear-spaced weights on
    ## [-1,1].
    
    x$r <- x$r - median(x$r, na.rm = TRUE)
    x$r <- x$r + sign(x$r) * x$r %% 1
    
    x$r <- x$r / max(x$r, na.rm = TRUE)

    w.complex <- function(x){
          
      if(is.na(x)){
        return(NA)
      }
      
      stopifnot(-1 <= x && x <= 1)
      s <- sign(x)
      x <- abs(x)
      
      if(0 <= x && x < 0.5){
        val <- 6 * x + 1
      }
      else if(0.5 <= x && x < 0.77){
        val <- 77.8 * x - 34.9
      }
      else if(0.77 <= x && x < 0.91){
        val <- 25
      }
      else if(0.91 <= x && x < 0.96){
        val <- 500 * x - 430
      }
      else if(0.96 <= x && x <= 1){
        val <- 50
      }
      else{
        stop("x out of range")
      }
      return(val * s)
    }

    x$rel.weight <- sapply(x$r, w.complex)
    
  }
  else if(type == "relative"){
    x$rel.weight <- x[[in.var]]
  }

  ## Now that we've arrived at a set of relative weights, decompose
  ## any of the filler work we may have done, and remove weights that
  ## aren't in line with our side directives.  Note that we don't do
  ## any of this if we're working with relative weights to start.

  if(type != "relative"){
    x <- x[x$orig,]

    if(length(sides) == 1){
      if(sides == "long"){
        is.na(x$rel.weight) <- x$rel.weight < 0
      }
      else if(sides == "short"){
        is.na(x$rel.weight) <- x$rel.weight > 0
      }
    }
  }
  
  ## Turn relative into absolute weights.
  
  x$abs.weight <- NA
  abs.sum.pos <-
    abs(sum(x$rel.weight[x$rel.weight > 0], na.rm = TRUE))
  abs.sum.neg <-
    abs(sum(x$rel.weight[x$rel.weight < 0], na.rm = TRUE))
  
  x$abs.weight <- ifelse(x$rel.weight > 0, x$rel.weight / abs.sum.pos,
                         x$abs.weight)
  x$abs.weight <- ifelse(x$rel.weight < 0, x$rel.weight / abs.sum.neg,
                         x$abs.weight)

  ## Finally, if weight.var is not null, apply these weights to the
  ## resulting weight vector, much like an override.

  if(!is.null(weight.var) && length(weight.var) == 1){
    stopifnot(weight.var %in% names(x))

    x$abs.weight <- ifelse(!is.na(x[[weight.var]]), x[[weight.var]],
                           x$abs.weight)
  }

  x$abs.weight
}
