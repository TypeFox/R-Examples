vec2vgc <- function (x, steps=200, stepsize=NA, m.max=0)
{
  x <- as.vector(x)    # token vector representing a random or observed sample
  .N <- length(x)      # size of this sample

  if (! (is.numeric(m.max) && length(m.max) == 1 && 0 <= m.max && m.max <= 9))
    stop("'m.max' must be an integer in the range 1 ... 9")
  if (!missing(steps) && !missing(stepsize)) 
    stop("please specify either 'steps' or 'stepsize', but not both")

  if (!missing(stepsize)) {
    .N.steps <- rev(seq(.N, 1, -stepsize)) # make sure last entry is for full sample
  }
  else {
    # (more or less) equally spaced integers, first step at N=1, last step at N=full sample
    .N.steps <- floor(0:(steps-1) * (.N - .5) / (steps-1)) + 1
    .N.steps <- unique(.N.steps)        # avoid duplicates if steps >= .N (or close to it)
  }
  
  idx.first <- !(duplicated(x))     # idx.first marks first occurrences of types
  .V <- cumsum(idx.first)[.N.steps] # vocabulary size = number of first occurrences up to this point
  .Vm.list <- list()
  if (m.max > 0) {
    idx.m <- idx.first
    for (m in 1:m.max) {
      
      ## idx.m marks m-th occurrences, idx.m.1 marks (m+1)-th occurrences
      x[idx.m] <- NA # 1st, 2nd, ..., m-th occurrences now set to NA (cumulatively)
      idx.m.1 <- !duplicated(x) # marks (m+1)-th occurrences, except for first NA
      idx.m.1[1] <- FALSE # first element is always first a occurrence, so set to NA in first iteration

      ## every m-th occurrence increases V_m, every (m+1)-th occurrence decreases it
      .Vm <- cumsum(idx.m) - cumsum(idx.m.1)
      .Vm.list <- c(.Vm.list, list(.Vm[.N.steps]))

      idx.m <- idx.m.1                  # for next iteration
    }
  }

  vgc(N=.N.steps, V=.V, Vm=.Vm.list)
}


## usage example:
## x <- readLines(choose.file())
## vec2spc(x)
## vec2vgc(x, m.max=5)

## for whitespace-spearated tokens (you don't want to do that):
## x <- scan(file, what=character(0), quote=""))
