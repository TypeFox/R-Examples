"bounds" <-
function (nums, dup.inds=FALSE, pre.sort=TRUE) {

  if(dup.inds) {
    ## bounds of concetive duplicated numbers 
    s.ind <- which(!duplicated(nums))
    e.ind   <- c(s.ind[-1]-1, length(nums))
    return(cbind(1:length(s.ind),"start"=s.ind,"end"=e.ind,
                 "length"=(e.ind - s.ind + 1)))

  } else {
    if (!is.numeric(nums)) 
        stop("must supply a numeric vector")
    if (length(nums)==0)
      return(nums)

    if(pre.sort) {
      ## should we pre-sort...
      nums <- sort(unique(nums))
    }

    if (length(nums) == 1) {
      bounds <- c(nums, nums, 1)
      names(bounds) <- c("start", "end", "length")
      ## Edit here to return matrix (not vector) 
      ## following dssp bug report from Yun Liu
      ## Fri, Apr 29, 2011
     return( t(as.matrix(bounds)) )
    }

    bounds <- nums[1]
    nums.start <- nums[1]
    diff.i <- 1
    for (i in 2:length(nums)) {
        if ((nums[i] - diff.i) != nums.start) {
            bounds <- c(bounds, nums[i - 1], nums[i])
            nums.start <- nums[i]
            diff.i <- 1
        }
        else {
            diff.i <- diff.i + 1
        }
    }
    bounds <- c(bounds, nums[length(nums)])
    bounds <- matrix(bounds, ncol = 2, byrow = TRUE,
                     dimnames = list(c(1:(length(bounds)/2)),
                       c("start", "end")))
    bounds <- cbind(bounds, length = (bounds[, 2] - bounds[,1]) + 1)
    return(bounds)
  }
}

