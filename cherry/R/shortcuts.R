pickFisher <- function(p, select = seq_along(p), alpha=0.05, silent=FALSE) {

  if (any(is.na(p))) stop("missing values in input p-values")
  rej <- p[select]
  if (any(is.na(rej))) stop("invalid selection or NA in p-values")
  nr <- setdiff(p, rej)
  nr <- nr[nr > min(rej)]
  lrej <- -2*log(sort(rej, decreasing=TRUE))
  lnr <- -2*log(sort(nr, decreasing=TRUE))
  cum.r <- cumsum(lrej)
  cum.nr <- c(0, cumsum(lnr))
  crit.val <- sapply(1:length(cum.r), function(st) {
    max(qchisq(1-alpha, df=2*(0:(length(cum.nr)-1) + st), lower.tail=TRUE) - cum.nr)
  })
  out <- max(c(0,which(cum.r <= crit.val)))
  if (!silent) {
    cat(length(rej), " hypotheses selected. At confidence level ", 1-alpha, ":\n", sep="")
    cat("False null-hypotheses >= ", length(rej)-out, "; ", sep="")
    cat("True null-hypotheses <= ", out, ".\n", sep="")
    invisible(length(rej)-out)
  } else
    length(rej)-out
}


curveFisher <- function(p, select = seq_along(p), order, alpha=0.05, plot = TRUE) {

  if (any(is.na(p))) stop("missing values in input p-values")
  selected <- !missing(select) || missing(order)
  ordered <- !missing(order)
  if (ordered & selected) 
    stop("please provide either select or order, but not both")
  if (selected)
    ranks <- sort(rank(p, ties.method="first")[select])
  else
    ranks <- rank(p, ties.method="first")[order]
  if (length(ranks)==0 || any(is.na(ranks)))
    stop("invalid selection or NA in p-values")
  others <- setdiff(length(p):1, ranks)
  lpv <- -2*log(sort(p))
  res <- numeric(length(ranks))
  chisqs <- qchisq(1-alpha, df=2*1:length(lpv), lower.tail=TRUE)
  st <- 1
  for (ed in 1:length(ranks)) {
    if (selected)
      ins <- ranks[seq(st,ed)]
    else
      ins <- sort(ranks[1:ed])[seq(st,ed)]
    outs <- setdiff(length(lpv):min(ins), ins)
    cr.v <- max(chisqs[ed-st+1+0:length(outs)] - cumsum(c(0,lpv[outs])))
    rej <- (sum(lpv[ins]) >= cr.v)
    if (rej)
      st <- st+1
    res[ed] <- st-1
  }
  names(res) <- names(ranks)
  if (plot) {
    false <- c(0, res)
    xs <- 1:length(false)-.5
    tots <- 0:length(res)
    plot(xs, tots, type="S", xlab="number of hypotheses", ylab="number of false null-hypotheses", lty=2)
    lines(xs, false, type="S")
    legend("topleft", c(paste("false null-hypotheses (", 100*(1-alpha), "% conf.)", sep=""),"others"), lty=1:2)
    invisible(res)
  } else
    res

}


setClass("hommel",
         representation(
           pvalues = "numeric",         # store vector of p-values (possibly with names)
           jvalues = "numeric",         # stores vector of jvalues
           jumpalpha = "numeric",       # stores vector of alphas where j(alpha) jumps
           adjusted = "numeric",        # stores adjusted p-values 
           simes = "logical"            # stores whether adjusted p-values ar calculated based on a Simes test (if TRUE) or on Hommel's test from 1983 (if FALSE)
         )
)

setGeneric("pvalue", function(object, ...) standardGeneric("pvalue"))
setMethod("pvalue", "hommel", function(object, indicator) {
  if(missing(indicator))
    indicator <- 1: length(object@adjusted)
  object@adjusted[indicator]
})


#own ceil function, needed because small rounding errors can be influential 
ceil <- function(x)
{
  y <- ceiling(x)
  ifelse(isTRUE(all.equal(y-1,x)),y-1,y)
}

#own smaller function, needed because small rounding errors can be influential 
smaller <- function(x,y)
{
  x<y && !(isTRUE(all.equal(x,y)))
}


#given matrix:  find middle column, find minimum, store minimum in array, split in two matrices. 
#               for each sub matrix: give left, right, lower, upper bound. 

#recursive procedure, stops when no splits are possible anymore

#function that finds minimum in certain column
findMin <- function(pvalues, col, lower, upper)
{
  minrow <- max(col,lower) #only take pvalues larger than col and lower
  pval <- pvalues[minrow:upper]
  values <- pval/(minrow:upper - col + 1)
  
  location <- which.min(values)
  min <- values[location]
  location <- minrow + location -1 #true location in full matrix
  
  return(list(min=min,location=location))
}


split <- function(pvalues, left,right,lower, upper)
{
  if(left > right)
  {
    numeric(0) # empty list
  } 
  else
  {
    #find middle column of sub matrix
    mid <- left + floor((right-left)/2) 
    #find minimum of this column + location
    result <- findMin(pvalues, mid, lower, upper)
    
    c(split(pvalues, left, mid-1, lower, result$location), 
      result$min,
      split(pvalues, mid+1, right, result$location, upper))
  }
}


# calculates jumps in function j(alpha)
jumps <- function(alphas) {
  
  n <- length(alphas)
  jumpalpha <- numeric(n) # alpha
  jvalues <- numeric(n) # j(alpha)
  
  index <- 1
  jumpalpha[index] <- alphas[1]
  jvalues[index] <- n
  
  if(n >= 2)
  {
    for(i in 2:n) {
      #if(alphas[i] > jumpalpha[index]) #new jump
      if(smaller(jumpalpha[index],alphas[i])) #new jump
      {
        index <- index+1
        jumpalpha[index] <- alphas[i]
        jvalues[index] <- n+1-i
      }
    }
  }
  
  list(jumpalpha=jumpalpha[1:index], jvalues=jvalues[1:index])
}

# now calculate smallest alphas for each pvalue by looping simultaneously over jumps and pvalues
adjpvalues <- function(pvalues, jumpalpha, jvalues)
{
  n <- length(pvalues)
  adjusted <- numeric(n)
  
  j <- 1
  for(i in 1:n)
  {
    while(j <= length(jvalues) && smaller(jumpalpha[j], pvalues[i] * jvalues[j]))
      j <- j + 1
    
    if(j <= length(jvalues))
      adjusted[i] <- pvalues[i] * jvalues[j]
    if(j > 1) #make sure alpha is in the right interval by taking the maximum over alpha found and the minimal alpha of the corresponding interval
      adjusted[i] <- max(adjusted[i], jumpalpha[j-1]) 
    #NB: if j > length(jvalues), i.e. j=length(jvalues) + 1, adjusted[i] = jumpalpha[length(jvalues)], i.e. the alpha on which j(alpha) goes to 0
    
  }
  adjusted
}

#hommel = F means standard Simes test, hommel = T means test of Hommel 1983
#hommel = F means standard Simes test, hommel = T means test of Hommel 1983
hommelFast <- function(pvalues, simes = TRUE)
{
  if(any(is.na(pvalues))) 
    stop("missing values in input p-values")
  
  if(length(pvalues)==0)
    stop("the pvalues vector is empty")    
  
  #save names in right order
  names <- names(pvalues)
  
  perm <- order(pvalues)
  pvalues <- pvalues[perm]
  
  n <- length(pvalues)
  alphas <- split(pvalues, 1, n, 1, n) #full matrix as input
  # find true values for minima by multiplying with the right factor. this should depend on hommel1983 or simes
  sums <- rep(1,n)
  if(!simes) #extra factor needed
  {
    if(n >= 2)
    {
      for(i in 2:n)
        sums[i] <- sums[i-1] + 1/i      
    }
  }  
  
  alphas <- alphas*(n:1)*rev(sums) 
  result <- jumps(alphas)
  
  jumpalpha <- result$jumpalpha
  jvalues <- result$jvalues
  
  adjusted <- adjpvalues(pvalues, jumpalpha, jvalues)
  adjusted[perm] <- adjusted
  pvalues[perm] <- pvalues #names are still in permuted order here
  
  names(pvalues) <- names
  names(adjusted) <- names  
  
  out <- new("hommel",
             pvalues = pvalues,
             jvalues = jvalues,
             jumpalpha = jumpalpha,
             adjusted = adjusted,
             simes = simes)
  
  return(out)
  
}

findRejections <- function(cats, j) {
  
  n <- length(cats)
  #maximum value of category that can ever lead to a rejection
  maxcats <- min(n,max(cats))
  if(j>0)
    maxcats <- min(maxcats,j)
  
  numRej <- rep(0,n)
  
  parent <- integer(maxcats)
  rank <- integer(maxcats)
  leftmost <- 1:maxcats
  
  find <- function(cat) {
    if(parent[cat] != 0)
      parent[cat] <<- find(parent[cat])
    else
      cat
  }
  
  merge <- function(root1, root2) {
    if(rank[root1] > rank[root2])
      merge(root2, root1)
    
    if(rank[root1] <= rank[root2])
    {
      parent[root1] <<- root2
      leftmost[root2] <<- min(leftmost[root1], leftmost[root2])
      if(rank[root1] == rank[root2])
        rank[root2] <<- rank[root2]+1
    }
  }
  
  for(i in 1:n) {
    numRej[i] <- ifelse(i > 1, numRej[i-1], 0)
    
    if(cats[i]<= maxcats) 
    {
      root1 <- find(cats[i])
      if(leftmost[root1] == 1)
        numRej[i] <- numRej[i]+1
      else
      {
        root2 <- find(leftmost[root1] - 1)
        merge(root1, root2)
      }      
    }
  }
  
  numRej
}

makeCats <- function(jumpalpha, jvalues, alpha, pvalues, simes)
{
  #find j for this alpha, can be made faster by using binary search for example if ever necessary
  i=1
  while(i <= length(jumpalpha) && !smaller(alpha, jumpalpha[i])) 
    i <- i + 1
  
  #if i = length(jumpalpha) + 1, j must be 0. Elseway, j = jvalues[i]
  j <- ifelse(i == (length(jumpalpha) + 1), 0, jvalues[i]) 
  
  #make categories, based on j
  cats <- rep(j+1,length(pvalues)) # j+1 is highest possible value
  
  if(j == 0) #everything could get rejected
    cats <- rep(1, length(pvalues))
  else
  {
    if(!simes)
      scaling <- sum(1/(1:j))
    
    for(i in 1:length(pvalues))
    {
      if(!simes)
        r <- ceil(pvalues[i]*j*scaling/alpha)
      else
        r <- ceil(pvalues[i]*j/alpha)
      if(r < cats[i])
        cats[i] <- max(r,1) #to avoid problems with p-values that are exactly zero      
    }    
  }
  
  return(list(cats=cats, j=j))
}



pickSimes <- function(hommel, select, alpha=0.05, silent=FALSE) {
  
  n <- length(hommel@pvalues)
  
  if(missing(select))
    select <- 1:n
  
  pvalues <- hommel@pvalues[select]
  jvalues <- hommel@jvalues
  jumpalpha <- hommel@jumpalpha
  simes <- hommel@simes
  
  res <- makeCats(jumpalpha, jvalues, alpha, pvalues, simes)
  j <- res$j
  cats <- res$cats  
  
  numRej <- findRejections(cats,j)
  corRej <- numRej[length(numRej)] #total number of rejections
  
  if (!silent) {
    cat(length(pvalues), " hypotheses selected. At confidence level ", 1-alpha, ":\n", sep="")
    cat("False null-hypotheses >= ", corRej, "; ", sep="")
    cat("True null-hypotheses <= ", length(pvalues) - corRej, ".\n", sep="")
    invisible(corRej)
  } else
    corRej
}


#TODO: check wheter select really goes well, seems to be an error with the names.. 
curveSimes <- function(hommel, select, order, alpha=0.05, plot = TRUE)
{
  if (!missing(order) & !missing(select)) 
    stop("please provide either select or order, but not both")
  
  #  if(missing(order))
  #    stop("No order specified.")
  
  n <- length(hommel@pvalues)
  
  if(missing(order) && missing(select))
    select = 1:n
  
  if(!missing(select)) #find the order based on increasing p-values
  {
    p <- hommel@pvalues[select]
    perm <- base::order(p, decreasing = FALSE)
    order <- select[perm]
  }
  
  pvalues <- hommel@pvalues[order]
  jvalues <- hommel@jvalues
  jumpalpha <- hommel@jumpalpha
  simes <- hommel@simes
  
  res <- makeCats(jumpalpha, jvalues, alpha, pvalues, simes)
  j <- res$j
  cats <- res$cats  
  
  res <- findRejections(cats,j)
  names(res) <- names(pvalues)
  
  if (plot) {
    false <- c(0, res)
    xs <- 1:length(false)-.5
    tots <- 0:length(res)
    plot(xs, tots, type="S", xlab="number of hypotheses", ylab="number of false null-hypotheses", lty=2)
    lines(xs, false, type="S")
    legend("topleft", c(paste("false null-hypotheses (", 100*(1-alpha), "% conf.)", sep=""),"others"), lty=1:2)
    invisible(res)
  } else
    res
}



pickSimes_old <- function(p, select = seq_along(p), alpha=0.05, hommel=FALSE, silent=FALSE) {
  
  if (any(is.na(p))) stop("missing values in input p-values")
  ranks <- sort(rank(p, ties.method="first")[select])
  p <- sort(p)
  others <- setdiff(length(p):1, ranks)
  st <- 1
  ed <- length(ranks)
  ready <- FALSE
  while (!ready) {
    ins <- seq_along(p) %in% ranks[seq(st,ed)]
    outs <- (!logical(length(p))) & (cummax(ins)==1) & (!ins)
    participate <- numeric(length(p))
    participate[ins] <- 1+sum(outs)
    participate[outs] <- seq_len(sum(outs))
    maxlag <- cumsum(outs)
    rej <- TRUE
    i <- 0
    while (rej && (i <= sum(outs))) {
      bottom.is <- (participate > i)
      K <- sum(bottom.is)
      if (hommel)
        lag <- floor(1:K - p[bottom.is]/(alpha/(K*sum(1/1:K))))
      else
        lag <- floor(1:K - p[bottom.is]/(alpha/K))
      if (any(lag >= 0 & lag >= maxlag[bottom.is] - i & ins[bottom.is]))
        i <- Inf
      else if (any(lag >= 0))
        i <- i + 1 + max(pmin(lag, maxlag[bottom.is] - i))
      else
        rej <- FALSE
    }
    if (rej) {
      st <- st+1
      ready <- st > ed
    } else
      ready <- TRUE
  }
  out <- ed-st+1
  if (!silent) {
    cat("Rejected ", length(ranks), " hypotheses. At confidence level ", 1-alpha, ":\n", sep="")
    cat("Correct rejections >= ", length(ranks)-out, "; ", sep="")
    cat("False rejections <= ", out, ".\n", sep="")
    invisible(length(ranks)-out)
  } else
    length(ranks)-out
}


curveSimes_old <- function(p, select = seq_along(p), order, alpha=0.05, hommel=FALSE, plot = TRUE) {

  if (any(is.na(p))) stop("missing values in input p-values")
  selected <- !missing(select) || missing(order)
  ordered <- !missing(order)
  if (ordered & selected) 
    stop("please provide either select or order, but not both")
  if (selected) {
    ranks <- sort(rank(p, ties.method="first")[select])
    endpoint <- pickSimes_old(p, select, alpha, hommel, silent=TRUE)
  } else {
    ranks <- rank(p, ties.method="first")[order]
    endpoint <- pickSimes_old(p, order, alpha, hommel, silent=TRUE)
  }
  if (length(ranks)==0 || any(is.na(ranks)))
    stop("invalid selection or NA in p-values")
  p <- sort(p)
  others <- setdiff(length(p):1, ranks)
  res <- numeric(length(ranks))
  st <- 1
  for (ed in 1:length(ranks)) {
    if (selected)
      ins <- seq_along(p) %in% ranks[seq(st,ed)]
    else
      ins <- seq_along(p) %in% sort(ranks[1:ed])[seq(st,ed)]
    outs <- (!logical(length(p))) & (cummax(ins)==1) & (!ins)
    participate <- numeric(length(p))
    participate[ins] <- 1+sum(outs)
    participate[outs] <- seq_len(sum(outs))
    maxlag <- cumsum(outs)
    rej <- TRUE
    i <- 0
    while (rej && (i <= sum(outs))) {
      bottom.is <- (participate > i)
      K <- sum(bottom.is)
      if (hommel)
        lag <- floor(1:K - p[bottom.is]/(alpha/(K*sum(1/1:K))))
      else
        lag <- floor(1:K - p[bottom.is]/(alpha/K))
      if (any(lag >= 0 & lag >= maxlag[bottom.is] - i & ins[bottom.is]))
        i <- Inf
      else if (any(lag >= 0))
         i <- i + 1 + max(pmin(lag, maxlag[bottom.is] - i))
      else
        rej <- FALSE
    }
    if (rej)
      st <- st+1
    res[ed] <- st-1
    if (st > endpoint) {
      res[ed:length(ranks)] <- endpoint
      break
    }
  }
  names(res) <- names(p[ranks])
  if (plot) {
    false <- c(0, res)
    xs <- 1:length(false)-.5
    tots <- 0:length(res)
    plot(xs, tots, type="S", xlab="number of rejections", ylab="number of rejections", lty=2)
    lines(xs, false, type="S")
    legend("topleft", c(paste("correct rejections (", 100*(1-alpha), "% conf.)", sep=""),"others"), lty=1:2)
    invisible(res)
  } else
    res
}









