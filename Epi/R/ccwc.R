ccwc <-
function(entry=0, exit, fail, origin=0, controls=1, match=list(),
         include=list(), data=NULL, silent=FALSE)
{
# Check arguments
  entry <- eval(substitute(entry), data)
  exit <- eval(substitute(exit), data)
  fail <- eval(substitute(fail), data)
  origin <- eval(substitute(origin), data)

  n <- length(fail)
  if (length(exit)!=n)
    stop("All vectors must have same length")
  if (length(entry)!=1 && length(entry)!=n)
    stop("All vectors must have same length")
  if (length(origin)==1) {
    origin <- rep(origin, n)
  }
  else {
    if (length(origin)!=n)
      stop("All vectors must have same length")
  }
# Transform times to correct scale
  t.entry <- as.numeric(entry - origin)
  t.exit <- as.numeric(exit - origin)
# match= argument
  marg <- substitute(match)
  if (mode(marg)=="name") {
    match <- list(eval(marg, data))
    names(match) <- as.character(marg)
  }
  else if (mode(marg)=="call" && marg[[1]]=="list") {
    mnames <- names(marg)
    nm <- length(marg)
    if (!is.null(mnames)) {
      if (nm>1) {
        for (i in 2:nm) {
          if (mode(marg[[i]])=="name")
            mnames[i] <- as.character(marg[[i]])
          else
            stop("illegal argument (match)")
        }
      }
      else {
        for (i in 2:nm) {
          if (mode(marg[[i]])=="name")
            mnames[i] <- as.character(marg[[i]])
          else
            stop("illegal argument (match)")
        }
        mnames[1] <= ""
      }
    }
    names(marg) <- mnames
    match <- eval(marg, data)
  }
  else {
    stop("illegal argument (match)")
  }
  m <- length(match)
  mnames <- names(match)
  if (m>0) {
    for (i in 1:m) {
      if (length(match[[i]])!=n) {
        stop("incorrect length for matching variable")
      }
    }
  }
                                        # include= argument
  iarg <- substitute(include)
  if (mode(iarg)=="name") {
    include <- list(eval(iarg, data))
    names(include) <- as.character(iarg)
  }
  else if (mode(iarg)=="call" && iarg[[1]]=="list") {
    ni <- length(iarg)
    inames <- names(iarg)
    if (ni>1) {
      if (!is.null(inames)) {
        for (i in 2:ni) {
          if (mode(iarg[[i]])=="name")
            inames[i] <- as.character(iarg[[i]])
          else
            stop("illegal argument (include)")
        }
      }
      else {
        for (i in 2:ni) {
          if (mode(iarg[[i]])=="name")
            inames[i] <- as.character(iarg[[i]])
          else
            stop("illegal argument (include)")
        }
        inames[1] <= ""
      }
    }
    names(iarg) <- inames
    include <- eval(iarg, data)
  }
  else {
    stop("illegal argument (include)")
  }
  ni <- length(include)
  inames <- names(include)
  if (ni>0) {
    for (i in 1:ni) {
      if (length(include[[i]])!=n) {
        stop("incorrect length for included variable")
      }
    }
  }
                                        # create group codes using matching variables
  grp <- rep(1,n)
  pd <- 1
  if (m>0) {
    for (im in 1:m) {
      v <- match[[im]]
      if (length(v)!=n)
        stop("All vectors must have same length")
      if (!is.factor(v))
        v <- factor(v)
      grp <- grp + pd*(as.numeric(v) - 1)
      pd <- pd*length(levels(v))
    }
  }
                                        # Create vectors long enough to hold results
  nn <- (1+controls)*sum(fail!=0)
  pr <- numeric(nn)
  sr <- numeric(nn)
  tr <- vector("numeric", nn)
  fr <- numeric(nn)
  nn <- 0
                                        # Sample each group
  if (!silent) {
    cat("\nSampling risk sets: ")
  }
  set <- 0
  nomatch <- 0
  incomplete <- 0
  ties <- FALSE
  fg <- unique(grp[fail!=0])
  for (g in fg) {
                                        # Failure times
    ft <- unique( t.exit[(grp==g) & (fail!=0)] )
                                        # Create case-control sets
    for (tf in ft) {
      if (!silent) {
        cat(".")
      }
      set <- set+1
      case <- (grp==g) & (t.exit==tf) & (fail!=0)
      ncase <- sum(case)
      if (ncase>0)
        ties <- TRUE
      noncase <- (grp==g) & (t.entry<=tf) &
      (t.exit>=tf) & !case
      ncont <- controls*ncase
      if (ncont>sum(noncase)) {
        ncont <- sum(noncase)
        if (ncont>0) incomplete <- incomplete + 1
      }
      if (ncont>0) {
        newnn <- nn+ncase+ncont
        sr[(nn+1):newnn] <- set
        tr[(nn+1):newnn] <- tf
        fr[(nn+1):(nn+ncase)] <- 1
        fr[(nn+ncase+1):newnn] <- 0
        pr[(nn+1):(nn+ncase)] <- (1:n)[case]
        pr[(nn+ncase+1):(newnn)] <-
          sample((1:n)[noncase], size=ncont)
        nn <- newnn
      }
      else {
        nomatch <- nomatch + ncase
      }
    }
  }
  if (!silent) {
    cat("\n")
  }
  res <- vector("list", 4+m+ni)
  if (nn>0) {
    res[[1]] <- sr[1:nn]
    res[[2]] <- map <- pr[1:nn]
    res[[3]] <- tr[1:nn] + origin[map]
    res[[4]] <- fr[1:nn]
  }
  if (m>0) {
    for (i in 1:m) {
      res[[4+i]] <- match[[i]][map]
    }
  }
  if (ni>0) {
    for (i in 1:ni) {
      res[[4+m+i]] <- include[[i]][map]
    }
  }
  names(res) <- c("Set", "Map", "Time", "Fail", mnames, inames)
  if (incomplete>0)
    warning(paste(incomplete, "case-control sets are incomplete"))
  if (nomatch>0)
    warning(paste(nomatch, "cases could not be matched"))
  if (ties)
    warning("there were tied failure times")
  data.frame(res)
}
