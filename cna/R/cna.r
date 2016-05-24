cna <- function(x, ordering = NULL, strict = FALSE, con = 1, cov = 1, notcols = NULL, maxstep = 5,
                suff.only = FALSE, what="mac") {
  
  if (inherits(x, "truthTab")) {tt <- x}
  else {tt <- truthTab(x)}
  
  if ((! is.null(notcols)) && notcols == "all")
  { 
    colnames(tt)<-chartr("qwertzuiopasdfghjklyxcvbnmQWERTZUIOPASDFGHJKLYXCVBNM","QWERTZUIOPASDFGHJKLYXCVBNMqwertzuiopasdfghjklyxcvbnm",colnames(tt))
    for (i in 1:nrow(tt))
    {
      for (j in 1:ncol(tt))
      {
        if (tt[i,j]==0) {tt[i,j]<-1}
        else if (tt[i,j]==1) {tt[i,j]<-0}
      }
    }
  }
  else if (! is.null(notcols))
  {
    for (i in 1:length(notcols))
    {
      colnumber<-which( colnames(tt)==notcols[i])
      colnames(tt)[colnumber]<-chartr("qwertzuiopasdfghjklyxcvbnmQWERTZUIOPASDFGHJKLYXCVBNM","QWERTZUIOPASDFGHJKLYXCVBNMqwertzuiopasdfghjklyxcvbnm",colnames(tt)[colnumber])
      for (j in 1:nrow(tt))
      {
        if (tt[j,colnumber]==0) {tt[j,colnumber]<-1}
        else if (tt[j,colnumber]==1) {tt[j,colnumber]<-0}
      }        
    }
  }
#  if (inherits(x, "truthTab")) {x <- tt}
x <-tt
  cl <- match.call()
  
  if (nrow(tt) <= 1)
    stop("Truth table must have at least two rows.")
  if (ncol(tt) < 2 || ncol(tt) > 26)
    stop("Truth table must have between 2 and 26 columns.")
  check.ordering(ordering, tt)
  if (maxstep < 1) suff.only <- FALSE
  
  f <- as.vector(attr(tt, "n"))
  sol <- vector("list", length(tt))
  effect.names <- colnames(tt)
  names(sol) <- effect.names
  tt.logi <- as.data.frame(lapply(tt, "mode<-", "logical"))
  
  for (zname in effect.names){
    z <- tt.logi[[zname]]
    xz <- tt.logi[potential.effects(tt, zname, ordering, strict = strict)]
    if (ncol(xz) == 0) next
    
    # Identify and minimize sufficient conditions
    # -------------------------------------------
    xsuff <- sufficient(xz, z, f, con = con)
    
    # initiate list of minimally sufficient conditions
    min.suff.conditions <- vector("list", ncol(xz))
    consist <- vector("list", ncol(xz))
    
    step <- 1L
    if (nrow(xsuff) > 0) repeat{
      if (step >= ncol(xz)){
        min.suff.conditions[[step]] <- xsuff
        consist[[step]] <- attr(xsuff, "consistency")
        break
      }
      
      # reduce sufficient conditions from previous step by removing one component
      reduced <- lapply(seq_along(xsuff),
                        function(i) xsuff[!is.na(xsuff[i]), -i, drop = FALSE])
      reduced.df <- do.call(Rbind, reduced)
      reduced.unique <- unique(reduced.df)   # Eliminate duplicates
      
      # determine the reduced conditions that are sufficient
      xsuff.new <- sufficient(xz, z, f, cond = reduced.unique, con = con, nms = names(xz))
      
      # identify mimally suffcient conditions among the sufficient conditions from the previous step
      row.origin <- unlist(lapply(xsuff, function(x) which(!is.na(x))),
                           use.names = FALSE)
      combinations.grouped <- split(as.integer(combine(reduced.df)), row.origin)
      is.min.suff <- sapply(combinations.grouped,
                            function(x) !any(x %in% attr(xsuff.new, "which.sufficient")))
      min.suff <- xsuff[is.min.suff, , drop = FALSE]
      
      # add to list
      min.suff.conditions[[step]] <- min.suff
      consist[[step]] <- attr(xsuff, "consistency")[is.min.suff]
      
      if (nrow(xsuff.new) == 0) break
      step <- step + 1L
      xsuff <- xsuff.new
    }
    
    minSuff <- do.call(rbind, min.suff.conditions)
    if (length(minSuff) == 0) next
    vc <- verify.conditions(minSuff, xz)
    msc <- data.frame(condition = label.conditions(minSuff),
                      consistency = unlist(consist),
                      coverage = apply(vc, 2, coverage, z, f),
                      stringsAsFactors = FALSE)
    msc <- msc[order(-with(msc, consistency * coverage), msc$condition), ]
    sol[[zname]] <- list(msc = msc)
    
    # Identify and minimize necessary conditions
    # -------------------------------------------
    
    # initiate list of atomic solution formulas
    
    if (!necessary(minSuff, xz, z, f, cov = cov) || suff.only) next
    
    cc <- vc[z, , drop = FALSE]
    fz <- f[z]
    sum.fz <- sum(fz)
    min.nec.conditions <- cov.min.nec <- cons.min.nec <- list()
    step <- 1L
    repeat{
      CondInd <- t(combn(seq_len(ncol(cc)), step))
      # eliminate conditions that can not be minimally necessary
      for (i in seq_along(min.nec.conditions)){
        elimCond <- contains(CondInd, min.nec.conditions[[i]])
        CondInd <- CondInd[!elimCond, , drop = FALSE]
        if (nrow(CondInd) == 0) break
      }
      cover <- apply(CondInd, 1,
                     function(x) sum(fz[apply(cc[, x, drop = FALSE], 1, any)]) / sum.fz)
      nec <- which(cover >= cov)
      # consistency of necessary conditions                                              
      necConds <- apply(CondInd[nec, , drop = FALSE], 1, function(nci) minSuff[nci, , drop = FALSE])
      cons.necConds <- sapply(necConds,
                              function(nc){
                                logvect <- apply(verify.conditions(nc, xz), 1, any)
                                consistency(logvect, z, f)
                              })
      # check whether necessary conditions are also sufficient
      if (length(nec) && con < 1){
        also.consistent <- cons.necConds >= con
        nec <- nec[also.consistent]
        cons.necConds <- cons.necConds[also.consistent]
      }
      if (length(nec)){
        min.nec.conditions <- c(min.nec.conditions, list(CondInd[nec, , drop = FALSE]))
        cov.min.nec <- c(cov.min.nec, list(cover[nec]))
        cons.min.nec <- c(cons.min.nec, list(cons.necConds))
      }
      if (step >= maxstep || step >= ncol(cc)) break
      step <- step + 1L
    }
    
    if (length(min.nec.conditions) > 0){
      lbls <- label.conditions(minSuff)
      asf0 <- lapply(min.nec.conditions,
                     function(x) apply(x, 1, function(r) paste(sort(lbls[r]), collapse = " + ")))
      sol.frame <- data.frame(condition = unlist(asf0),
                              consistency = unlist(cons.min.nec),
                              coverage = unlist(cov.min.nec),
                              stringsAsFactors = FALSE)
      sol.frame <-
        sol.frame[order(-with(sol.frame, consistency * coverage), sol.frame$condition), ]
      rownames(sol.frame) <- NULL
      sol[[c(zname, "asf")]] <- sol.frame
    }
  }
  
  out <- structure(list(), class = "cna")
  out$call <- cl
  out$x <- x
  out$ordering <- ordering
  out$truthTab <- tt
  #  names(sol) <- toupper(names(sol))
  out$solution <- sol
  out$what <- what
  return(out)
}

# print method for class cna
print.cna <- function(x, what=x$what , digits = 3, nsolutions = 5, 
                      row.names = FALSE, show.cases=FALSE, ...){
  cat("--- Coincidence Analysis (CNA) ---\n")
 # cat("\nFunction call:\n", deparse(x$call), "\n", sep = "")
  
  what <- tolower(what)
  if (what == "all") whatl <- rep(TRUE, 4)
  else whatl <- !is.na(match(c("t", "m", "a", "c"), unlist(strsplit(what, ""))))
  names(whatl) <- c("t", "m", "a", "c")
  
  if (whatl["t"]){
    cat("\nTruth table:\n")
    if (show.cases==TRUE)
      {print(x$truthTab,show.cases=TRUE)}
    else
      {print(x$truthTab)}
  }
  if (!is.null(x$ordering))
    cat("\nCausal ordering",
        if(!is.null(x$call$strict) && eval(x$call$strict)) " (strict)", ":\n",
        do.call(paste, c(lapply(x$ordering, paste, collapse = ", "),
                         sep = " < ")),
        "\n", sep = "")
  else cat("\nFactors:", paste(names(x$truthTab), collapse = ", "), "\n")
  
  if (whatl["m"]){
    msc.df <- msc(x)
    cat("\nMinimally sufficient conditions:\n",
        "--------------------------------", sep = "")
    if (nrow(msc.df) == 0) cat("\n*none*\n")
    else for (msc1 in split(msc.df, msc.df$outcome)){
      cat("\nOutcome ", msc1$outcome[1], ":\n", sep = "")
      if (short <- ((nsol <- nrow(msc1)) > nsolutions))
        msc1 <- msc1[seq_len(nsolutions), , drop = FALSE]
      print(msc1[c("condition", "consistency", "coverage")],
            digits = digits, row.names = row.names, ...)
      if (short)
        cat(" ... (total no. of conditions: ", nsol, ")\n", sep = "")
    }
  }
  
  if (any(whatl[c("a", "c")]))
    asf.df <- asf(x)
  
  if (whatl["a"]){
    cat("\nAtomic solution formulas:\n",
        "-------------------------", sep = "")
    if (nrow(asf.df) == 0) cat("\n*none*\n")
    else for (asf1 in split(asf.df, asf.df$outcome)){
      cat("\nOutcome ", asf1$outcome[1], ":\n", sep = "")
      if (short <- ((nsol <- nrow(asf1)) > nsolutions))
        asf1 <- asf1[seq_len(nsolutions), , drop = FALSE]
      print(asf1[c("condition", "consistency", "coverage")],
            digits = digits, row.names = row.names, ...)
      if (short)
        cat(" ... (total no. of formulas: ", nsol, ")\n", sep = "")
    }
  }
  
  if (whatl["c"]){
    csf1 <- csf(asfx = asf.df)
    cat("\nComplex solution formulas:\n",
        "--------------------------\n", sep = "")
    if (nrow(csf1) == 0) cat("*none*\n")
    else {
      if (short <- ((nsol <- nrow(csf1)) > nsolutions))
        csf1 <- csf1[seq_len(nsolutions), , drop = FALSE]
      print(csf1, digits = digits, row.names = row.names, ...)
      if (short)
        cat(" ... (total no. of formulas: ", nsol, ")\n", sep = "")
    }
  }  
  
  invisible(x)
}
