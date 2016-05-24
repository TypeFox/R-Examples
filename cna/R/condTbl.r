
# extract msf from cna-object
msc <- function(x){
  stopifnot(inherits(x, "cna"))
  msc.list <- lapply(x$solution, "[[", "msc")
  if (all(sapply(msc.list, is.null)))
    return(data.frame(outcome = character(0),
                      condition = character(0),
                      consistency = numeric(0),
                      coverage = numeric(0)))
  all.msc <- do.call(rbind, msc.list)
  outcome <- rep(names(msc.list),
                 sapply(msc.list, function(fr) if(is.null(fr)) 0 else nrow(fr)))
  out <- 
    data.frame(outcome = outcome,
               condition = paste(all.msc$condition, "->", outcome),
               consistency = all.msc$consistency,
               coverage = all.msc$coverage,
               row.names = NULL, stringsAsFactors = FALSE)
  structure(out, class = c("condTbl", "data.frame"))
  }

# extract msf from cna-object
asf <- function(x){
  stopifnot(inherits(x, "cna"))
  asf.list <- lapply(x$solution, "[[", "asf")
  if (all(sapply(asf.list, is.null)))
    return(data.frame(outcome = character(0),
                      condition = character(0),
                      consistency = numeric(0),
                      coverage = numeric(0)))
  all.asf <- do.call(rbind, asf.list)
  outcome <- rep(names(asf.list),
                 sapply(asf.list, function(fr) if(is.null(fr)) 0 else nrow(fr)))
  out <-                
    data.frame(outcome = outcome,
               condition = paste(all.asf$condition, "<->", outcome),
               consistency = all.asf$consistency,
               coverage = all.asf$coverage,
               row.names = NULL, stringsAsFactors = FALSE)
  structure(out, class = c("condTbl", "data.frame"))
  }

# Extract csf from cna object
csf <- function(x, asfx = asf(x)){
  if (nrow(asfx) == 0)
    return(data.frame(solution = character(0),
      consistency = numeric(0), coverage = numeric(0)))
  splitasf <- lapply(split(asfx, asfx$outcome), "[", -1)
  splitasf <- lapply(splitasf, as.data.frame)
  rnms.grid <- do.call(expand.grid,
    c(lapply(splitasf, function(x) seq_len(nrow(x))), KEEP.OUT.ATTRS = FALSE))
  splitasf.expanded <-
    mapply(function(asf, gr) asf[gr, , drop = FALSE],
           splitasf, rnms.grid, SIMPLIFY = FALSE)
  X <- do.call(cbind, splitasf.expanded)
  csfCon <- apply(X[grepl("\\.consistency$", names(X))], 1, min)
  csfCov <- apply(X[grepl("\\.coverage$", names(X))], 1, min)
  csfOrder <- order(csfCon * csfCov, decreasing = TRUE)
  conds <- X[grepl("\\.condition$", names(X))]
  conds[] <- paste("(", unlist(conds), ")", sep = "")
  conds[] <- lapply(conds, format, justify = "right")
  csfName <- apply(conds, 1, paste, collapse = "  *  ")
  out <- data.frame(condition = csfName,
    consistency = csfCon, coverage = csfCov)[csfOrder, ]
  rownames(out) <- NULL
  class(out) <- c("condTbl", "data.frame")
  out
  }
  
# print method for class condTbl
print.condTbl <- function(x, digits = 3, ...){
  arrow <- regexpr("->", x$condition)
  valid <- arrow > 2
  x$condition <- 
    paste(format(mapply(substring, x$condition[valid], 1, arrow[valid] - 2)),
          mapply(substring, x$condition[valid], arrow[valid] - 1), sep = "")
  if (!is.null(x$consistency)) x$consistency <- 
    formatC(x$consistency, format = "f", digits = digits)
  if (!is.null(x$coverage)) x$coverage <- 
    formatC(x$coverage, format = "f", digits = digits)
  print.data.frame(x, ...)
  invisible(NULL)
  }

# as.condTbl: builds a "condTbl"-object from a list of "cond"-objects 
as.condTbl <- function(condlst, ...){
  stopifnot(is.list(condlst), all(sapply(condlst, inherits, "cond")))
  outcome <- sapply(strsplit(names(condlst), "-> "), "[", 2)
  outcome[is.na(outcome)] <- "(No outcome)"
  out <- data.frame(outcome = outcome, condition = names(condlst),
                    stringsAsFactors = FALSE)
  if (any(!is.na(con <- getAttribute(condlst, "consistency", use.names = FALSE))))
    out$consistency <- con
  if (any(!is.na(cov <- getAttribute(condlst, "coverage", use.names = FALSE))))
    out$coverage <- cov
  if (any(!is.na(f <- getAttribute(condlst, "freq", use.names = FALSE))))
    out$rel.freq <- f
  class(out) <- c("condTbl", "data.frame")
  out
  } 

# getAttribute: 
# auxiliary function to extract attributes from a list of "cond"-objects
getAttribute <- function(condlst, attrName, use.names = TRUE){
  val <- lapply(condlst, attr, attrName)
  val[sapply(val, is.null)] <- NA
  if (!use.names) val <- unname(val)
  unlist(val)
  }
  
