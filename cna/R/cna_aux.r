# Auxiliary function for cna
# --------------------------

# function combine: creates a factor of all combinations of x[rows],
# with the order of the levels corresponding to the order as thy occur in the x
combine <- function(x, rows = rownames(x)){
  rows <- do.call(paste, c(x, sep = ""))
  factor(rows, levels = unique(rows))
  }

# function Rbind: modified version of rbind
Rbind <- function(..., fill = NA){
  dots <- list(...)
  nms <- sort(unique(unlist(lapply(dots, names))))
  dots <- dots[sapply(dots, nrow) > 0]
  if (length(dots) == 0){
    out <- as.data.frame(rep(list(integer(0)), length(nms)))
    names(out) <- nms
    return(out)
    }
  for (i in seq_along(dots))
      dots[[i]][setdiff(nms, names(dots[[i]]))] <- fill
  out <- do.call(rbind, dots)
  out[names(out)]
  }


# function label.conditions
label.conditions <- function(cond){
  if (length(cond) == 0) return (NULL)
  condm <- as.matrix(cond)
  ch <- matrix("", ncol = ncol(cond), nrow = nrow(cond))
  # mettre en maj. ou ne rien faire
    ch[!is.na(condm) & condm] <- 
    #  toupper(names(cond))[col(condm)][!is.na(condm) & condm]
    names(cond)[col(condm)][!is.na(condm) & condm]
  # mettre en min. ou mettre en pol. opp.
  ch[!is.na(condm) & !condm] <- 
chartr("qwertzuiopasdfghjklyxcvbnmQWERTZUIOPASDFGHJKLYXCVBNM","QWERTZUIOPASDFGHJKLYXCVBNMqwertzuiopasdfghjklyxcvbnm",names(cond))[col(condm)][!is.na(condm) & !condm]
# tolower(names(cond))[col(condm)][!is.na(condm) & !condm]
  
  apply(ch, 1, function(x) paste(x[nchar(x) > 0], collapse = "*"))
  }

  
# function check.ordering
check.ordering <- function(ordering, tt){
  if (!is.null(ordering)){
    if (!all(colnames(tt) %in% unlist(ordering)))
      stop("Not all factors are included in the ordering.")
    if (anyDuplicated(unlist(ordering)))
      stop("At least one factors appears twice in the ordering.")
    }
  invisible(NULL)            
  }            
  
# function potential.effects
potential.effects <- function(x, zname, ordering, strict = FALSE){
  if (is.null(ordering))
    poteff <- setdiff(names(x), zname)
  else {  
    rownames(x) <- NULL
    ordLevelZ <- 
      which(sapply(ordering, function(x) any(zname %in% x)))
    poteff <- unlist(ordering[seq_len(ordLevelZ - strict)])
    poteff <- setdiff(poteff, zname)
    }
  if (!is.null(poteff)) sort(poteff)
  }

# function verify.conditions
verify.conditions <- function(cond, data){
  eq.or.NA <- function(da, co) ifelse(da == co | is.na(co), TRUE, FALSE)
  logMatList <- lapply(names(cond), function(colnm) outer(data[[colnm]], cond[[colnm]], eq.or.NA))
  matrix(do.call(mapply, c(all, logMatList, SIMPLIFY = TRUE)), nrow(data))
  }              

# function consistency
consistency <- function(selection, z, f) sum(f[selection & z]) / sum(f[selection])

# function sufficient
sufficient <- function(x, z, f, con = 1, cond = NULL, nms = names(x)){
  stopifnot(nrow(x) == length(z), length(z) == length(f), all(names(x) %in% nms))
  if (!is.null(cond)) stopifnot(all(names(cond) %in% names(x)))
  if (is.null(cond))
    cond <- x[which(!duplicated(combine(x))), , drop= FALSE]
  cc <- verify.conditions(cond, x)
  consist <- apply(cc, 2, consistency, z, f)
  isSuff <- consist >= con
  if (!any(isSuff)){
    out <- matrix(integer(0), nrow = 0, ncol = length(nms))
    colnames(out) <- nms
    out <- as.data.frame(out)
    attr(out, "which.sufficient") <- which(isSuff)
    attr(out, "consistency") <- consist[isSuff]
    return(out)
    }
  out <- cond[isSuff, , drop = FALSE]
  if (length(additional.nms <- setdiff(nms, names(out))) > 0)
    out[additional.nms] <- NA_integer_
  out <- out[sort(names(out))]
  attr(out, "which.sufficient") <- which(isSuff)
  attr(out, "consistency") <- consist[isSuff]
  out
  }

# function necessary
coverage <- function(selection, z, f, sumf = sum(f[z]))
  sum(f[selection & z]) / sumf

# function necessary
necessary <- function(cond, x, z, f, cov = 1){
  stopifnot(nrow(x) == length(z), length(z) == length(f), 
            all(names(cond) %in% names(x)))
  cc <- verify.conditions(cond, x)
  anyCond <- apply(cc, 1, any) 
  cover <- coverage(anyCond, z, f, sumf = sum(f[z]))
  out <- (cover >= cov)
  attr(out, "coverage") <- cover
  out
  }

# function reduce.nec.cond
reduce.nec.cond <- function(cond)
  lapply(seq_len(nrow(cond)), function(i) cond[-i, , drop = FALSE])


# function contains: 
# returns a logical vector of length nrow(longer) with nth entry TRUE 
# if there is a row in shorter with all values contained in the nth row of longer
contains <- function(longer, shorter){
  a <- sapply(seq_len(nrow(shorter)), function(i)
    apply(matrix(apply(longer, 1, match, x = shorter[i, ], nomatch = 0) > 0,
                 nrow = ncol(shorter)),
          2, all))
  apply(matrix(a, nrow = nrow(longer)), 1, any)
  }
  
