############################################################################################
## package 'secr'
## split.traps.R
## last changed 2012 09 14, 2012-12-18 (usage)
############################################################################################

split.traps <- function (x, f, drop = FALSE, prefix='S', byoccasion = FALSE, ...) {
  if (!inherits(x, 'traps'))
      stop ("argument to split.traps should have class 'traps'")
  if (ms(x))
      stop ("'split.traps' is not suitable for multi-session traps")

  options(warn=-1)
  f <- factor(f)

  if (any(!is.na(as.numeric(levels(f))))) {
      f <- factor(paste (prefix,f,sep=''))
      sp <- paste(prefix, levels(polyID(x)), sep='')
  }
  else {
      sp <- levels(polyID(x))
  }
  options(warn=0)

  if (byoccasion) {
      usg <- usage(x)
      usg[] <- as.integer(usg)
      if (is.null(usg))
          stop ("byoccasion requires usage codes")
      if (length(f) != ncol(usg))
          stop ("byoccasion requires length(f) == ncol(usage(x))")
      splitusgt <- split(data.frame(t(usg)), f)
      sessiontraps <- function(u) {
          tmp <- x
          usage(tmp) <- t(u)
          used <- apply(u,2,sum)>0
          tmp <- subset(tmp, used)
          u <- usage(tmp)
          dimnames(usage(tmp)) <- list(1:nrow(u), 1:ncol(u))
          tmp
      }
      out <- lapply(splitusgt, sessiontraps)
      class(out) <- c('list','traps')
      out
  }

  else {

      if (detector(x) %in% .localstuff$polydetectors) {
          if (length(unique(f)) > length(levels(polyID(x))))
              warning ("split factor does not match traps object")
      }
      out <- list()
      for (i in levels(f)) {
          if (detector(x) %in% .localstuff$polydetectors) {
              temp <- subset (x, subset = (sp == i), ...)
          }
          else
              temp <- subset (x, subset = (f == i), ...)
          if (!drop | (nrow(temp)>0))
              out[[i]] <- temp
          if (detector(x) %in% .localstuff$pointdetectors) {
              spacing(out[[i]]) <- spacing(out[[i]], recalculate = TRUE)
          }
      }
      class (out) <- c('list', 'traps')
      out
  }
}

###############################################################################

