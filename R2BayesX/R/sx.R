sx <- function(x, z = NULL, bs = "ps", by = NA, ...)
{
  by <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  call <- match.call()

  available.terms <- c(
    "rw1", "rw2",
    "season",
    "ps", "psplinerw1", "psplinerw2", "pspline",
    "te", "pspline2dimrw2", "te1", "pspline2dimrw1",
    "kr", "kriging",
    "gk", "geokriging",
    "gs", "geospline",
    "mrf", "spatial",
    "bl", "baseline",
    "factor",    
    "ridge", "lasso", "nigmix",
    "re", "ra", "random",
    "cs", "catspecific",
    "offset",
    "generic",
    "rsps", "hrandom_pspline"
  )
  if(!bs %in% available.terms) stop(paste("basis type", sQuote(bs), "not supported by BayesX"))

  if(bs %in% c("rsps", "hrandom_pspline")) {
    bs <- "rsps"
    x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
    rcall <- paste("r(x = ", by, ", bs = ", sQuote(bs), ", by = ", x, ", ...)", sep = "")
    rval <- eval(parse(text = rcall))
  } else {
    if(length(grep("~", term <- deparse(call$x))) && bs %in% c("re", "ra", "random")) {
      x <- deparse(substitute(x), backtick = TRUE, width.cutoff = 500)
      rcall <- paste("r(x = ", x, ", by = ", by, ", ...)", sep = "")
      rval <- eval(parse(text = rcall))
    } else {
      k <- -1
      m <- NA
      xt <- list(...)
      if("m" %in% names(xt))
        stop("argument m is not allowed, please see function s() using this specification!")
      if("k" %in% names(xt))
        stop("argument k is not allowed, please see function s() using this specification!")
      if(!is.null(xt$xt))
        xt <- xt$xt
      warn <- getOption("warn")
      options("warn" = -1)
      if(by != "NA" && is.vector(by) && length(by) < 2L && !is.na(as.numeric(by))) {
        xt["b"] <- by
        by <- "NA"
      }
      options("warn" = warn)
      if(bs %in% c("pspline2dimrw1", "pspline2dimrw2", "te",
        "gs", "geospline", "kr", "gk", "kriging", "geokriging")) {
        if(by != "NA")
          stop(paste("by variables are not allowed for smooths of type bs = '", bs, "'!", sep = ""))
      }
      if(bs == "te1") bs <- "pspline2dimrw1"
      if(!is.null(xt$knots))
        xt$nrknots <- xt$knots
      if(bs %in% c("ps", "te", "psplinerw1", "psplinerw2", "pspline",
        "pspline2dimrw2", "pspline2dimrw1", "gs", "geospline")) {
        if(!is.null(xt$degree))
          m <- xt$degree
        if(!is.null(xt$order)) {
          if(is.na(m))
            m <- c(3L, xt$order)
          else
            m <- c(m[1L], xt$order)
        }
        if(length(m) < 2L && (bs %in% c("gs", "geospline")))
          m <- c(3L, 1L)
        if(length(m) < 2L && is.na(m))
          m <- c(3L, 2L)   
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 2L)
        if(is.null(xt$order) && length(m) < 2L)
          m <- c(m, 1L)
        m[1L] <- m[1L] - 1L
        if(!is.null(xt$nrknots))
          k <- xt$nrknots + m[1L]
        else {
          if(bs %in% c("ps", "psplinerw1", "psplinerw2", "pspline"))
            k <- 20L + m[1L]
          else
            k <- 10L + m[1L]
        }
      }
      if(bs %in% c("kr", "gk", "kriging", "geokriging")) {
        m <- c(1L, 1L)
        if(!is.null(xt$nrknots)) {
          k <- xt$nrknots
        } else k <- -1L
      }
      if(!is.null(xt$map))
        xt$map.name <- as.character(call$map)
      xt[c("degree", "order", "knots", "nrknots")] <- NULL
      if(!length(xt))
        xt <- NULL
      if(!is.null(call$z)) 
        term <- c(term, deparse(call$z))
      rval <- mgcv::s(x, z, k = k, bs = bs, m = m, xt = xt)
      rval$term <- term
      rval$by <- by
      rval$label <- paste("sx(", paste(term, collapse = ",", sep = ""), ")", sep = "")
    }
  }

  return(rval)
}

