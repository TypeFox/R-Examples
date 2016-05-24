write.term.info <- function(file, terms, data, object = NULL, contrasts.arg = NULL, 
  xlev = NULL, intcpt = TRUE, rdafile)
{
  warn <- getOption("warn")
  options(warn = -1L)
  nt <- length(terms)

  reclass <- function(x) {
    if(x == "ps.smooth.spec")
      x <- "sm.bayesx"
    if(x == "rw2.smooth.spec")
      x <- "sm.bayesx"
    if(x == "rw1.smooth.spec")
      x <- "sm.bayesx"
    if(x == "bl.smooth.spec")
      x <- "sm.bayesx"
    if(x == "cs.smooth.spec")
      x <- "sm.bayesx"
    if(x == "gk.smooth.spec")
      x <- "geo.bayesx"
    if(x == "gs.smooth.spec")
      x <- "geo.bayesx"
    if(x == "mrf.smooth.spec")
      x <- "mrf.bayesx"
    if(x == "re.smooth.spec")
      x <- "random.bayesx"
    if(x == "random.smooth.spec")
      x <- "random.bayesx"
    x
  }

  if(nt > 0L) {
    for(k in 1L:nt) {
      if(is.sm(terms[k])) {
        a <- parse(text = terms[[k]])[[1L]]
        a_has_xt <- "xt" %in% names(a)
        map <- paste("NULL")
        if(a_has_xt) {
          if("map" %in% names(a$xt))
            map <- paste("\'", a$xt$map, "\'", sep = "")
          if("polys" %in% names(a$xt))
            map <- paste("\'", a$xt$polys, "\'", sep = "")
        }
        te <- eval(parse(text = terms[k]))
        if(!is.null(te$map.name))
          map <- paste("\'", te$map.name, "\'", sep = "")
        if(!is.null(te$xt$map.name))
          map <- paste("\'", te$xt$map.name, "\'", sep = "")
        fby <- FALSE
        if(te$by != "NA") {
          if(!is.character(data)) {
            by <- data[[te$by]]
            if(is.factor(by)) {
              fby <- TRUE
              fnv <- paste("c(", paste("\'", te$by, levels(by), "\'", sep = "", collapse = ","), 
                ")", sep = "")
            }
          }
        }
        israndom <- FALSE
        if(class(te) == "ra.smooth.spec" || class(te) == "re.smooth.spec")
          israndom <- TRUE
        if(fby) {
          te$label <- gsub(")", paste(",by=", te$by, ")", sep = ""), te$label)
          info <- paste("list(term=\'", te$label, "\',pos=", k, ",by=\'", te$by,
            "\',isFactor=FALSE", ",isFactorBy=", fby, ",isFactorByNames=", fnv, 
            ",map=", map, ",israndom=", israndom, ",class=\'", reclass(class(te)),
            "\',bs=\'", gsub(".smooth.spec", "", class(te)),
            "\',call=\'", terms[k], "\')", sep = "")
        } else {
          info <- paste("list(term=\'", te$label, "\',pos=", k, ",by=\'", te$by,
            "\',isFactor=FALSE", ",isFactorBy=", fby, ",map=", map, 
            ",israndom=", israndom, ",class=\'", reclass(class(te)),
            "\',bs=\'", gsub(".smooth.spec", "", class(te)),
            "\',call=\'", terms[k], "\')", sep = "")
        }
      } else {
        sp <- FALSE
        if(grepl(":", terms[k]))
          sp <- TRUE
        if(!is.character(data) && !sp)
          x <- data[[terms[k]]]
        if(is.factor(x) && !sp) {
          m <- model.matrix(as.formula(paste("~ -1 +", terms[k])), data, contrasts.arg, xlev)
          realname <- colnames(m)
          fn <- rmf(realname)
          realname <- paste("c(\'", paste(realname, sep = "", collapse = "\', \'"), "\')", sep = "")
          fnv <- "c("
          nf <- length(fn)
          if(nf > 1L)
            for(i in 1L:(nf - 1L))
              fnv <- paste(fnv, "\'", fn[i], "\',", sep = "")
          fnv <- paste(fnv, "\'", fn[nf], "\')", sep = "")
          xl <- paste(levels(x), collapse = "\',\'")
          xl <- paste("c(\'", xl , "\')", sep = "")
          info <- paste("list(term=\'", terms[k], "\',pos=" , k, 
            ",isFactor=TRUE", ",names=", fnv, ",levels=", xl, ",realname=", realname,
            ",class=\'linear.bayesx\'", ")", sep = "")
        } else {
          info <- paste("list(term=\'", rmf(terms[k]), "\',pos=", k, ",isFactor=FALSE, realname=", 
            paste("\'", terms[k], "\'", sep = ""), ",class=\'linear.bayesx\'", ")", sep = "")
        }
      }
      info <- paste(info,"\n")
      if(k < 2L)
        cat(info, file = file, append = FALSE)
      else
        cat(info, file = file, append = TRUE)
    }
    if(!is.null(object)) {
      if(!is.null(object$YLevels)) {
        YLevels <- paste(object$YLevels, collapse = "\\',\\'")
        YLevels <- paste("c(\\'", YLevels, "\\')", sep = "")
      }
      if(!is.null(object$nYLevels)) {
        nYLevels <- paste(object$nYLevels, collapse = "\\',\\'")
        nYLevels <- paste("c(\\'", nYLevels, "\\')", sep = "")
      }
      if(!is.null(object$order)) {
        ooo <- paste(object$order, collapse = ",")
        ooo <- paste("\'c(", ooo, ")\'", sep = "")
      }
      f <- object$oformula
      save(f, file = rdafile)
      f <- paste(f[[2L]], " ~ ", paste(attr(terms(f),
        "term.labels"), collapse= " + "), sep = "") 
      info <- paste("list(formula=\'", f, "\',", sep = "")
      if(!is.null(object$hlevel))
        object$method <- "HMCMC"
      info <- paste(info, "method=\'", object$method, "\',", sep = "")
      info <- paste(info, "family=\'", object$family, "\',", sep = "")
      info <- paste(info, "iterations=\'", object$iterations, "\',", sep = "")
      info <- paste(info, "step=\'", object$step, "\',", sep = "")
      if(!is.null(object$YLevels))
        info <- paste(info, "YLevels=\'", YLevels, "\',", sep = "")
      if(!is.null(object$nYLevels))
        info <- paste(info, "nYLevels=\'", nYLevels, "\',", sep = "")
      if(!is.null(object$order))
        info <- paste(info, "order=", ooo, ",", sep = "")
      info <- paste(info, "model.name=\'", object$model.name, "\')\n", sep = "")
      cat(info, file = file, append = TRUE)
    }
  }
  options(warn = warn)

  return(invisible(NULL))
}

