
data.frame2txt <- function(x, datafile = "x.txt",
                           r.codefile = "x.R",
                           sps.codefile = "x.sps",
                           df.name = "x",
                           user.missing){
  x.names <- names(x)

  systemfile <- sub("\\....$", "", datafile)

  sink(r.codefile)
  cat(df.name, ' <- read.delim("', datafile, '", quote = "", as.is = TRUE)\n\n',
      sep = "")
  for(column in x.names){
    xx <- x[[column]]
    if(is.factor(xx)){
      xx.levels <- gsub('"', '\\\\"', levels(xx))
      n.levels <- length(xx.levels)
      cat(df.name, "$", column, " <- factor(", df.name, "$", column,
          ", levels = 1:", n.levels, ',\n  labels = c("', sep = "")
      cat(xx.levels[1], '"', sep = "")
      if(n.levels > 1) cat(", ")
      i <- 2
      len <- 2
      n.levels1 <- n.levels - 1
      while(i < n.levels){
	len <- len + nchar(xx.levels[i]) + 4
	if(len > 80){
	  cat("\n    ")
	  len <- nchar(xx.levels[i]) + 6
	}
	cat('"', xx.levels[i], '", ', sep = "")
	i <- i + 1
      }
      if(len > 80) cat("\n  ")
      if(n.levels > 1) cat('"', xx.levels[n.levels], '"', sep = "")
      cat("))\n")
    }
  }
  for(column in x.names){
    xx <- x[[column]]
    xx.label <- attr(xx, "label")
    if(!is.null(xx.label)){
      cat("attr(", df.name, "$", column, ', "label") <- "', xx.label,
          '"\n', sep = "")
    }
  }
  cat("save(", df.name, ", file = \"", systemfile, ".RData\")\n", sep = "")
  sink()

  sink(sps.codefile)
  cat("GET DATA\n")
  cat("  /TYPE=TXT\n")
  cat("  /FILE='", datafile, "'\n", sep = "")
  cat("  /DELCASE=LINE\n")
  cat("  /DELIMITERS=\"\\t\"\n")
  cat("  /ARRANGEMENT=DELIMITED\n")
  cat("  /FIRSTCASE=2\n")
  cat("  /IMPORTCASE=ALL\n")
  cat("  /VARIABLES=\n")
  for(column in x.names){
    cat("  ", column, " ", sep = "")
    xx <- x[[column]]
    if(is.character(xx)) cat("A", max(nchar(xx)), "\n", sep = "")
    else if(is.factor(xx)){
      nlevs <- length(levels(xx))
      if(nlevs < 10) cat("F1.0\n")
      else if(nlevs > 9 && nlevs < 100) cat("F2.0\n")
      else if(nlevs > 99) cat("F3.0\n")
    } else if(is.numeric(xx)){
        if(sum(grepl("(chron|dates|times)", class(xx))) > 0){
            cat("A", max(nchar(as.character(xx))), "\n", sep = "")
        } else {
            cat("F", max(nchar(as.character(xx))), ".0\n", sep = "")
        }
    } else cat("error: undefined type\n")
  }
  cat("  .\n")
  cat("EXECUTE.\n\n")

  for(column in x.names){
    xx <- x[[column]]
    xx.label <- attr(xx, "label")
    if(!is.null(xx.label))
      cat("VARIABLE LABELS ", column, ' "', xx.label, '" .\n', sep = "")
  }
  cat("\n")

  for(column in x.names){
    xx <- x[[column]]
    if(is.factor(xx)){
      cat("VALUE LABELS ", column, "\n", sep = "")
      xx.levels <- levels(xx)
      len <- length(xx.levels)
      for(i in 1:len){
	if(i < len){
	  cat("  ", i, ' "', xx.levels[i], '"\n', sep = "")
	} else {
	  cat("  ", i, ' "', xx.levels[i], '" .\n', sep = "")
	}
      }
      if(!missing(user.missing)){
	user.missing <- paste("^", user.missing, "$", sep = "")
	nmiss <- 0
	umiss <- numeric()
	i <- 1
	for(lmiss in user.missing){
	  idx <- grep(lmiss, xx.levels)
	  if(length(idx) == 1){
	    nmiss  <- nmiss + 1
	    umiss[i] <- idx
	    i <- i + 1
	  } else {
	    if(length(idx) > 1){
	      msg <- paste(gettext("Repeated labels in ", domain = "R-descr"),
                           column, ": ",  lmiss, sep = "")
	      stop(msg)
	    }
	  }
	}
	if(nmiss > 0){
	  cat("MISSING VALUES", column, "(")
	  cat(umiss, sep = ", ")
	  cat(").\n")
	}
      }
      cat("\n")
    }
  }
  cat("SAVE OUTFILE='", systemfile, ".sav'\n  /COMPRESSED.\n", sep = "")
  sink()

  for(column in x.names)
    if(is.factor(x[[column]])) x[[column]] <- as.numeric(x[[column]])

  write.table(x, file = datafile, quote = FALSE, sep = "\t", col.names = TRUE,
    row.names = FALSE, na = "")
}

