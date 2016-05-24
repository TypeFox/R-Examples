tex.table <- 
  function(dm, bare=FALSE, prec = if (bare) "NA" else 2, rnames = if (bare) "-1" else dimnames(dm)[[1]], cnames = if (bare) "-1" else dimnames(dm)[[2]], caption = NULL, label = NULL, tpos = "b", stretch = NULL, adjust = "r", file = NULL)
{  # tex.table.s
## 
## function converts a data matrix into LATEX code
##
## Arguments:
##
## dm: R data matrix
## prec: precision of rounding within the LATEX table
##            if NA, then no transformation to numeric is done
## rnames: row names
## cnames: column names
## caption: caption for LATEX table, default: no caption
## label: LATEX label for the table,  default: no lable
## tpos: position of captions: above or below table
##       "a" for above, "b" for below
## stretch: optional vector with two entries, giving the baselinestretch for the 
##          caption (stretch[1]) and the colums of the table (stretch[2])
##          default: no adjustment of baselinestretch
## adjust: adjusts the columns of the LATEX table, default: "r" (right), 
##          also possible: "l" (left) and "c" (centre) or user defined:
##          "adjust=c("l","c","r",...)" yields {l|cr...}
## file: output file, default: printout in console
##
## Value:
##
## LATEX-Code 
##
##
  exact1 <- function(x, prec) {
    if (!is.na(prec) & is.numeric(x)) {
      x <- round(x, prec)
      for (i in prec:1) {
        if ((i == prec) & ((x * 10^prec) %% (10^prec) == 0))
          return(cat(" & ", x, ".", rep(0, prec), sep = ""))
        else if(round((x * 10^prec) %% (10^i), 5) == 0)
          return(cat(" & ", x, rep(0, i), sep = ""))
      }
      return(cat(" & ", x, sep = ""))
    }
    return(cat(" & ", x, sep = ""))
  } # exact1

  printout1 <- function(dm, rnames, prec) {
    for (i in 1:nrow(dm)) {
      if (is.null(rnames)) cat("Row", i)
      else if (rnames[1]!="-1") cat(rnames[i])
        for (j in 1:ncol(dm)) exact1(dm[i, j], prec = prec)

      cat("\\\\ \n")
    }
  } # printout1
        
  if (!is.null(file)) sink(file = file)
  ncols <- ncol(dm)
  if (length(adjust) == 1) {
    if(adjust != "r" && adjust != "l" && adjust != "c")
        stop("adjust should be r, l or c")
    ausr <- rep(adjust, length = ncols + 1)
  } else {
    if (length(adjust) == ncols + 1) ausr <- adjust
    else 
      stop("user defined adjust is required to be of length ncol(daten) + 1")
  }

  cat("\\documentclass[a4paper,11pt,leqno]{scrartcl}\n\\usepackage{geometry}\n\\begin{document}\n\\begin{table}[ht]\n\\begin{center}\n")
  if (tpos == "a" && !is.null(caption)) {
    if(length(stretch) == 2)
      cat("\\renewcommand{\\baselinestretch}{", stretch[1], "} \\normalsize\n", sep = "")
    cat("\\caption{", sep = "")
    if (!is.null(label)) cat("\\label{", label, "}", sep = "")
    cat(caption, "}\n", sep = "")
  }
  if(length(stretch) == 2)
    cat("\\renewcommand{\\baselinestretch}{", stretch[2],
            "} \\normalsize\n\\vspace{0.3cm}\n", sep = "")
  cat("\\begin{footnotesize}\n\\begin{tabular}")
  cat("{",ausr[1],"|",ausr[2:length(ausr)],"}\n\\hline\n", sep = "")
  if (!(length(cnames)==1 && cnames=="-1")) {
    if (is.null(cnames)) cat("", paste(1:ncols), sep=" & Col ")
    else  for (i in 1:ncols) cat(" &", cnames[i])
      cat("\\\\ \\hline\n")
  }
  printout1(dm = dm, rnames = rnames, prec = prec)
  cat("\\hline\n\\end{tabular}\n\\end{footnotesize}\n")
  if (tpos == "b" && !is.null(caption)) {
    if(length(stretch) == 2)
      cat("\\renewcommand{\\baselinestretch}{", stretch[1], "} \\normalsize\n", sep = "")
    cat("\\caption{")
    if (!is.null(label)) cat("\\label{", label, "}", sep = "")
    cat(caption, "}\n", sep = "")
  }
  cat("\\end{center}\n\\end{table}\n\n\\end{document}\n")
  if (!is.null(file)) sink()
  invisible()
} # tex.table

tex.tab0 <- 
  function(dm, prec = 2, rnames = NULL,
    cnames = NULL, caption = NULL, label = NULL,
    tpos = "b", stretch = NULL, adjust = "r", file = NULL)
{ # tex.tab0
  exact2 <- function(x, prec) {
    x <- round(x, prec)  
    for (i in prec:1) {
      if ((i == prec) & ((x * 10^prec) %% (10^prec) == 0)) 
        return(cat(" & ", x, ".", rep(0, prec), sep = ""))
    else if (round((x * 10^prec) %% (10^i), 5) == 0) 
      return(cat(" & ", x, rep(0, i), sep = ""))
    }
    return(cat(" & ", x, sep = ""))
  } # exact2  
    
  printout2 <- function(dm, rnames, prec) {
    for (i in 1:nrow(dm)) {
      if (is.null(rnames)) cat("Row", i)
      else cat(rnames[i])
        for (j in 1:ncol(dm)) exact2(dm[i, j], prec = prec)
      cat("\\\\ \n")
    }
  } #printout2

 if (!is.null(file)) sink(file = file)
 ncols <- ncol(dm)
    if (length(adjust) == 1) {
        if(adjust != "r" && adjust != "l" && adjust != "c")
            stop("adjust should be r, l or c")
        ausr <- rep(adjust, length = ncols + 1)
    }
    else {
        if (length(adjust) == ncols + 1) 
            ausr <- adjust
        else 
            stop("user defined ajust is required to be of length ncol(daten) + 1")
    }
   
  cat("\\documentclass[a4paper,11pt,leqno]{scrartcl}\n\\usepackage{geometry}\n\\begin{document}\n\\begin{table}[ht]\n\\begin{center}\n")
    if (tpos == "a" && !is.null(caption)) {
        if(length(stretch) == 2)
            cat("\\renewcommand{\\baselinestretch}{", stretch[1], 
                "} \\normalsize\n", sep = "")
        cat("\\caption{", sep = "")
        if(!is.null(label)) cat("\\label{", label, "}", sep = "")
        cat(caption, "}\n", sep = "")
    }
    if(length(stretch) == 2)
        cat("\\renewcommand{\\baselinestretch}{", stretch[2],
            "} \\normalsize\n\\vspace{0.3cm}\n", sep = "")
    cat("\\begin{tabular}")
    cat("{",ausr[1],"|",ausr[2:length(ausr)],"}\n", sep = "")
    if (is.null(cnames))
        cat("", paste(1:ncols), sep=" & Col ")
    else
        for (i in 1:ncols) cat(" &", cnames[i])
    cat("\\\\ \\hline\n")
    printout2(dm = dm, rnames = rnames, prec = prec)
    cat("\\end{tabular}\n")
    if (tpos == "b" && !is.null(caption)) {
        if(length(stretch) == 2)
            cat("\\renewcommand{\\baselinestretch}{", stretch[1], 
                "} \\normalsize\n", sep = "")
        cat("\\caption{")
        if(!is.null(label)) cat("\\label{", label, "}", sep = "")
        cat(caption, "}\n", sep = "")
    }
    cat("\\end{center}\n\\end{table}\n\\end{document}\n")
    if (!is.null(file)) sink()
    invisible()
}  # tex.tab0
