#' A simpler latex output of the latex.anova.rms 
#' 
#' The original problem is that the anova default function
#' output is very detailed and cause a complaint in Sweave/knitr that
#' \\hbox is overfull. It basically changes capitalized TOTAL,
#' TOTAL INTERACTION and TOTAL NONLINEAR INTERACTION into lower
#' case letters. It also deletes the (Factor + Higher Order Factors).
#'   
#' @param anova_output An object from the \code{\link{anova}()} function 
#' @param subregexps A 2 column matrix with sub() regular expressions
#'   to search for and their substitutions. The regular expression
#'   should be located in column 1 and the substitution in column
#'   2.
#' @param digits Number of digits in using the round
#' @param rowlabel The label of the rows
#' @param pval_lim.sig The threshold before setting "<", default is < 0.0001 
#' @param ... Passed on to latex() or htmlTable
#' @return void See the latex() function 
#' 
#' @example inst/examples/simpleRmsAnova_example.R
#' 
#' @import htmlTable
#' @rdname SimpleRmsAnova
#' @export
simpleRmsAnova <- function(anova_output, subregexps, digits=4, 
                           pval_lim.sig = 10^-4, rowlabel="", ...){
  
  if (!inherits(anova_output, "anova.rms"))
    if (inherits(anova_output, "rms")){
      anova_output <- anova(anova_output) 
    }else{
      stop("You must provide either an rms-regression object or an anova rms output for this to work")
    }
  
  rnames <- names(attr(anova_output, "which"))
  if (!missing(subregexps)){
    if (is.matrix(subregexps) && NCOL(subregexps) == 2){
      for (i in 1:NROW(subregexps))
        rnames <- sub(subregexps[i, 1], subregexps[i, 2], rnames)
    }else if (length(subregexps) == 2){
      rnames <- sub(subregexps[1], subregexps[2], rnames)
    }else{
      stop("Regression substitution through the subregexps should be in",
           " matrix format with two columns or in vector format with length of 2")
    }
  }
  rnames <- sub("TOTAL", "Total", rnames)
  rnames <- sub("INTERACTION", "interaction", rnames)
  rnames <- sub("NONLINEAR", "nonlinear", rnames)
  rnames <- sub("\\(Factor\\+Higher Order Factors\\)", "", rnames)
  
  mtrx <- as.matrix(anova_output)
  pvals <- mtrx[,ncol(mtrx)]
  mtrx <- mtrx[,-ncol(mtrx)]
  # The digits differ per column and we need
  # to handle NA:s
  mtrx <- apply(mtrx, MARGIN=2, 
    FUN=function(x, digits) {
      ret <- c()
      for(val in x){
        if (is.numeric(val) == FALSE){
          ret <- append(ret, val)
        }else if (is.na(val)){
          ret <- append(ret, "")
        }else if (round(val, digits) == 0){
          ret <- append(ret, sprintf(sprintf("%%.%df", digits), 0))
        }else{
          ret <- append(ret, format(val, digits=digits))
        }
      }
      return(ret)
    },
    digits=digits)
  pvals <- txtPval(pvals, lim.sig = pval_lim.sig)
  mtrx <- cbind(mtrx, pvals)
  rownames(mtrx) <- rnames
  number_of_total_rows <- length(grep("TOTAL", names(attr(anova_output, "which"))))
  attr(mtrx, "n.rgroup") <- c(NROW(mtrx)-number_of_total_rows, number_of_total_rows)
  attr(mtrx, "rgroup") <- c("Variables", "Total")
  attr(mtrx, "rowlabel") <- rowlabel
  attr(mtrx, "other") <- list(...)
  class(mtrx) <- c("simpleRmsAnova", class(mtrx))
  return(mtrx)
}

setClass("simpleRmsAnova", contains = "matrix")

#' @param x The output object from the SimpleRmsAnova function 
#' @rdname SimpleRmsAnova
#' @method print simpleRmsAnova
#' @param html If HTML output through the htmlTable should be used 
#'   instead of traditional latex() function
#' @export
#' @import htmlTable
#' @importFrom Gmisc fastDoCall
#' @keywords internal
print.simpleRmsAnova <- function(x, html=TRUE, ...){
  dots <- list(...)
  html = TRUE
  if ("html" %in% names(dots)){
    html <- dots[["html"]]
    dots[["html"]] <- NULL
    
  }else if (length(attr(x, "html")) > 0){
    html <- attr(x, "html")
    attr(x, "html") <- NULL
  }
  
  call_list <- list(x = x, 
    n.rgroup=attr(x, "n.rgroup"), 
    rgroup=attr(x, "rgroup"),
    rowlabel=attr(x, "rowlabel"))
  
  if (html){
    call_list[["rnames"]] <- sub("^ ", "&nbsp;&nbsp;", rownames(x))
    call_list[["x"]] <- gsub("<", "&lt;", call_list[["x"]])
  }else{
    call_list[["rowname"]] <- sub("^ ", "\\\\hspace{3 mm}", rownames(x))
  }
  
  if (length(attr(x, "other")) > 0){
    other <- attr(x, "other")
    for (option in names(other))
      if (nchar(option) > 0) call_list[option] <- other[[option]]
  }
  
  if (length(dots) > 0){
    for (option in names(dots))
      if (nchar(option) > 0) call_list[option] <- dots[[option]]
  }
  
  fastDoCall(ifelse(html, "htmlTable", "latex"), call_list) %>%
    print
}
