independence_table <- function(x, frequency = c("absolute", "relative")) {
  if (!is.array(x))
    stop("Need array of absolute frequencies!")
  frequency <- match.arg(frequency)

  n <- sum(x)
  x <- x / n
  d <- dim(x)

  ## build margins
  margins <- lapply(1:length(d), function(i) apply(x, i, sum))

  ## multiply all combinations & reshape
  tab <- array(apply(expand.grid(margins), 1, prod), d, dimnames = dimnames(x))
  
  if (frequency == "relative") tab else tab * n
}

mar_table <- function(x) {
  if(!is.matrix(x))
    stop("Function only defined for 2-way tables.")
  tab <- rbind(cbind(x, TOTAL = rowSums(x)), TOTAL = c(colSums(x), sum(x)))
  names(dimnames(tab)) <- names(dimnames(x))
  tab
}

table2d_summary <- function(object,
                            margins = TRUE,
                            percentages = FALSE,
                            conditionals = c("none", "row", "column"),
                            chisq.test = TRUE,
                            ...
                            )
{
  ret <- list()
  if (chisq.test)
      ret$chisq <- summary.table(object, ...)
  
  if(is.matrix(object)) {
    
    conditionals <- match.arg(conditionals)
  
    tab <- array(0, c(dim(object) + margins, 1 + percentages + (conditionals != "none")))

    ## frequencies
    tab[,,1] <- if(margins) mar_table(object) else object

    ## percentages
    if(percentages) {
      tmp <- prop.table(object)
      tab[,,2] <- 100 * if(margins) mar_table(tmp) else tmp
    }

    ## conditional distributions
    if(conditionals != "none") {
      tmp <- prop.table(object, margin = 1 + (conditionals == "column"))
      tab[,,2 + percentages] <- 100 * if(margins) mar_table(tmp) else tmp
    }

    ## dimnames
    dimnames(tab) <- c(dimnames(if(margins) mar_table(object) else object),
                       list(c("freq",
                              if(percentages) "%",
                              switch(conditionals, row = "row%", column = "col%")
                              )
                            )
                       )

    ## patch row% / col% margins
    if(conditionals == "row") 
      tab[dim(tab)[1],,2 + percentages] <- NA
    
    if(conditionals == "column")
      tab[,dim(tab)[2],2 + percentages] <- NA
    
    ret$table <- tab
  }    

  class(ret) <- "table2d_summary"
  ret
}

print.table2d_summary <- 
function (x, digits = max(1, getOption("digits") - 3), ...) 
{
  if (!is.null(x$table))
    if(dim(x$table)[3] == 1)
      print(x$table[,,1], digits = digits, ...)
    else
      print(ftable(aperm(x$table, c(1,3,2))), 2, digits = digits, ...)
  
  cat("\n")
  
  if (!is.null(x$chisq))
    print.summary.table(x$chisq, digits, ...)
  invisible(x)
}

