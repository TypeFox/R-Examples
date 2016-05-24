crossTable <- function(..., deparse.level = 2){
  table <- table(..., deparse.level = deparse.level)
  class(table) <- c("CrossTable", "table")
  return(table)
}

summary.CrossTable <- function(object, digits=3, latex=FALSE, ...){
  vcd    <- suppressWarnings(require(vcd, quietly=TRUE))
  x      <- object
  class(x) <- "table"
  sep    <- ifelse(latex, "&", " ")
  twoDimTable <- function(x, digits=3, width=6){
    output <- NULL
    dim    <- dim(x)
    dimnames <- dimnames(x)
    varnames <- names(dimnames)
    if(latex) varnames[2] <- sprintf("\\multicolumn{%s}{c}{%s}", dim[2], varnames[2])
    x <- addmargins(x, margin=1)
    p <- prop.table(x, margin=1) * 100
    x <- addmargins(x, margin=2)
    p <- addmargins(p, margin=2)
    p[is.nan(p)] <- 0
    
    
    rowcat <- paste(c(dimnames[[1]], " "), " ", sep="\t", collapse="\t")
    rowcat <- strsplit(rowcat, "\t")[[1]]
    rowcat[length(rowcat)-1] <- "Total"
    rowcat <- format(c(" ", " ", varnames[1], rowcat), justify="left")
    #    rowvar <- c(" ", " ", " ", varnames[1], rep(" ", length(rowcat)-4))
    #    rowvar[length(rowvar)-1] <- "Total"
    #    rowvar <- format(rowvar)
    
    for(i in seq_len(dim[2])){
      count   <- x[, i]
      percent <- format(p[, i], digits=digits)
      if(latex)
        percent <- paste(percent, "\\%", sep="")
      else
        percent <- paste(percent, "%", sep="")
      
      col <- paste(count, percent, sep="\t", collapse="\t")
      col <- strsplit(col, "\t")[[1]]
      col <- format(col, justify="right", width=width)
      col <- c(dimnames[[2]][i], col)
      if(latex) col[1] <- sprintf("\\multicolumn{1}{c}{%s}", col[1])
      col <- format(col, justify="centre")
      if(is.null(output))
        output <- col
      else
        output <- paste(output, col, sep=sep)
    }
    i <- dim[2]+1
    count   <- x[, i]
    percent <- format(p[, i], digits=digits)
    if(latex)
      percent <- paste(percent, "\\%", sep="")
    else
      percent <- paste(percent, "%", sep="")
    
    col <- paste(count, percent, sep="\t", collapse="\t")
    col <- strsplit(col, "\t")[[1]]
    col <- format(col, justify="right", width=width)
    if(latex){
      col <- c(" ", " ", "\\multicolumn{1}{c}{Total}", col)      
    } else {
      col <- c(" ", " ", "Total", col)      
    }
    col <- format(col, justify="centre")
    
    nchar  <- nchar(output[1], type="width")
    line1  <- paste(rep("-", nchar), collapse="")
    output <- format(c(varnames[2], line1, output), justify="centre")
    output <- paste(output, col, sep=sep)
    output <- paste(rowcat, output, sep=sep)
    nchar  <- nchar(output[1], type="width")
    
    #    output <- paste(rowvar, output, sep=sep)
    nchar  <- nchar(output[1], type="width")
    if(latex) {
      output <- paste(output, "\\\\")
      line1 <- "\\midrule"
      line2 <- "\\toprule"
      line3 <- "\\bottomrule"      
    } else {
      line1  <- paste(rep("-", nchar), collapse="")
      line2  <- paste(rep("=", nchar), collapse="")
      line3  <- paste(rep("=", nchar), collapse="")      
    }
    
    output <- c(line2, output[1:3], line1, output[4:length(output)], line3)
    output <- c(output[1:(length(output)-3)], line1,
                output[(length(output)-2):length(output)])
    return(output)
  }
  
  dim <- dim(x)
  varnames <- names(dimnames(x))
  if(length(dim) == 2) {
    # Two dimensional
    output <- twoDimTable(x)
    if(latex) output[3] <- sprintf("\\cline{%s-%s}", 2, 2+dim[2]-1)
    output <- paste(output, collapse="\n")
    if(latex){
      
      output <- sprintf("\\begin{table}[htbp]
\\centering
\\caption{%s $\\times$ %s}
\\begin{tabular}{l%s}
%s
\\end{tabular}
\\end{table}",
                      varnames[1], varnames[2],
                        paste(rep("r",dim[2]+1), collapse=""), output)
    }
    
    
    cat(output, fill=TRUE)
    cat("\n")
    cat("Chi-Square Test for Independence", fill=TRUE)
    if(vcd){
      print(summary(assocstats(x)))
    } else {
      cat("\n")
      print(summary.table(x))
    }
    
    
  } else {
    # Three Dimensional
    stratumcat <- dimnames(x)[[1]]
    stratumvar <- varnames[1]
    stratumcat <- format(c(stratumvar, stratumcat, "Total"), justify="left")
    stratumvar <- stratumcat[ 1]
    stratumcat <- stratumcat[-1]
    output <- list()
    col    <- list()
    width  <- nchar(as.character(max(x)))
    width  <- ifelse(width > 6, width, 6)
    for(i in seq_len(dim[1])) {
      x.tmp <- as.table(x[i, , ])
      output[[i]] <- twoDimTable(x.tmp, width=width)
    }
    total <- margin.table(x, c(2, 3))
    output[[dim[1]+1]] <- twoDimTable(total, width=width)
    
    output.header <- output[[1]][2:4]
    if(latex) output.header[2] <- sprintf("\\cline{%s-%s}", 3, 3+dim[3]-1)
    output.header[1] <- paste(
      paste(rep(" ", nchar(stratumvar)+2), collapse=""),
      output.header[1], sep=sep)
    output.header[2] <- paste(
      paste(rep(" ", nchar(stratumvar)+2), collapse=""),
      output.header[2], sep=sep)
    output.header[3] <- paste(stratumvar, output.header[3], sep=sep)
    
    output <- lapply(output, function(x) return(x[ -c(1:5, length(x))]))
    for(i in seq_along(output)) {
      col         <- c(stratumcat[i], rep(" ", length(output[[i]])-1))
      col         <- format(col, justify="left")
      output[[i]] <- paste(col, output[[i]], sep=sep)
      nchar  <- nchar(output[[i]][1], type="width")
      if(latex)
        line <- "\\midrule"
      else
        line <- paste(rep("-", nchar), collapse="")
      
      output[[i]] <- c(output[[i]], line)
    }
    output <- unlist(output)
    output <- output[-length(output)]
    #    col    <- c(stratumvar, rep(" ", length(output)-1))
    #    col    <- format(col)
    #    output <- paste(col, output, sep=sep)
    
    nchar  <- nchar(output[1], type="width")
    if(latex) {
      line1 <- "\\midrule"
      line2 <- "\\toprule"
      line3 <- "\\bottomrule"      
    } else {
      line1  <- paste(rep("-", nchar), collapse="")
      line2  <- paste(rep("=", nchar), collapse="")
      line3  <- paste(rep("=", nchar), collapse="")      
    }
    output <- c(line2, output.header, line1, output, line3)
    if(latex){
      output <- gsub("&\\\\cline",   "\\\\cline", output)
      output <- gsub("&\\\\midrule", "\\\\midrule", output)      
    }
    output <- paste(output, collapse="\n")
    if(latex){
      
      output <- sprintf("\\begin{table}[htbp]
\\centering
\\caption{%s $\\times$ %s $\\times$ %s}
\\begin{tabular}{ll%s}
%s
\\end{tabular}
\\end{table}",
                      varnames[1], varnames[2], varnames[3],
                        paste(rep("r",dim[3]+1), collapse=""),output)
    }
    
    cat(output, fill=TRUE)
    cat("\n")
    cat("Chi-Square Test for Independence", fill=TRUE)
    cat("\n")
    for(i in seq_len(dim[1])) {
      x.tmp <- as.table(x[i, , ])
      cat(sprintf("%s : %s", names(dimnames(x))[1], stratumcat[i]), fill=TRUE)
      
      if(vcd){
        print(summary(assocstats(x.tmp)))
      } else {
        cat("\n")
        print(summary.table(x.tmp))
      }
      cat("\n")
    }
    cat("Total", fill=TRUE)

    if(vcd){
      print(summary(assocstats(margin.table(x, c(2, 3)))))
    } else {
      cat("\n")
      print(summary.table(margin.table(x, c(2, 3))))
    }
    cat("\n")
  }
  if (!vcd) {
    message("Please install vcd package to output Cramer's V.")
  }
}
