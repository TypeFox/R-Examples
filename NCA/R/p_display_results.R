p_display_results <-
function (results.nca, loop.data, prefix="out", pdf="FALSE") {
  data  <- c(sapply(results.nca, p_pretty, "ceiling", FALSE, "NA", "auto"),
             rep(p_pretty_number(loop.data$scope, "NA", "auto"), length(results.nca)),
             sapply(results.nca, p_ceiling_line),
             sapply(results.nca, p_pretty_accuracy, loop.data),
             sapply(results.nca, p_pretty, "effect", FALSE, "", 2),
             sapply(results.nca, p_pretty, "ineffs$x", TRUE),
             sapply(results.nca, p_pretty, "ineffs$y", TRUE))
  table <- matrix(data, ncol=7)
  
  # Add col/row names, also remove the OLS row
  rownames(table) <- names(results.nca)
  m <- match("ols", names(results.nca), 0)
  if (m > 0) {
    table <- matrix(table[-m,], ncol=7)
    rownames(table) <- names(results.nca)[-m]
  }
  colnames(table) <- p_TABLERESULT
  
  # Nothing to print
  if (nrow(table) == 0) {
    return()
  }
  
  if (pdf) {
    p_new_pdf(prefix, "result", loop.data$id.x, loop.data$id.y)
    for (i in 1:length(table[,7])) {
      colNames <- colnames(table)
      rowName <- rownames(table)[i]
      tmpTable <- matrix(table[i,], ncol=1)
      colnames(tmpTable) <- c(rowName)
      rownames(tmpTable) <- colNames
      
      textplot (tmpTable, cex=1, halign="left", valign="top",
                mar=c(0,0,4,0), cmar=0.5, rmar=1)
      title(paste("X :", loop.data$names[1], "\nY :", loop.data$names[2]), cex.main=1) 
    }
    dev.off()
  } else {
    width <- p_calc_width(table)
    p_new_window(title="Result table", width=width)
    par(mfrow=c(1, 1))
    par(family="")    
    textplot(t(table), cex=1, halign="left", valign="top", 
             mar=c(0,0,4,0), cmar=0.5, rmar=1)
    title(paste("X :", loop.data$names[1], "\nY :", loop.data$names[2]), cex.main=1)
  }
}

p_pretty_accuracy <- 
function (nca, loop.data) {  
  nObservations <- min(length(loop.data$x), length(loop.data$y))
  tmp <- 100 * (nObservations - nca[["above"]]) / nObservations
  s = p_pretty_number(tmp, "", 1)
  if (substr(s, nchar(s), nchar(s)) == "0") {
    s <- substr(s, 1, nchar(s)-2)
  }
  return(sprintf("%s%%", s))
}

p_pretty <- 
function (nca, field, percentage=FALSE, default="NA", digits=0) {
  if (length(grep("$", field, fixed=TRUE)) > 0) {
    parts <- unlist(strsplit(field, "[$]"))    
    value <- nca[[parts[1]]][[parts[2]]]
  } else {
    value <- nca[[field]]    
  }
  
  if (is.na(value)) {
    return(default)
  }
  if (percentage) {
    return(sprintf("%s%%", p_pretty_number(value, default, digits)))
  } else {
    return(sprintf("%s", p_pretty_number(value, default, digits)))
  }  
}

p_ceiling_line <- 
function (nca) {
  slope <- p_pretty_number(nca[["slope"]], "NA", "auto")
  intercept <- p_pretty_number(abs(nca[["intercept"]]), "NA", "auto")
      
  if (slope == "NA" || is.na(intercept)) {
    return("NA")
  }
      
  if (nca[["intercept"]] < 0) {
    return(sprintf("Yc = %s Xc - %s", slope, intercept))
  } else {
    return(sprintf("Yc = %s Xc + %s", slope, intercept))
  }
}

p_calc_width <-
function (table) {
  width <- 1.5
  
  if (is.null(dim(table))) {
    table <- matrix(table, nrow=1)
  }
  
  for (i in 1:length(table[,7])) {
    tmp <- table[i,7]
    if (tmp == "NA") {
      tmp <- rownames(table)[i]
    }
    width <- width + 0.24 + 1.2 * strwidth(tmp, font = 12, units = 'in')
  }
        
  return( width )
}
