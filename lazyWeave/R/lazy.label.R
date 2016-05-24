#' @describeIn lazy.ref
#' @export lazy.label
#' 

lazy.label <- function(label){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html", "markdown")) stop("option(\"lazyReportFormat\") must be either 'latex', 'html', or 'markdown'")
  
  
  #*** Construct the comment with the function call
  comment.char <- if (reportFormat == "latex") {
    if (getOption("lazyWeave_latexComments") == "latex") c("%%", "") else c("<!-- ", " -->")
  }  else if (reportFormat == "html") c("<!--", "-->")
  
  fncall <- paste(comment.char[1], paste(deparse(match.call()), collapse=" "), comment.char[2], "\n")
  
  
  if (reportFormat == "latex") return(paste("\\label{", label, "}", sep=""))
  
  if (reportFormat == "html") return(paste(fncall, "<a name='", label, "'></a>\n", sep=""))
  
  if (reportFormat == "markdown"){ 
    warning("labelling and referencing are not currently available in markdown")
    return("")
  }
}
