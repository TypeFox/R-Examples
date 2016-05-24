.loadHTMLTemplate <- function(name, params=list()) {
  sourceFile <- system.file("extdata",sprintf("htmlTemplates/%s.txt", name) , package="exreport")
  template <- readChar(sourceFile, file.info(sourceFile)$size)  
  template <- .replaceVarsForContent(template,params)
  #for ( var in names(params))
  #  template <- gsub( sprintf("$%s", var), params[[var]], template, fixed=TRUE)
  
  template
}

.loadLatexTemplate <- function(name, params=list()) {
  sourceFile <- system.file("extdata",sprintf("latexTemplates/%s.txt", name) , package="exreport")
  template <- readChar(sourceFile, file.info(sourceFile)$size)  
  template <- .replaceVarsForContent(template,params)
  # For latex code, we need to replace single and double quotes for `' and ``"
  template <- .sanitizeLatexCode(template)
  
  template
}

####################################
## print2Report implementation.
#
#

.print2exreport <- function(element, id, file, path, target, ...) UseMethod(".print2exreport")

.print2exreport.experiment <- function(element, id, file, path, target) {

  # Format experiment data:
  methods     <- paste(unique(element$data[[element$method]]), collapse = ', ')
  problems    <- paste(unique(element$data[[element$problem]]), collapse = ', ')
  
  parameters <- ""
  
  # Print the parameters list if any
  params <- c()
  if (length(element$parameters) != 0) 
    for (p in element$parameters)
      params <- c(params, paste0(p, ' [', paste0(levels(element$data[[p]]), collapse = ","), ']'))
  
  if (length(element$configuration) != 0) 
    params <- c(params, element$configuration)
  
  if(target == "html")
    parameters <- .nestedList2HTML(params, numbered=F)
  else
    parameters <- .nestedList2Latex(params)
  
  outputs     <- paste(element[["outputs"]], collapse = ', ')
  
  if(target == "html")
    history <- .nestedList2HTML(element$h)
  else
    history <- .nestedList2Latex(element$h)
  
  templateParams <- c(id = id,
                      element$tags,
                      methods = methods, 
                      problems = problems, 
                      parameters = parameters, 
                      outputs = outputs, 
                      history = history)
  
  if(target == "html")
    content <- .loadHTMLTemplate("experiment", templateParams)
  else
    content <- .loadLatexTemplate("experiment", templateParams)
  
  cat( content, file = file)
}

.print2exreport.testMultiple <- function(element, id, file, path, target) {
  
  # TODO: Make proper function
  
  .print2exreport(element$friedman, id, file, path, target)
}

.print2exreport.testFriedman <- function(element, id, file, path, target) {
  
  templateParams <- c(id = id,
                      element$tags)
  if(target=="html")
    content <- .loadHTMLTemplate("friedman", templateParams)
  else
    content <- .loadLatexTemplate("friedman", templateParams)
  
  cat( content, file = file)
}

.print2exreport.testPaired <- function(element, id, file, path, target) {
  
  templateParams <- c(id = id,
                      element$tags,
                      worstMethod = element$worstMethod,
                      bestMethod = element$bestMethod)
  if(target=="html")
    content <- .loadHTMLTemplate("wilcoxon", templateParams)
  else
    content <- .loadLatexTemplate("wilcoxon", templateParams)
    
  cat( content, file = file)
}

.loadHTMLTemplate <- function(name, params=list()) {
  sourceFile <- system.file("extdata",sprintf("htmlTemplates/%s.txt", name) , package="exreport")
  template <- readChar(sourceFile, file.info(sourceFile)$size)  
  template <- .replaceVarsForContent(template,params)
  #for ( var in names(params))
  #  template <- gsub( sprintf("$%s", var), params[[var]], template, fixed=TRUE)
  
  template
}

.print2exreport.exTabular <- function(element, id, file, path, target) {
  
  # Configure tabular structures:
  tables <- element$tables
  formats <- element$formats
  
  # Split tables by maximum number of columns
  colHeader <- tables[[1]][, 1, drop=FALSE]
  colFormatHeader <- formats[[1]][, 1, drop=FALSE]
  maxCol <- ncol(tables[[1]])
  colIndex <- 2
  
  if (element$tableSplit < 1)
    tableFolds <- 1
  
  # Just to tell how many tables still need to be created
  tableFolds <- element$tableSplit
  
  htmlTable  <- ""
  latexTable <- ""
  
  while (colIndex <= maxCol)
  {
    # The number of columns of the actual table (not including the first one)
    numColsTable <- round( (maxCol-colIndex+1) / tableFolds )
    endIndex <- colIndex+numColsTable-1
    auxTables <- lapply(tables, FUN = function(tab){
      cbind(colHeader,tab[,colIndex:endIndex,drop=FALSE])
    })
    auxFormats <- lapply(formats, FUN = function(tab){
      cbind(colFormatHeader,tab[,colIndex:endIndex,drop=FALSE])
    })
    if(target=="html")
      htmlTable  <- paste0(htmlTable, "\n<br/><br/>\n", .formatDataFrame(tables = auxTables, formats = auxFormats, src = "html"))
    
    newLatexTable <- .formatDataFrame(tables = auxTables, formats = auxFormats, src = "latex")
    if(target=="pdf")
      # In this case, we append the code for centering and scaling the table
      newLatexTable <- paste0("\\exTable{",newLatexTable,"}")
    latexTable <- paste0(latexTable, "\n\n", newLatexTable)
    
    colIndex   <- colIndex + numColsTable
    tableFolds <- tableFolds - 1
  }
  
  if(target=="html"){
    templateParams <- c(id = id,
                        element$tags,
                        htmlTable = htmlTable,
                        latexTable = latexTable
                        )
    # Generate final output, it is conditioned on tableType
    if (element$tableType=="plain")
      content <- .loadHTMLTemplate("exTabular_plain", templateParams)
    else if (element$tableType=="phtest")
      content <- .loadHTMLTemplate("exTabular_phtest", templateParams)
  }
  else{
    templateParams <- c(id = id,
                        element$tags,
                        latexTable = latexTable
                        )
    # Generate final output, it is conditioned on tableType
    if (element$tableType=="plain")
      content <- .loadLatexTemplate("exTabular_plain", templateParams)
    else if (element$tableType=="phtest")
      content <- .loadLatexTemplate("exTabular_phtest", templateParams)
  }
  cat( content, file = file )
}

.print2exreport.exPlot <- function(element, id, file, path, target) {
  
  figPath <- paste(path, "/img/", sep="")
  
  # Generate the pictures
  # The first three only if target is HTML. For PDF we only use .eps
  if(target=="html"){
    pdf(paste(figPath, id, ".pdf", sep=""), width=11.69, height=8.27)
    print(element$plot)
    a <- dev.off()
    
    png(paste(figPath, id, ".png", sep=""), width=11.69, height=8.27, units = "in", res=72)
    print(element$plot)
    a <- dev.off()
    
    svg(paste(figPath, id, ".svg", sep=""), width=11.69, height=8.27)
    print(element$plot)
    a <- dev.off()
  }
  
  setEPS()
  postscript(paste(figPath, id, ".eps", sep=""), width=11.69, height=8.27)
  print(element$plot)
  a <- dev.off()
  
  templateParams <- c(id = id,
                      element$tags
                     )
  
  if(target=="html")
    content <- .loadHTMLTemplate("exPlot", templateParams)
  else
    content <- .loadLatexTemplate("exPlot", templateParams)
  
  cat( content, file = file)
}