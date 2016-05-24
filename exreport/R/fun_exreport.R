#' Create a new exreport document
#'
#' This function inits a new exreport document to start adding elements for 
#' later rendering.
#'
#' @export
#' @param title A string representing a short title for this document
#' @return an empty exreport document
#' 
#' @seealso exreportRender, exreportAdd
#' 
#' 
exreport <- function(title) {
  r <- .exReport(title)
  r
}

#' Add elements to an existing exreport document
#'
#' This function allows to add one or more reportable objects to an exisiting 
#' exreport document.
#'
#' @export
#' @param rep an exreport object in which the elem will be added
#' @param elem a reportable object or a list of them
#' @return an extended exreport document
#' 
#' 
#' @examples
#' 
#' # Create an empty document:
#' report <- exreport("Test document")
#' 
#' # Create a reportable object (an experiment)
#' experiment <- expCreate(wekaExperiment, name="test-exp", parameter="fold")
#' 
#' # Add this object to the document
#' exreportAdd(report, experiment)
#' 
exreportAdd <- function(rep, elem) {
  
  if ( !is.exReport(rep) )
    stop(.callErrorMessage("wrongParameterError", "rep", "exReport"))
  
  if (!is.reportable(elem) & !(is.list(elem) & all(sapply(elem, FUN=function(x) "reportable" %in% class(x))) ) )
    stop(.callErrorMessage("reportableError"))
  
  if (is.reportable(elem))
    rep$content[[length(rep$content)+1]] <- elem
  else
    rep$content <- append(rep$content, elem)
    
  
  rep
}


#' Render an exreport document
#'
#' This function renders an existing exreport object to a given file and format.
#'
#' @export
#' @param rep The exreport object to be rendered
#' @param destination Path to the rendered file. If NULL, it uses a temporary directory
#' @param target The format of the target rendering. HTML and PDF are allowed.
#' @param safeMode Denies or allows (TRUE or FALSE) output files overwriting
#' @param visualize Visualize the generated output or not
#' @return an experiment object
#' 
exreportRender <- function(rep, destination=NULL, target="html", safeMode=TRUE, visualize=TRUE) {
  
  if ( !is.exReport(rep) )
    stop(.callErrorMessage("wrongParameterError", "rep", "exReport"))
  
  target <- tolower(target)
  if ( !(target %in% c("html", "pdf")) )
    stop(.callErrorMessage("wrongParameterError", "target", "[\"html\" or \"pdf\"]"))
  
  # If destination is NULL, it uses a temporary directy
  if(is.null(destination))
    path <- tempdir()
  else{
    path <- destination
    # If destination ends in a slash, we remove it (as we add it when it is necessary).
    if ( substr(path,start=nchar(path),stop=nchar(path))=='/' )
      path <- substr(path,start=1,stop=nchar(path)-1)
  }
  
  # If the folder does not exist, stop
  if ( !file.exists(path) )
    stop(.callErrorMessage("invalidFolderError", "rep", "exReport"))
  
  # If the folder 'exreport_output' exists inside path folder,
  # and safeMode equals TRUE, and destination!=NULL, stop
  if( file.exists(sprintf("%s/exreport_output", path)) && safeMode && !is.null(destination))
    stop(.callErrorMessage("directoryAlreadyExists", sprintf("%s/exreport_output", path)))
  # Else creates the directory
  else
    dir.create( file.path(path, "exreport_output"), showWarnings = FALSE )
  
  # Now we will work inside the directory path/exreport_output
  path <- sprintf("%s/exreport_output",path)
  
  # Now generate the propper file depending on target
  if (target == "html")
    genFile <- .exreportRenderHTML(rep, path)
  else
    genFile <- .exreportRenderPDF(rep, path)
  
  # Show the content in the default browser if visualize=TRUE
  if(visualize)
    browseURL(genFile)
}

# Auxiliar function that complements the exreportRender one, but specifically
# for HTML target format
.exreportRenderHTML <- function(rep, path){
  
  # open an output file connection
  f <- file( sprintf("%s/report.html", path), "w" )
  
  # make dir structure and copy required files:
  dir.create( file.path(path, "img"), showWarnings = FALSE )
  dir.create( file.path(path, "assets"), showWarnings = FALSE )
  sourceAssets <- system.file("extdata","htmlTemplates/assets/" , package="exreport")
  file.copy( sprintf("%s", sourceAssets), sprintf("%s/", path), recursive=TRUE )
  
  
  # Add the head
  template <- .loadHTMLTemplate( "main_start" , list(title = rep$title, date = toString(Sys.Date())))
  cat( template, file = f )
  
  counter <- 1
  for ( element in rep$content ) {
    id <- sprintf( "elem_%d", counter )
    .print2exreport( element, id, file = f, path = path, target = "html" ) 
    counter <- counter + 1
  }
  
  # Add the ending
  
  # first, generate the index
  indexString <- ""
  counter <- 1
  pattern <- '<li><a href="#elem_%d_hook"><strong>%s:</strong>%s</a></li>'
  for ( element in rep$content ) {
    indexString <- paste(indexString, sprintf(pattern, 
                                        counter,
                                        element$tags$alias, 
                                        ifelse(!is.null(element$tags$target),
                                               element$tags$target,
                                               element$tags$context)), 
                   sep="\n")
    counter <- counter + 1
  }
  
  template <- .loadHTMLTemplate( "main_end", list(r_version = R.version.string, index = indexString) )
  cat( template, file = f )
  
  close(f)
  
  genFile <- sprintf("%s/report.html", path)
  genFile
}


# Auxiliar function that complements the exreportRender one, but specifically
# for PDF target format
.exreportRenderPDF <- function(rep, path){
  
  # open an output file connection
  f <- file( sprintf("%s/report.tex", path), "w" )
  
  # make dir structure and copy required files:
  dir.create( file.path(path, "img"), showWarnings = FALSE )
  dir.create( file.path(path, "assets"), showWarnings = FALSE )
  sourceAssets <- system.file("extdata","latexTemplates/assets/" , package="exreport")
  file.copy( sprintf("%s", sourceAssets), sprintf("%s/", path), recursive=TRUE )
  
  
  # Add the head
  template <- .loadLatexTemplate( "main_start" , list(title = rep$title, date = toString(Sys.Date())))
  cat( template, file = f )
  
  counter <- 1
  for ( element in rep$content ) {
    id <- sprintf( "elem_%d", counter )
    .print2exreport( element, id, file = f, path = path, target = "pdf" ) 
    counter <- counter + 1
  }
  
  # Add the ending
  template <- .loadLatexTemplate( "main_end", list(r_version = R.version.string) )
  cat( template, file = f )
  
  close(f)
  
  # Now we need to compile the latex code to generate the pdf
  prevwd <- getwd()
  setwd(path)
  tools::texi2pdf("report.tex")
  setwd(prevwd)
  
  genFile <- sprintf("%s/report.pdf", path)
  genFile
}
