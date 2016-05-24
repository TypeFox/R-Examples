#' @name lazy.build
#' @export lazy.build
#' 
#' @title Compile a PDF or HTML Report
#' @description Executes a command to build a pdf file based on a .tex or 
#'   .html file.  For HTML files, compiles the figure into a subfolder 
#'   and places all of the contents into a zip file
#'   
#' @param filename Character string giving the location of the file to be built.
#'   Usually a .tex or .html file.
#' @param pdf.zip filename where the completed document should be saved. Usually a 
#'   .pdf or .zip file.
#' @param quiet Sets the system call flag for a quiet (non-verbose) run.
#' @param clean Sets the system call for cleaning up the folder after completion.
#'   If \code{TRUE}, extra files generated during processing will be deleted.
#' @param replace when \code{pdf.zip} is not \code{NULL}, this determines if 
#'   the new file will overwrite an existing file.
#' @param ... Additoinal arguments to be passed to \code{tools::texi2dvi}
#' 
#' @details For TEX files, a call is made using \code{tools::texi2dvi} to 
#'   compile a PDF file.
#'   
#'   For HTML files, the referenced figures are gathered, copied into a
#'   subdirectory and the HTML document and the figures are placed into a 
#'   zip folder for convenient transfer.  All of the image links in the HTML
#'   code are modified to reflect the new location.  Lastly, a text file 
#'   is added with instructions for unzipping the files for convenient viewing
#'   (but don't worry, no one ever reads this).
#'   
#' @author Benjamin Nutter
#' 
lazy.build <- function(filename, pdf.zip=NULL, quiet=TRUE, clean=TRUE, 
                       replace=TRUE, ...){
  #*** retrieve the report format
  reportFormat <- getOption("lazyReportFormat")
  if (!reportFormat %in% c("latex", "html")) stop("option(\"lazyReportFormat\") must be either 'latex' or 'html'")
  
  #*** the following block changes the extension of 'filename' to match options("lazyReportFormat")
  file <- unlist(strsplit(filename, "[.]"))
  fileout <- if (!is.null(pdf.zip)) unlist(strsplit(pdf.zip, "[.]")) else NULL
  file.ext <- utils::tail(file, 1)
  if (reportFormat == "latex" && file.ext %in% c("html", "htm")){ 
    filename <- paste(c(file[-length(file)], "tex"), collapse=".")
    pdf.zip <- if (!is.null(pdf.zip)) paste(c(fileout[-length(fileout)], "pdf"), collapse=".") else NULL
  }
  if (reportFormat == "html" && file.ext == "tex"){
    filename <- paste(c(file[-length(file)], "html"), collapse=".")
    pdf.zip <- if (!is.null(pdf.zip)) paste(c(fileout[-length(fileout)], "zip"), collapse=".") else NULL
  }
  

  #*** Build the PDF and copy to new location if requested.
  if (reportFormat == "latex"){
    #*** The next three blocks generate objects for the file paths of the .tex file and the .pdf file
    path.comp <- unlist(strsplit(filename, .Platform$file.sep))
    path <- paste(path.comp[-length(path.comp)], collapse=.Platform$file.sep)
    if (path == "") path <- getwd()

    outfile <- gsub("[.]tex", ".pdf", path.comp[length(path.comp)])
    outfile <- paste(getwd(), outfile, sep=.Platform$file.sep)
                                
  
    if (is.null(pdf.zip)) pdf.zip <- gsub("[.]tex", ".pdf", filename)
    path.pdf <- unlist(strsplit(pdf.zip, .Platform$file.sep))
    path.pdf <- paste(path.pdf[-length(path.pdf)], collapse=.Platform$file.sep)
    if (path.pdf == "") path.pdf <- getwd()
    
    tools::texi2dvi(filename, pdf=TRUE, quiet=quiet, clean=clean, ...)
   
    if (!(pdf.zip %in% path)){
      if(path.pdf != getwd()){
        file.copy(outfile, pdf.zip, overwrite=replace)
        file.remove(outfile)
      }
    }
  }
  
  
  if (reportFormat == "html"){
    if (is.null(pdf.zip)){ 
      pdf.zip <- unlist(strsplit(filename, "[.]"))
      pdf.zip <- paste(c(pdf.zip[-length(pdf.zip)], "zip"), collapse=".")
    }
    
    #*** Remove extension from pdf.zip
    pdf.zip <- gsub("[.]zip", "", pdf.zip)
    
    #*** designate the filename for the HTML file
    filename.to <- sapply(strsplit(filename, .Platform$file.sep), utils::tail, 1)
    
    #*** read the HTML code and identify rows with figures.
    code <- readLines(filename)
    figure.pos <- grep("[<]img src", code)
    
    #*** identify figure files
    figures <- code[figure.pos]
    figures <- gsub("[[:print:]]+[<]img src='", "", figures)
    figures <- gsub("' height[[:print:]]+", "", figures)
    
    #*** create file paths for copies of figures
    weblinks <- grep("http[://]", figures)
    #igures.to <- figures
    figures.to <- sapply(strsplit(figures, .Platform$file.sep), utils::tail, 1)  
    
    if (length(weblinks) > 0) figures.to[weblinks] <- figures
    
    if (length(figure.pos) > 0)
      for(i in 1:length(figure.pos)) code[figure.pos[i]] <- gsub(figures[i], figures.to[i], code[figure.pos[i]])
    
    #*** if the HTML file exists in the working directory, copy the file to prevent deletion
    if (file.exists(filename.to)){
      rename <- TRUE
      file.rename(filename, "Temporary_Name_Change.html")
    }
    else rename <- FALSE
    
    #*** Rewrite the HTML with abbreviated figure file paths
    write(code, file.path(getwd(), filename.to))
    
    #*** Copy files to the working directory
    if (length(weblinks) > 0) 
      if (length(figures.to[-weblinks]) > 0) file.copy(figures[-weblinks], file.path(getwd(), figures.to[-weblinks]))
    else file.copy(figures, file.path(getwd(), figures.to))
    
    #*** Troubleshoot file
    write(paste("Depending on your system, you may need to unzip/extract all of the",
                "files in the .zip folder in order to view images and figures.",
                "",
                "",
                "This HTML file may not display as intended in some versions of Internet Explorer.",
                "If available, it is recommend that you open the HTML files in Mozilla Firefox or Google Chrome.",
                "If these are not available, you may right-click on the file and open in MS Word.", sep="\n"),
          "00_Troubleshooting.txt")
    
    #*** Zip the files and clean up
    if (length(figures.to) > 0)
      utils::zip(pdf.zip, c(if (length(weblinks) > 0) figures.to[-weblinks] else figures.to, filename.to, "00_Troubleshooting.txt"))
    else 
      utils::zip(pdf.zip, c(filename.to, "00_Troubleshooting.txt"))
    unlink(c(filename.to, figures.to, "00_Troubleshooting.txt"))
    if (rename) file.rename("Temporary_Name_Change.html", filename)  
  }
}
