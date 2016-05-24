#' @title Render and Export DOT Graphs in R
#' @usage dot(DOT, file = NULL, return = "auto")
#'
#' @keywords graphics plot diagram literate visualization
#' @description Graph Description Language (DOT) is a simplified and intuitive plain text graphical language. The \code{dot()} function renders the DOT markup language in R and also provides the possibility to export the graphs in PostScript and SVG (Scalable Vector Graphics) formats. In addition, it supports literate programming packages such as Knitr and R2HTML. Visit \url{http://haghish.com/dot} for downloading examples of creating algorithms and several graphs for \code{Rmarkdown} and \code{R HTML} to generate dynamic procedural graphs in dynamic documents using the \code{DOT} package.
#'
#' @param DOT This argument takes a string containing the DOT script. It is advised to use single quotation mark for the DOT string since the script often includes double quotations which can crash the function.
#'
#' @param file defines the file name for exporting the graph. The acceptable file extensions are \code{"ps"} for PostScript and \code{"svg"} for SVG format (see examples below).
#'
#' @param return specifies if PS or SVG script should be printed in R console. The acceptable values are \code{"auto"}, \code{"cat"}, \code{"verbatim"}, and \code{NULL}. The default value is \code{"auto"} which does not return anything unless the \code{dot()} function is called within \code{'knitr'} or \code{'rmarkdown'} packages, which require concatenated graphical script. The \code{"cat"} returns concatenated PS or SVG script, printed in multiple lines in the R console. The \code{"verbatim"} returns a single string which is assignable to an object.
#'
#' @return By default, the function only renders and loads the DOT plot in RStudio but does not return any PS or SVG script. If the \code{return} argument is specified, it returns PostScript or SVG script. Note that for assigning the script returned from \code{dot()} to an object, only \code{"verbatim"} value can be used to create a string object.
#'
#' @author E. F. Haghish \cr
#' Medical Informatics and Biostatistics (IMBI) \cr
#' University of Freiburg, Germany \cr
#' \email{haghish@imbi.uni-freiburg.de} \cr
#' \cr
#' Department of Mathematics and Computer Science \cr
#' University of Southern Denmark \cr
#' \email{haghish@imada.sdu.dk}
#'
#' @examples

#' #create a simple DOT graph and load it in RStudio
#' dot("digraph {A -> B;}")
#'
#' #to produce a dynamic document including a diagram in 'rmarkdown'
#' \dontrun{
#' ```{r echo=FALSE, results='asis'}
#' library(DOT)
#' dot("digraph {A -> B;}", return = "cat")
#' ```}
#'
#'
#' #create a DOT graph and export a SVG file to the working directory
#' dot("digraph {A -> B; B -> C; B -> D;}", file = "myfile.svg")
#'
#' #export the example above in PostScript format
#' dot("digraph {A -> B; B -> C; B -> D;}", file = "myfile.ps")
#'
#' #create a DOT graph and save the script in a string object in R
#' myString <- dot("digraph {A -> B;}", return = "verbatim")



#' @export
#' @import V8
#' @importFrom tools file_ext
#' @importFrom utils packageVersion

dot <- function(DOT, file = NULL, return = "auto") {

    # ---------------------------------------------------------
    # SYNTAX PROCESSING
    #   - create temporary directory
    #   - obtain the full path to the JS library in the package
    #   - check the file extension
    #   - load V8 if there is a demand for exporting graphics
    # =========================================================

    if(!requireNamespace("tools", quietly = TRUE)) stop("tools package is required.", call. = FALSE)
    #requireNamespace("tools", quietly = TRUE)

    extension <- ""
    autoReturn <- ""

    if (!is.null(file)) {
        extension <- file_ext(file)
    }

    if (!is.null(file) & extension != "svg" & extension != "ps") {
        err <- paste("\n", extension, "is not a valid file extension.",
            "DOT recognizes 'svg' and 'ps' file extensions only" , sep = " ")
        stop(err)
    }


    if(!requireNamespace("V8", quietly = TRUE)) stop("V8 package is required", call. = FALSE)
    #requireNamespace("V8", quietly = TRUE)
    if (!(packageVersion("V8") >= "1.0")) {
        stop('your "V8" package is too old. The DOT package requires V8
        package version 1.0 or newer', call. = FALSE)
    }

    # ---------------------------------------------------------
    # CHECK THE PARENT FUNCTION
    #   - is dot() called by 'knitr' or 'rmarkdown' ?
    # =========================================================
    if (return == "auto") {
        parentName <- try(deparse(sys.calls()[[1]]), silent = TRUE)
        PN <- try(unlist(strsplit(parentName, "::"))[[1]], silent = TRUE)
        if (PN == "rmarkdown") {
            autoReturn = "cat"
        }
        else {
            PN <- substring(parentName, 1, 4)[[1]]
            if (PN == "knit") {
                autoReturn = "cat"
            }
        }
    }


    # ---------------------------------------------------------
    # DOT PROCESSING
    #   - remove the "\n" from DOT source
    #   - detecting DOT engine (digraph vs graph)
    # =========================================================
    scriptDOT <- gsub("[\r\n]", " ", DOT)
    List <- strsplit(DOT, " ")[[1]]
    engine <- List[[1]]
    if (engine == "graph") {
        DOT <- paste("'", scriptDOT, "'", ', { engine: "neato" }', sep="")
    }
    else {
        DOT <- paste("'", scriptDOT, "'", sep="")
    }

    # ---------------------------------------------------------
    # CREATING THE GRAPH
    #   - obtain the full path to the JS library in the package
    #   - create temporary directory and temporary SVG file
    #   - obtain the full path to the JS library in the package
    #   - create the dot.svg file
    #   - load in the viewer or web-browser
    # =========================================================


    tempDir <- tempfile()
    dir.create(tempDir)
    dotFile <- file.path(tempDir, "dot.svg")
    file.create(dotFile)
    JS <- system.file("lib/viz.js", package = "DOT")
    ENV <- v8();
    ENV$source(JS)
    call <- paste("image = Viz(", DOT, ")", sep = "")
    content <- ENV$eval(call)
    #print(content)

    #call2 <- paste("Viz.svgXmlToPngImageElement(Viz(", DOT, "))", sep = "")
    #content2 <- ENV$eval(call2)

    #print(content2)

    write(content, file=dotFile, append=TRUE)

    viewer <- getOption("viewer")
    if (!is.null(viewer)) {
        viewer(dotFile)
    }
    #else {
    #    utils::browseURL(dotFile)
    #}

    svgContent <- content

    # ---------------------------------------------------------
    # EXPORTING THE GRAPH
    # =========================================================
    if (!is.null(file)) {
        if (extension == "ps") {
            scriptDOT <- paste("'", scriptDOT, "'", ' , { format: "ps" }', sep = "")
            if (engine == "graph") {
                DOT <- paste(scriptDOT, ' , { engine: "neato" }', sep="")
            }
            else {
                DOT <- paste(scriptDOT, sep="")
            }
            call <- paste("Viz(", DOT, ")", sep = "")
            content <- ENV$eval(call)
        }

        export <- file.path(getwd(), file)
        file.create(export)
        write(content, file=export, append=TRUE)
    }

    if (!is.null(return)) {
        if (return == "cat") {
            cat(content)
        }
        else if (return == "verbatim") {
            return(content)
        }
        else if (return == "auto" & autoReturn == "cat") {
            cat(svgContent)       #always return SVG
        }
    }
}






