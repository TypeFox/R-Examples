## Copyright 2013 Stefan Widgren and Maria Noremark,
## National Veterinary Institute, Sweden
##
## Licensed under the EUPL, Version 1.1 or - as soon they
## will be approved by the European Commission - subsequent
## versions of the EUPL (the "Licence");
## You may not use this work except in compliance with the
## Licence.
## You may obtain a copy of the Licence at:
##
## http://ec.europa.eu/idabc/eupl
##
## Unless required by applicable law or agreed to in
## writing, software distributed under the Licence is
## distributed on an "AS IS" basis,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
## express or implied.
## See the Licence for the specific language governing
## permissions and limitations under the Licence.

html_summary_table <- function(contacts, direction) {
    if(identical(direction, "in")) {
        arrow <- "&larr;"
    } else {
        arrow <- "&rarr;"
    }

    lines <- "<table border=0>"

    for(i in seq_len(nrow(contacts))) {
        lines <- c(lines, "<tr>")

        ## Pad with empty cells to the left
        lines <- c(lines,
                   paste(rep("<td>&nbsp;</td>", 2L*(contacts$distance[i]-1L)), collapse = ""),
                   sprintf('<td align="right">%s</td>', contacts$lhs[i]),
                   sprintf('<td align="right"><a href="#%s-%s-%s">%s</a></td>',
                           direction, contacts$lhs[i], contacts$rhs[i], arrow),
                   sprintf('<td align="right">%s</td>', contacts$rhs[i]))

        ## Pad with empty cells to the right
        lines <- c(lines,
                   paste(rep("<td>&nbsp;</td>",
                             2L*(max(contacts$distance-1L) - (contacts$distance[i]-1L))),
                         collapse = ""))

        lines <- c(lines, "</tr>")
    }

    lines <- c(lines, "</table>")

    lines
}

html_detailed_table <- function(contacts, direction) {
    if(identical(direction, "in")) {
        arrow <- "&larr;"
    } else {
        arrow <- "&rarr;"
    }

    contacts <- contacts[order(contacts$t, contacts$id, decreasing = FALSE),]

    contacts$id <- as.character(contacts$id)
    contacts$id[is.na(contacts$id)] <- "&nbsp;"

    contacts$n <- as.character(contacts$n)
    contacts$n[is.na(contacts$n)] <- "&nbsp;"

    contacts$category <- as.character(contacts$category)
    contacts$category[is.na(contacts$category)] <- "&nbsp;"

    as.character(unlist(by(contacts, sprintf("%s - %s", contacts$lhs, contacts$rhs), function(x) {
        lines <- "<p>"

        ## Create the lhs to/from rhs title and add the name of the href link
        lines <- c(lines,
                   sprintf('<h3><a name="%s-%s-%s">%s %s %s</a></h3>', direction,
                           x$lhs[1], x$rhs[1], x$lhs[1], arrow, x$rhs[1]),
                   "<table border=1>",
                   "<tr><th>Date</th><th>Id</th><th>N</th><th>Category</th><th>Source</th><th>Destination</th></tr>")

        for(i in seq_len(nrow(x))) {
            lines <- c(lines,
                       "<tr>",
                       sprintf('<td align="right">%s</td>', x$t[i]),
                       sprintf('<td align="right">%s</td>', x$id[i]),
                       sprintf('<td align="right">%s</td>', x$n[i]),
                       sprintf('<td align="right">%s</td>', x$category[i]),
                       sprintf('<td align="right">%s</td>', x$source[i]),
                       sprintf('<td align="right">%s</td>', x$destination[i]),
                       "</tr>")
        }

        lines <- c(lines,
                   "</table>",
                   "</p>")

        lines
    })))
}

html_report <- function(x)
{
    lines <- c("<html>",
               "<head>",
               sprintf("<title>%s</title>", x@root),
               "</head>",
               "<body>",
               "<h1 align='center'>Contact Tracing</h1>",
               sprintf("<h1 align='center'>Root: %s</h1>", x@root),
               "<h3 align='center'>EpiContactTrace</h3>",
               sprintf("<h3 align='center'>Version: %s</h3>", packageVersion("EpiContactTrace")),
               "<hr/>",
               "<h2>Summary ingoing contacts</h2>",
               "<p>",
               "<table border=0>",
               sprintf("<tr><td>In begin date:</td><td>%s</td></tr>", x@ingoingContacts@tBegin),
               sprintf("<tr><td>In end date:</td><td>%s</td></tr>", x@ingoingContacts@tEnd),
               sprintf("<tr><td>In days:</td><td>%i</td></tr>", x@ingoingContacts@tEnd - x@ingoingContacts@tBegin),
               sprintf("<tr><td>In degree:</td><td>%i</td></tr>", InDegree(x@ingoingContacts)),
               sprintf("<tr><td>Ingoing contact chain:</td><td>%i</td></tr>", IngoingContactChain(x@ingoingContacts)),
               "</table>",
               "</p>")

    if(length(x@ingoingContacts@source) > 0L) {
        ## Get network structure. The distance is used for indentation.
        df <- NetworkStructure(x@ingoingContacts)

        ## Add lhs and rhs, with respect to direction
        df$rhs <- df$source
        df$lhs <- df$destination

        lines <- c(lines, html_summary_table(df, "in"))
    } else {
        lines <- c(lines, "<p>No ingoing contacts during the search period.</p>")
    }

    lines <- c(lines,
               "<hr/>",
               "<h2>Summary outgoing contacts</h2>",
               "<p>",
               "<table border=0>",
               sprintf("<tr><td>Out begin date:</td><td>%s</td></tr>", x@outgoingContacts@tBegin),
               sprintf("<tr><td>Out end date:</td><td>%s</td></tr>", x@outgoingContacts@tEnd),
               sprintf("<tr><td>Out days:</td><td>%i</td></tr>", x@outgoingContacts@tEnd - x@outgoingContacts@tBegin),
               sprintf("<tr><td>Out degree:</td><td>%i</td></tr>", OutDegree(x@outgoingContacts)),
               sprintf("<tr><td>Outgoing contact chain:</td><td>%i</td></tr>", OutgoingContactChain(x@outgoingContacts)),
               "</table>",
               "</p>")

    if(length(x@outgoingContacts@source) > 0L) {
        ## Get network structure. The distance is used for indentation.
        df <- NetworkStructure(x@outgoingContacts)

        ## Add lhs and rhs, with respect to direction
        df$lhs <- df$source
        df$rhs <- df$destination

        lines <- c(lines, html_summary_table(df, "out"))
    } else {
        lines <- c(lines, "<p>No outgoing contacts during the search period.</p>")
    }

    if(length(x@ingoingContacts@source) > 0L) {
        lines <- c(lines,
                   "<hr/>",
                   "<h2>Direct ingoing contacts</h2>")

        df <- as(x@ingoingContacts, "data.frame")

        # Add lhs and rhs, with respect to direction
        df$lhs <- df$destination
        df$rhs <- df$source

        lines <- c(lines,
                   html_detailed_table(df[df$destination == df$root,], "in"))

        df <- df[df$destination != df$root,]
        if(nrow(df) > 0) {
            lines <- c(lines,
                       "<hr/>",
                       "<h2>Indirect ingoing contacts</h2>",
                       html_detailed_table(df, "in"))
        }
    }

    if(length(x@outgoingContacts@source) > 0L) {
        lines <- c(lines,
                   "<hr/>",
                   "<h2>Direct outgoing contacts</h2>")

        df <- as(x@outgoingContacts, "data.frame")

        ## Add lhs and rhs, with respect to direction
        df$lhs <- df$source
        df$rhs <- df$destination

        html_detailed_table(df[df$source == df$root,], "out")

        df <- df[df$source != df$root,]
        if(nrow(df) > 0) {
            lines <- c(lines,
                       "<hr/>",
                       "<h2>Indirect outgoing contacts</h2>",
                       html_detailed_table(df, "out"))
        }
    }

    lines <- c(lines,
               "<hr/>",
               sprintf("<h5>Generated %s by EpiContactTrace version %s</h5>\n",
                       Sys.time(), packageVersion("EpiContactTrace")),
               "</body>",
               "</html>")

    lines
}

##' Generate a contact tracing \code{Report}
##'
##' EpiContatTrace contains report templates to generate pdf- or html reports
##' for the farm specific contacts. These reports can be useful for hands-on
##' disease tracing in the field. The templates are used by Sweave and can be
##' adapted by the end user. However, in the default setting the report has the
##' following layout; first the contacts are visualised graphically in a plot,
##' as to give an immediate signal to the reader of the report of the number of
##' contacts. In the following, the contact data are presented with different
##' levels of detail split by ingoing and outgoing contacts. The first includes
##' collapsed data and the sequential contact structure at group level (i.e. no
##' information on individuals or dates). In this summary, the sequential
##' structure of each part of the chain is included, and a holding that appears
##' in several different parts of the chain can therefore be included more than
##' once in the summary. The reason for this is to facilitate sequential
##' tracing and getting an overview of each part of the chain. After the
##' summary all details of all contacts included in the contact chains is
##' presented, i.e. date of contact and data on individual level when
##' available. To generate pdf files a TeX installation must exist to compile the
##' latex file. The report is saved in the working directory with the name of
##' the root as filename.
##'
##'
##' @rdname Report-methods
##' @aliases Report Report-methods Report,ContactTrace-method
##' Report,list-method
##' @docType methods
##' @section Methods: \describe{
##'
##'   \item{\code{signature(object = "ContactTrace")}}{
##'     Generate a report for a \code{ContactTrace} object.
##'   }
##'
##'   \item{\code{signature(object = "list")}}{
##'     Generate reports for a list of \code{ContactTrace} objects.
##'   }
##' }
##' @param object the object
##' @param ... Additional arguments to the method
##' @param format the format to use, can be either 'html' or 'pdf'. The default
##'        is 'html'
##' @param template the Sweave template file to use. If none is provided, the default
##'        is used.
##' @references \itemize{
##'   \item Dube, C., et al., A review of network analysis terminology
##'     and its application to foot-and-mouth disease modelling and policy
##'     development. Transbound Emerg Dis 56 (2009) 73-85, doi:
##'     10.1111/j.1865-1682.2008.01064.x
##'
##'   \item Noremark, M., et al., Network analysis
##'     of cattle and pig movements in Sweden: Measures relevant for
##'     disease control and riskbased surveillance.  Preventive Veterinary
##'     Medicine 99 (2011) 78-90, doi: 10.1016/j.prevetmed.2010.12.009
##'
##'   \item Friedrich Leisch. Sweave: Dynamic generation of
##'     statistical reports using literate data analysis. In Wolfgang
##'     Hardle and Bernd Ronz, editors, Compstat 2002 - Proceedings in
##'     Computational Statistics, pages 575-580. Physica Verlag,
##'     Heidelberg, 2002. ISBN 3-7908-1517-9.
##'     \url{http://www.statistik.uni-muenchen.de/~leisch/Sweave/}
##' }
##' @seealso Sweave, texi2pdf.
##' @note Plots are not supported in version 0.8.6 since igraph0 has
##' been archived. We intend to resolve the issue in a future
##' version. Install version 0.8.5 and igraph0 manually from the
##' archive if plots are required. See section 6.3 in 'R Installation
##' and Administration' on how to install packages from source.
##' @keywords methods
##' @importFrom tools texi2pdf
##' @importFrom utils packageVersion
##' @importFrom utils Sweave
##' @examples
##' \dontrun{
##'
##' ## Load data
##' data(transfers)
##'
##' ## Perform contact tracing
##' contactTrace <- Trace(movements=transfers,
##'                       root=2645,
##'                       tEnd='2005-10-31',
##'                       days=90)
##'
##' ## Generate an html report showing details of the contact tracing for
##' ## root 2646.
##' ## Note: Creates the files 2645.html and 2645.png in the working
##' ## directory.
##' Report(contactTrace)
##'
##' ## It's possible to generate reports for a list of ContactTrace objects.
##' ## Perform contact tracing for ten of the included herds
##' root <- sort(unique(c(transfers$source, transfers$destination)))[1:10]
##'
##' ## Perform contact tracing
##' contactTrace <- Trace(movements=transfers,
##'                       root=root,
##'                       tEnd='2005-10-31',
##'                       days=90)
##'
##' ## Generate reports
##' ## Note: Creates the files 1.html, 2.html, ..., 10.html and
##' ## 1.png, 2.png, ..., 10.png in the working directory
##' Report(contactTrace)
##' }
##'
setGeneric("Report",
           signature = "object",
           function(object, ...) standardGeneric("Report"))

##' @rdname Report-methods
##' @export
setMethod("Report",
          signature(object = "ContactTrace"),
          function(object, format = c("html", "pdf"), template = NULL)
      {
          format <- match.arg(format)

          if (!is.null(.ct_env$ct)) {
              stop("Unable to create report. The ct object already exists")
          }

          ## Make sure the added object is removed.
          on.exit(.ct_env$ct <- NULL)

          ## Add the ContactTrace object to the .ct_env environment
          .ct_env$ct <- object

          if(identical(format, "html")) {
              writeLines(html_report(object), con = sprintf("%s.html", object@root))
          } else {
              if(is.null(template)) {
                  template <- system.file("Sweave/speak-latex.rnw", package="EpiContactTrace")
              }

              utils::Sweave(template, syntax="SweaveSyntaxNoweb")
              tools::texi2pdf(sub("rnw$", "tex", basename(template)), clean=TRUE)
              file.rename(sub("rnw$", "pdf", basename(template)), sprintf("%s.pdf", object@root))
              unlink(sub("rnw$", "tex", basename(template)))
          }

          invisible(NULL)
      }
)

##' @rdname Report-methods
##' @export
setMethod("Report",
          signature(object = "list"),
          function(object, format = c("html", "pdf"), template = NULL)
      {
          format <- match.arg(format)

          if(!all(sapply(object, function(x) length(x)) == 1)) {
              stop("Unexpected length of list")
          }

          if(!all(sapply(object, function(x) class(x)) == "ContactTrace")) {
              stop("Unexpected object in list")
          }

          lapply(object, function(x) Report(x,
                                            format=format,
                                            template=template))

          invisible(NULL)
      }
)

## In order to communicate with Sweave add an environment to store the
## current report object
.ct_env <- new.env()

##' Get current \code{ContactTrace} object when generating a report
##'
##' EpiContatTrace contains report templates to generate pdf- or html reports
##' for the farm specific contacts. These reports can be useful for hands-on
##' disease tracing in the field. The templates are used by Sweave and can be
##' adapted by the end user. This method enables communication of the current
##' \code{ContactTrace} object to the report.
##'
##'
##' @return The current \code{ContactTrace} object when generating a report
##' @export
ReportObject <- function()
{
    return(.ct_env$ct)
}
