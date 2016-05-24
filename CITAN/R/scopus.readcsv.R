## This file is part of the CITAN package for R
##
## Copyright 2011-2015 Marek Gagolewski
##
##
## CITAN is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## CITAN is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with CITAN. If not, see <http://www.gnu.org/licenses/>.


#' Reads bibliography entries from a UTF-8 encoded CSV file.
#'
#'
#'
#' The  \code{\link{read.csv}} function is used to read the bibliography.
#' You may therefore freely modify its behavior
#' by passing further arguments (\code{...}), see the manual page
#' of \code{\link{read.table}} for details.
#'
#' The CSV file should consist at least of the following columns.
#' \enumerate{
#' \item \code{Authors}: Author name(s) (surname first; multiple names are comma-separated,
#' e.g. \dQuote{Smith John, Nowak G. W.}),
#' \item \code{Title}: Document title,
#' \item \code{Year}: Year of publication,
#' \item \code{Source.title}: Source title, e.g. journal name,
#' \item \code{Volume}: Volume number,
#' \item \code{Issue}: Issue number,
#' \item \code{Page.start}: Start page number,
#' \item \code{Page.end}: End page number,
#' \item \code{Cited.by}: Number of citations received,
#' \item \code{Link}: String containing unique document identifier, by default of the form ...id=\emph{\strong{UNIQUE_ID}}&... (see \code{alternativeIdPattern} parameter),
#' \item \code{Document.Type}: Document type, one of: \dQuote{Article}, \dQuote{Article in Press},
#'        \dQuote{Book}, \dQuote{Conference Paper}, \dQuote{Editorial},
#'        \dQuote{Erratum}, \dQuote{Letter}, \dQuote{Note}, \dQuote{Report},
#'        \dQuote{Review}, \dQuote{Short Survey}, or \code{NA}
#'        (other categories are treated as \code{NA}s),
#' \item \code{Source}: Data source identifier, must be the same as the
#'        \code{dbIdentifier} parameter value. It is used for parse errors detection.
#' }
#'
#' The CSV file to be read may, for example, be created by \emph{SciVerse Scopus}
#' (Export format=\emph{comma separated file, .csv (e.g. Excel)},
#' Output=\emph{Complete format} or \emph{Citations only}).
#' Note that the exported CSV file sometimes needs to be corrected by hand
#' (wrong page numbers, single double quotes in character strings instead of two-double quotes etc.).
#' We suggest to make the corrections in a \dQuote{Notepad}-like application
#' (in plain text). The function tries to indicate line numbers causing
#' potential problems.
#'
#' @title Import bibliography entries from a CSV file.
#' @param filename the name of the file which the data are to be read from, see \code{\link{read.csv}}.
#' @param stopOnErrors logical; \code{TRUE} to stop on all potential parse errors or just warn otherwise.
#' @param dbIdentifier character or \code{NA}; database identifier, helps detect parse errors, see above.
#' @param alternativeIdPattern character; regular expression used to extract AlternativeId, \code{NA} to get the id as is,
#' @param ... further arguments to be passed to \code{read.csv}.
#' @return A \code{data.frame} containing the following 11 columns:
#' \tabular{ll}{
#' \code{Authors} \tab	Author name(s), comma-separated, surnames first.\cr
#' \code{Title} \tab	Document title.\cr
#' \code{Year} \tab	Year of publication.\cr
#' \code{AlternativeId} \tab	Unique document identifier.\cr
#' \code{SourceTitle} \tab	Title of the source containing the document.\cr
#' \code{Volume} \tab	Volume.\cr
#' \code{Issue} \tab	Issue.\cr
#' \code{PageStart} \tab	Start page; numeric.\cr
#' \code{PageEnd} \tab	End page; numeric.\cr
#' \code{Citations} \tab	Number of citations; numeric.\cr
#' \code{DocumentType} \tab	Type of the document; see above.\cr
#' }
#' The object returned may be imported into a local bibliometric storage via \code{\link{lbsImportDocuments}}.
#' @export
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' data <- Scopus_ReadCSV("db_Polish_MATH/Poland_MATH_1987-1993.csv");
#' lbsImportDocuments(conn, data, "Poland_MATH");
#' ## ...
#' lbsDisconnect(conn);}
#' @seealso \code{\link{Scopus_ASJC}}, \code{\link{Scopus_SourceList}},
#' \code{\link{lbsConnect}},
#' \code{\link{Scopus_ImportSources}},\cr
#' \code{\link{read.table}}, \code{\link{lbsImportDocuments}}
Scopus_ReadCSV <- function(filename, stopOnErrors=TRUE, dbIdentifier='Scopus', alternativeIdPattern="^.*\\id=|\\&.*$", ...)
{
   datafile <- read.csv(filename, header = T, encoding="UTF-8", fileEncoding="UTF-8", stringsAsFactors=FALSE, ...);

   if (!is.na(dbIdentifier) && is.null(datafile$Source)) stop("Column not found: `Source'.");
   if (is.null(datafile$Authors)) stop("Column not found: `Authors'.");
   if (is.null(datafile$Title)) stop("Column not found: `Title'.");
   if (is.null(datafile$Year)) stop("Column not found: `Year'.");
   if (is.null(datafile$Source.title)) stop("Column not found: `Source.title'.");
   if (is.null(datafile$Volume)) stop("Column not found: `Volume'.");
   if (is.null(datafile$Issue)) stop("Column not found: `Issue'.");
   if (is.null(datafile$Page.start)) stop("Column not found: `Page.start'.");
   if (is.null(datafile$Page.end)) stop("Column not found: `Page.end'.");
   if (is.null(datafile$Cited.by)) stop("Column not found: `Cited.by'.");
   if (is.null(datafile$Link)) stop("Column not found: `Link'.");
   if (is.null(datafile$Document.Type)) stop("Column not found: `Document.Type'.");


   if (!is.na(dbIdentifier) && any(datafile$Source != dbIdentifier))
   {
      msg <- (sprintf("source database does not match 'dbIdentifier'. This may possibly indicate a parse error. Check records: %s.",
         paste(which(datafile$Source != dbIdentifier), collapse=", ")));

      if (stopOnErrors) stop(msg) else warning(msg);
   }


   if (!is.na(alternativeIdPattern))
   {
      datafile$AlternativeId <- gsub(alternativeIdPattern, "", datafile$Link); # REG EXP
   } else {
      datafile$AlternativeId <- datafile$Link; # AS IS
   }


   naAlternativeId <- which(is.na(datafile$AlternativeId));
   if (length(naAlternativeId) > 0)
   {
      msg <- (sprintf("some documents do not have unique identifiers. Check line %s (or its neighborhood). \
   Perhaps somethings is wrong with the end page (check for ', ' nearby).",
         naAlternativeId[1]+1));

      if (stopOnErrors) stop(msg) else warning(msg);
   }

   checkAlternativeId <- unique(datafile$AlternativeId, incomparables=NA);
   if (length(checkAlternativeId) != nrow(datafile))
   {
      msg <- (sprintf("non-unique document identifiers at rows: %s.",
         paste((1:nrow(datafile))[-checkAlternativeId], collapse=", ")));

      if (stopOnErrors) stop(msg) else warning(msg);
   }


   datafile$Cited.by[!is.na(gsub("^([[:digit:]]+)$", NA, datafile$Cited.by))] <- NA;
   datafile$Cited.by <- as.numeric(datafile$Cited.by);
   checkCitations <- which(datafile$Cited.by < 0 | datafile$Cited.by>100000);
   if (length(checkCitations) > 0)
   {
      msg <- (sprintf("something is wrong with citation counts at rows: %s.",
         paste((1:nrow(datafile))[-checkCitations], collapse=", ")));

      if (stopOnErrors) stop(msg) else warning(msg);
   }


   datafile$Page.start[!is.na(gsub("^([[:digit:]]+)$", NA, datafile$Page.start))] <- NA;
   datafile$Page.end  [!is.na(gsub("^([[:digit:]]+)$", NA, datafile$Page.end  ))] <- NA;

   datafile$Page.start <- as.numeric(datafile$Page.start);
   datafile$Page.end   <- as.numeric(datafile$Page.end);

   checkPages <- which((datafile$Page.start<0) | (datafile$Page.end<datafile$Page.start) | (datafile$Page.end-datafile$Page.start>10000));
   if (length(checkPages) > 0)
   {
      msg <- (sprintf("some documents seem to have incorrect page numbers. Check line %s (or its neighborhood).",
         checkPages[1]+1));

      if (stopOnErrors) stop(msg) else warning(msg);
   }



   datafile <- data.frame(Authors=as.character(datafile$Authors),
                      Title=as.character(datafile$Title),
                      Year=datafile$Year,
                      AlternativeId=datafile$AlternativeId,
                      SourceTitle=as.character(datafile$Source.title),
                      Volume=datafile$Volume,
                      Issue=datafile$Issue,
                      PageStart=datafile$Page.start,
                      PageEnd=datafile$Page.end,
                      Citations=datafile$Cited.by,
                      DocumentType=datafile$Document.Type);

   attr(datafile, "filename") <- filename;
   return(datafile);
}
