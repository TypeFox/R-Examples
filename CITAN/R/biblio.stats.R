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


#' Performs preliminary analysis of data in a Local Bibliometric Storage
#' by creating some basic descriptive statistics (numeric and graphical).
#' Dataset may be restricted to any given document types
#' or a single survey.
#'
#' Plot types (accessed with \code{which}):
#' \itemize{
#' 	\item \code{1} --- "Document types",
#' 	\item \code{2} --- "Publication years",
#' 	\item \code{3} --- "Citations per document",
#' 	\item \code{4} --- "Citations of cited documents per type",
#' 	\item \code{5} --- "Number of pages per document type",
#' 	\item \code{6} --- "Categories of documents" (based od source categories),
#' 	\item \code{7} --- "Documents per author".
#' }
#'
#' Note that this user interaction scheme is similar in behavior
#' to the \code{\link{plot.lm}} function.
#'
#' @title Perform preliminary analysis of data in a Local Bibliometric Storage
#' @param conn connection object, see \code{\link{lbsConnect}}.
#' @param documentTypes character vector or \code{NULL}; specifies document types to restrict to;
#'    a combination of \code{Article}, \code{Article in Press}, \code{Book}, \code{Conference Paper},
#'    \code{Editorial}, \code{Erratum}, \code{Letter}, \code{Note}, \code{Report}, \code{Review},
#'    \code{Short Survey}. \code{NULL} means no restriction.
#' @param surveyDescription single character string or \code{NULL}; survey to restrict to, or \code{NULL} for no restriction.
#' @param which numeric vector with elements in 1,...,7, or \code{NULL}; plot types to be displayed.
#' @param main title for each plot.
#' @param ask logical; if \code{TRUE}, the user is asked to press return before each plot.
#' @param ... additional graphical parameters, see \code{\link{plot.default}}.
#' @param cex.caption controls size of default captions.
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' lbsDescriptiveStats(conn, surveyDescription="Scientometrics",
#'    documentTypes=c("Article", "Note", "Report", "Review", "Short Survey"));
#' ## ...
#' lbsDisconnect(conn);}
#'
#' @export
#'
#' @seealso
#' \code{\link{plot.default}},
#' \code{\link{lbsConnect}}
lbsDescriptiveStats <- function(conn,
   documentTypes=NULL,
   surveyDescription=NULL,
   which=(1L:7L),
   main="",
   ask = (prod(par("mfcol")) < length(which) && dev.interactive()),
   ...,
   cex.caption=1
)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   ## ---------------------- AUXILIARY FUNCTION -------------------------------

   #' /internal/
   .lbsDescriptiveStats_PrintStatsSubset <- function(conn, subQueryWhere)
   {
      query <- sprintf("
         SELECT COUNT(idDocument)
         FROM
         (
            SELECT DISTINCT Biblio_Documents.IdDocument
            FROM Biblio_Documents
            JOIN ViewBiblio_DocumentsSurveys ON ViewBiblio_DocumentsSurveys.IdDocument=Biblio_Documents.IdDocument
            WHERE %s
         );",
         subQueryWhere);
      res <- dbGetQuery(conn, query);
      cat(sprintf("Number of documents in the selected subset: %g.\n", res[1,1]));


      query <- sprintf("
         SELECT COUNT(IdAuthor)
         FROM
         (
            SELECT DISTINCT Biblio_AuthorsDocuments.IdAuthor
            FROM Biblio_AuthorsDocuments
            JOIN
            (
               SELECT Biblio_Documents.IdDocument
               FROM Biblio_Documents
               JOIN ViewBiblio_DocumentsSurveys ON ViewBiblio_DocumentsSurveys.IdDocument=Biblio_Documents.IdDocument
               WHERE %s
            ) AS Docs ON Docs.IdDocument=Biblio_AuthorsDocuments.IdDocument
         );",
         subQueryWhere);
      res <- dbGetQuery(conn, query);
      cat(sprintf("Number of authors in the selected subset:   %g.\n", res[1,1]));

      cat("\n");
   }

   ## ---------------------- AUXILIARY FUNCTION -------------------------------

   #' /internal/
   .lbsDescriptiveStats_PrintSurveyStats <- function(conn, subQueryWhere, split2filenames)
   {
      if (!split2filenames)
      {
         query <- sprintf("
            SELECT Description AS surveyDescription, COUNT(IdDocument) AS DocumentCount
            FROM
            (
               SELECT DISTINCT Description, ViewBiblio_DocumentsSurveys.IdDocument
               FROM ViewBiblio_DocumentsSurveys
               JOIN Biblio_Documents ON ViewBiblio_DocumentsSurveys.IdDocument=Biblio_Documents.IdDocument
               WHERE %s
            ) GROUP BY Description;",
            subQueryWhere
         );
         res <- dbGetQuery(conn, query);

         cat("Surveys:\n");
         print(res);
      } else
      {
         query <- sprintf("
            SELECT Filename, Timestamp, COUNT(IdDocument) AS DocumentCount
            FROM
            (
               SELECT DISTINCT Filename, Timestamp, ViewBiblio_DocumentsSurveys.IdDocument
               FROM ViewBiblio_DocumentsSurveys
               JOIN Biblio_Documents ON ViewBiblio_DocumentsSurveys.IdDocument=Biblio_Documents.IdDocument
               WHERE %s
            ) GROUP BY Filename;",
            subQueryWhere
         );
         res <- dbGetQuery(conn, query);

         cat("Source files in selected survey:\n");
         print(res);
      }

      cat("  * Note that a document may belong to many surveys/files.\n");
      cat("\n");
   }


   ## -----------------------------------------------------
   ## Basic stats

   res <- dbGetQuery(conn, "SELECT COUNT(idSource) FROM Biblio_Sources;");
   cat(sprintf("Number of sources in your LBS:           %g\n", res[1,1]));

   res <- dbGetQuery(conn, "SELECT COUNT(idDocument) FROM Biblio_Documents;");
   cat(sprintf("Number of documents in your LBS:         %g\n", res[1,1]));

   res <- dbGetQuery(conn, "SELECT COUNT(idAuthor) FROM Biblio_Authors;");
   cat(sprintf("Number of author records in your LBS:    %g\n", res[1,1]));

   res <- dbGetQuery(conn, "SELECT COUNT(*) FROM (SELECT AuthorGroup FROM Biblio_Authors GROUP BY AuthorGroup);");
   cat(sprintf("Number of author groups in your LBS:     %g\n", res[1,1]));

   res <- dbGetQuery(conn, "SELECT COUNT(*) FROM (SELECT IdAuthor FROM Biblio_Authors WHERE AuthorGroup IS NULL);");
   cat(sprintf("Number of ungrouped authors in your LBS: %g\n", res[1,1]));

   cat("\n");


   ## -----------------------------------------------------
   ## Data set restrictions & subset stats

   surveyDescription  <- .lbs_PrepareRestriction_SurveyDescription(conn, surveyDescription);
   documentTypesShort <- .lbs_PrepareRestriction_DocumentTypes(conn, documentTypes);



   # Get subQueryWhere
   if (length(documentTypesShort)>0)
   {
      subQueryWhere <- sprintf("(%s)",
         paste("Type", documentTypesShort, sep="=", collapse=" OR "));
   } else subQueryWhere <- "1";

   if (!is.null(surveyDescription))
      subQueryWhere <- paste(c(subQueryWhere, sprintf(" Description='%s'", surveyDescription)), collapse=" AND ");



   cat("You have chosen the following data restrictions:\n");
   cat(sprintf("\tSurvey:         %s.\n", ifelse(is.null(surveyDescription), "<ALL>", surveyDescription)));
   cat(sprintf("\tDocument types: %s.\n", ifelse(is.null(documentTypesShort), "<ALL>", paste(documentTypesShort, collapse=", "))));
   cat("\n");



   if (length(documentTypesShort)>0 || !is.null(surveyDescription))
      .lbsDescriptiveStats_PrintStatsSubset(conn, subQueryWhere);



   ## -----------------------------------------------------
   ## Survey(s) stats

   .lbsDescriptiveStats_PrintSurveyStats(conn, subQueryWhere, !is.null(surveyDescription));

   which <- as.integer(which);
   which <- which[which >=1 & which <= 7];
   if (length(which) < 1) return();

   ## -----------------------------------------------------
   ## UI
   # This user interaction scheme is based on the code of plot.lm() function

   show <- rep(FALSE, 7)
   show[which] <- TRUE

   if (ask) {
      oask <- devAskNewPage(TRUE);
      on.exit(devAskNewPage(oask));
   }

   captions <- c("Document types", "Publication years", "Citations per document",
      "Citations of cited documents per type", "Number of pages per document type",
      "Categories of documents", "Documents per author");



   ## -----------------------------------------------------
   ## 1:5

   if (any(show[1L:5L]))
   {
      res <- dbGetQuery(conn, sprintf("
         SELECT DISTINCT Biblio_Documents.IdDocument, Citations, Type, Year, Pages
         FROM Biblio_Documents
         JOIN ViewBiblio_DocumentsSurveys ON Biblio_Documents.IdDocument=ViewBiblio_DocumentsSurveys.IdDocument
         WHERE %s",
         subQueryWhere));
   } else res <- NULL;



   if (show[1L])
   {
      tab <- sort(table(res$Type), decreasing=TRUE);
      barplot(tab, main=main, ...);
      mtext(as.graphicsAnnot(captions[1]), 3, 0.25, cex=cex.caption);

      cat(sprintf("%s:\n", captions[1]));
      print(tab);
      cat("\n\n");
   }

   if (show[2L])
   {
      tab <- table(res$Year);
      barplot(tab, main=main, ...);
      mtext(as.graphicsAnnot(captions[2]), 3, 0.25, cex=cex.caption);

      cat(sprintf("%s:\n", captions[2]));
      print(tab);
      cat("\n\n");
   }

   if (show[3L])
   {
      tab <- table(res$Citations);
      barplot(tab, main=main, ...);
      mtext(as.graphicsAnnot(captions[3]), 3, 0.25, cex=cex.caption);

      cat(sprintf("%s:\n", captions[3]));
      print(tab);
      cat("\n\n");
   }

   if (show[4L])
   {
      boxplot(res$Citations[res$Citations>0]~res$Type[res$Citations>0], log="y", main=main, ...);
      mtext(as.graphicsAnnot(captions[4]), 3, 0.25, cex=cex.caption);
   }

   if (show[5L])
   {
      boxplot(res$Pages[res$Pages>0 & res$Pages<500]~res$Type[res$Pages>0 & res$Pages<500], log="y", main=main, ...);
      mtext(as.graphicsAnnot(captions[5]), 3, 0.25, cex=cex.caption);
   }



   # -----------------------------------------------------
   # 6

   if (any(show[6L:6L]))
   {
      query <- sprintf(
         "SELECT DISTINCT DocInfo.IdDocument, IdCategoryParent, DescriptionParent
         FROM ViewBiblio_DocumentsCategories
         JOIN
         (
            SELECT DISTINCT Biblio_Documents.IdDocument FROM Biblio_Documents
            JOIN ViewBiblio_DocumentsSurveys ON Biblio_Documents.IdDocument=ViewBiblio_DocumentsSurveys.IdDocument
            WHERE %s
         ) AS DocInfo ON DocInfo.IdDocument=ViewBiblio_DocumentsCategories.IdDocument;",
         subQueryWhere
      );
      res <- dbGetQuery(conn, query)
   } else res <- NULL;

   if (show[6L])
   {
      mergepercent <- 0.017;
      tab <- table(as.factor(res$DescriptionParent));

      tab2 <- tab[tab>mergepercent*sum(tab)];
      tab2 <- c(tab2, "Other"=sum(tab[tab<=mergepercent*sum(tab)]));

      otab <- order(tab2);
      for (i in 1:(length(otab)/2))
      {
         if (i %% 2 == 0)
         {
            r <- otab[length(otab)-i/2];
            otab[length(otab)-i/2] <- otab[i];
            otab[i] <- r;
         }
      }
      tab <- tab2[otab];

      pie(tab, main=main, ...);
      mtext(as.graphicsAnnot(captions[6]), 3, 0.25, cex=cex.caption);

      cat(sprintf("%s:\n", captions[6]));
      print(tab);
      cat("\n\n");
   }




   # -----------------------------------------------------
   # 7

   if (any(show[7L:7L]))
   {
      query <- sprintf(
         "SELECT IdAuthor, COUNT(IdDocument) AS Number
         FROM
         (
            SELECT DISTINCT Biblio_AuthorsDocuments.IdAuthor, Biblio_AuthorsDocuments.IdDocument
            FROM Biblio_AuthorsDocuments
            JOIN Biblio_Documents ON Biblio_Documents.IdDocument=Biblio_AuthorsDocuments.IdDocument
            JOIN ViewBiblio_DocumentsSurveys ON ViewBiblio_DocumentsSurveys.IdDocument=Biblio_Documents.IdDocument
            WHERE %s
         )
         GROUP BY (IdAuthor)",
         subQueryWhere
      );
      res <- dbGetQuery(conn, query)
   } else res <- NULL;

   if (show[7L])
   {
      tab <- table(res$Number);
      barplot(tab, main=main, ...);
      mtext(as.graphicsAnnot(captions[7]), 3, 0.25, cex=cex.caption);

      cat(sprintf("%s:\n", captions[7]));
      print(tab);
      cat("\n\n");
   }
}
