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


#' Indicates, by finding similarities between documents' titles,
#'  groups of documents that possibly should be merged.
#'
#' The function determines fuzzy similarity measures of the titles. Its
#' specificity is controlled by the \code{aggressiveness} parameter.
#'
#' Search results are presented in a convenient-to-use graphical dialog box.
#' The function tries to order the groups of documents according
#' to their relevance (**EXPERIMENTAL** algorithm).
#' Note that the calculation often takes a few minutes!
#'
#' The \code{ignoreTitles.like} parameter determines search patterns in an SQL \code{LIKE} format,
#' i.e. an underscore \code{_} matches a single character and a percent sign
#' \code{\%} matches any set of characters. The search is case-insensitive.
#'
#' @title Find documents to be merged (**EXPERIMENTAL**)
#' @param conn connection object, see \code{\link{lbsConnect}}.
#' @param surveyDescription character string or \code{NULL}; survey description to restrict to or \code{NULL}.
#' @param ignoreTitles.like character vector of SQL-LIKE patterns to match documents' titles to be ignored or \code{NULL}.
#' @param aggressiveness nonnegative integer; \code{0} for showing only exact matches;
#' the higher the value, the more documents will be proposed.
#' @return
#' A numeric vector of user-selected documents' identifiers to be removed.
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' listdoc <- lbsFindDuplicateTitles(conn,
#'    ignoreTitles.like=c("\%In this issue\%", "\%Editorial", "\%Introduction",
#'    "Letter to \%", "\%Preface"),
#'    aggressiveness=2);
#' lbsDeleteDocuments(conn, listdoc);
#' dbCommit(conn);
#' ## ...}
#'
#' @export
#'
#' @seealso
#' \code{\link{lbsDeleteDocuments}},
#' \code{\link{lbsFindDuplicateAuthors}},
#' \code{\link{lbsGetInfoDocuments}}
lbsFindDuplicateTitles <- function(conn,
   surveyDescription=NULL,
   ignoreTitles.like=NULL,
   aggressiveness=1)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   ## ---------------- auxiliary function -----------------------------------------------------
   .lbsFindDuplicateTitles_getDups <- function(conn, surveyDescription, ignoreTitles.like, aggressiveness)
   {
      surveyDescription  <- .lbs_PrepareRestriction_SurveyDescription(conn, surveyDescription);

      cat("Looking for documents with duplicate titles... ");
      dups <- list();
      res <- NULL;

      if (!is.null(ignoreTitles.like))
         ignoreTitles.like <- toupper(ignoreTitles.like);

      if (aggressiveness == 0)
      {
         res <- dbGetQuery(conn, sprintf(
            "SELECT DupRes.IdDocument, DupRes.Title, DupRes.Type, COUNT(IdAuthor) AS AuthorCount
            FROM
            (
               SELECT DISTINCT Biblio_Documents.IdDocument, Title, Type, Citations
               FROM Biblio_Documents
               JOIN ViewBiblio_DocumentsSurveys ON Biblio_Documents.IdDocument=ViewBiblio_DocumentsSurveys.IdDocument
               WHERE %s AND %s AND Title IN (
                  SELECT Title FROM (
                     SELECT Title, COUNT(IdDocument) AS cnt FROM Biblio_Documents
                     GROUP BY Title
                  ) WHERE cnt > 1
               )
            ) AS DupRes
            JOIN Biblio_AuthorsDocuments ON DupRes.IdDocument=Biblio_AuthorsDocuments.IdDocument
            GROUP BY DupRes.IdDocument
            ORDER BY Title ASC, Citations DESC;",
            ifelse(is.null(surveyDescription), "1", sprintf(" Description='%s'", surveyDescription)),
            ifelse(is.null(ignoreTitles.like), "1", paste("NOT Title LIKE", sprintf("'%s'", ignoreTitles.like), collapse=" AND ", sep=" "))
         ));

         n <- nrow(res);
         if (n <= 1) return(dups);


         window <- .gtk2.progressBar(0, n, info="Looking for documents with duplicate titles...");

         res$Title <- toupper(res$Title);

         k <- 0;
         i <- 1;
         while (i <= n)
         {
            j <- i+1;
            while (j <= n && (res$Title[i] == res$Title[j]))
               j <- j+1;

            k <- k+1;
            dups[[k]] <- i:(j-1);

            .gtk2.progressBar(j-1,n,each=1,window=window);
            i <- j;
         }
      } else if (aggressiveness >= 1)
      {
         query <- sprintf("
            SELECT IdDocument, Title, Type, Citations, COUNT(IdAuthor) AS AuthorCount
            FROM
            (
               SELECT DISTINCT Biblio_Documents.IdDocument, Title, Type, Citations, IdAuthor
               FROM Biblio_Documents
               JOIN ViewBiblio_DocumentsSurveys  ON Biblio_Documents.IdDocument=ViewBiblio_DocumentsSurveys.IdDocument
               LEFT JOIN Biblio_AuthorsDocuments ON Biblio_Documents.IdDocument=Biblio_AuthorsDocuments.IdDocument
               WHERE %s AND %s
            )
            GROUP BY IdDocument
            ORDER BY Title ASC, Citations DESC;",
            ifelse(is.null(ignoreTitles.like), "1", paste("NOT Title LIKE", sprintf("'%s'", ignoreTitles.like), collapse=" AND ", sep=" ")),
            ifelse(is.null(surveyDescription), "1", sprintf(" Description='%s'", surveyDescription))
         );

         res <- dbGetQuery(conn, query);

         n <- nrow(res);

         if (n <= 1) return(dups);

         res$Title <- toupper(res$Title);

         if (aggressiveness == 1)
         {
            res$TitleComp <- (gsub("[^[:alpha:]]*", "", res$Title));
         } else
         {
            window <- .gtk2.progressBar(i, n, info="Preparing data...");
            res$TitleComp <- character(n);
            mtch <- gregexpr("[[:alpha:]]+", res$Title);
            for (i in 1:n)
            {
               m <- mtch[[i]];
               idx <- which(attr(m, "match.length") >= aggressiveness);
               if (length(idx) == 0)
               {
                  res$TitleComp[i] <- (gsub("[^[:alpha:]]*", "", res$Title[i]));
               } else
               {
                  res$TitleComp[i] <- toupper(paste(
                     sapply(idx, function(x) {
                        substr(res$Title[i], m[x], m[x]+attr(m, "match.length")[x]-1)
                     }),
                     sep="", collapse=""));
               }

               .gtk2.progressBar(i, n, window=window);
            }
         }


         res <- res[order(res$TitleComp),];


         window <- .gtk2.progressBar(0, n, info="Looking for documents with duplicate titles...");

         k <- 0;
         i <- 1;
         while (i <= n)
         {
            j <- i+1;
            while (j <= n && (
                  substr(res$TitleComp[i],1,min(nchar(res$TitleComp[i]),nchar(res$TitleComp[j]))) ==
                  substr(res$TitleComp[j],1,min(nchar(res$TitleComp[i]),nchar(res$TitleComp[j])))))
               j <- j+1;

            if (j-1 > i)
            {
               k <- k+1;
               dups[[k]] <- i:(j-1);
            }

            .gtk2.progressBar(j-1,n,each=1,window=window);
            i <- j;
         }
      }

      cat("DONE.\n");


      k <- length(dups);
      if (k == 0) return(dups);


      # Now, re-order dups according to some heuristic utility function
      # Assumptions: article-in-press++
      # letter--, editorial--, erratum--, note--
      # similar titles and number of authors+
      # similar titles but not number of authors-
      # many results---
      # take only two fist elements into account (they are sorted by title)
      utility <- rep(0, k);

      for (i in 1:k)
      {
         check <- dups[[i]];
         stopifnot(length(check)>=2);

         if (gsub("[^[:alpha:]]*", "", res$Title[check[1]]) == gsub("[^[:alpha:]]*", "", res$Title[check[2]]))
            utility[i] <- utility[i] + 1.0;

         utility[i] <- utility[i] + ifelse((res$AuthorCount[check[1]] == res$AuthorCount[check[2]]), 1.5, -2.3);
         if (all(!is.na(res$Type[check[1:2]])))
         {
            if (any(res$Type[check[1:2]]=="le") ||
                  any(res$Type[check[1:2]]=="no") ||
                  any(res$Type[check[1:2]]=="ed") ||
                  any(res$Type[check[1:2]]=="er"))
               utility[i] <- utility[i]-1.3;

            if (any(res$Type[check[1:2]]=="ip")) utility[i] <- utility[i]+2.5;
         }

         utility[i] <- utility[i]-sqrt(length(check))+sqrt(2);
      }


      ord <- order(utility, decreasing=TRUE);
      dups <- dups[ord];
      for (i in 1:k)
      {
         dups[[i]] <- res$IdDocument[dups[[i]]];
      }


      return(dups);
   }
   ## ----------------------------------------------------------------------------------


   if (!is.numeric(aggressiveness) || length(aggressiveness)!=1 || aggressiveness<0)
      stop("incorrect 'aggressiveness'.");

   if (!is.null(ignoreTitles.like) && !is.character(ignoreTitles.like))
      stop("incorrect 'ignoreTitles.like'.");


   dups <- .lbsFindDuplicateTitles_getDups(conn, surveyDescription, ignoreTitles.like, aggressiveness);
   n <- length(dups);

   if(n == 0) return(NULL);

   removed <- integer(0);
   for (i in 1:n)
   {
# @TODO: Check whether .gtk2.selectDocuments can be moved as an internal function
      ret <- .gtk2.selectDocuments(conn, dups[[i]], sprintf("Select documents to remove (stage %g of %g)", i, n), remove=TRUE);

      if (is.null(ret))
         return(removed);

      removed <- c(removed, ret);
   }

   return(removed);
}

