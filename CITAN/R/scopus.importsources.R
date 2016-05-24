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


#' Imports \emph{SciVerse Scopus} covered titles and their ASJC codes to an empty Local Bibliometric Storage (\acronym{LBS}).
#'
#' This function should be called prior to importing any document information
#' to the LBS with the function \code{\link{lbsImportDocuments}}.
#'
#' Note that adding all the sources takes some time.
#'
#' Only elementary ASJC and \emph{SciVerse Scopus} source data
#' read from \code{\link{Scopus_ASJC}} and \code{\link{Scopus_SourceList}}
#' will be added to the LBS (\code{Biblio_Categories}, \code{Biblio_Sources}, \code{Biblio_SourcesCategories}).
#'
#'
#' @title Import SciVerse Scopus coverage information and ASJC codes to a Local Bibliometric Storage
#' @param conn a connection object, see \code{\link{lbsConnect}}.
#' @param verbose logical; \code{TRUE} to display progress information.
#' @return \code{TRUE} on success.
#' @export
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' lbsCreate(conn);
#' Scopus_ImportSources(conn);
#' ## ...
#' lbsDisconnect(conn);}
#' @seealso \code{\link{Scopus_ASJC}}, \code{\link{Scopus_SourceList}}, \code{\link{Scopus_ReadCSV}}, \code{\link{lbsConnect}}, \code{\link{lbsCreate}}
Scopus_ImportSources <- function(conn, verbose=T)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   ## ---- check if 3 tables used here are empty -----------------------------
   res <- dbGetQuery(conn, "SELECT COUNT(*) AS Count FROM Biblio_Categories");
   if (res[1,1] != 0) stop("table 'Biblio_Categories' is not empty.");

   res <- dbGetQuery(conn, "SELECT COUNT(*) AS Count FROM Biblio_Sources");
   if (res[1,1] != 0) stop("table 'Biblio_Sources' is not empty.");

   res <- dbGetQuery(conn, "SELECT COUNT(*) AS Count FROM Biblio_SourcesCategories");
   if (res[1,1] != 0) stop("table 'Biblio_SourcesCategories' is not empty.");
   ## ------------------------------------------------------------------------



   ## ----- auxiliary functions ----------------------------------------------

   # /internal/
   # Imports data from Scopus_ASJC to Biblio_Categories
   .Scopus_ImportSources_Categories <- function(conn, verbose)
   {
      if (verbose) cat("Importing Scopus ASJC codes... ");

      n <- nrow(Scopus_ASJC);


      ## ----- prepare queries ----------------------------------------------
      queries <- sprintf("INSERT INTO Biblio_Categories
                ('IdCategory', 'IdCategoryParent', 'Description')
         VALUES (%g,           %g,                '%s');",
            as.integer(Scopus_ASJC$ASJC),
            as.integer(Scopus_ASJC$ASJC_Parent),
            sqlEscapeTrim(Scopus_ASJC$Description)
         );
      ## --------------------------------------------------------------------

      dbBegin(conn);
      for (i in 1:n) dbExecQuery(conn, queries[i], TRUE);
      dbCommit(conn);

      if (verbose) cat(sprintf("Done, %g records added.\n", n));
   }



   ## ----- auxiliary function ----------------------------------------------

   # /internal/
   .Scopus_ImportSources_Sources <- function(conn, verbose)
   {
      n <- nrow(Scopus_SourceList);

      if (verbose) cat("Importing Scopus source list... ");


      ## ----- prepare queries ------------------------
      queries <- sprintf("INSERT OR FAIL INTO Biblio_Sources(
            'IdSource', 'AlternativeId', 'Title', 'IsActive', 'IsOpenAccess', 'Type',
            'Impact1', 'Impact2', 'Impact3', 'Impact4', 'Impact5', 'Impact6'
         )
         VALUES (%.0f, '%s', '%s', %s, %s, %s, %s, %s, %s, %s, %s, %s);",
         1:n,
         sqlEscapeTrim(as.character(Scopus_SourceList$SourceId)),
         sqlEscapeTrim(Scopus_SourceList$Title),
         as.integer(Scopus_SourceList$Status == "Active"),
         as.integer(Scopus_SourceList$OpenAccess != "Not OA"),
         sqlSwitchOrNULL(Scopus_SourceList$Type,
            .lbs_SourceTypesFull,
            .lbs_SourceTypesShort
         ),
         sqlNumericOrNULL(Scopus_SourceList$SJR_2009),
         sqlNumericOrNULL(Scopus_SourceList$SNIP_2009),
         sqlNumericOrNULL(Scopus_SourceList$SJR_2010),
         sqlNumericOrNULL(Scopus_SourceList$SNIP_2010),
         sqlNumericOrNULL(Scopus_SourceList$SJR_2011),
         sqlNumericOrNULL(Scopus_SourceList$SNIP_2011)
      );

      ## ------ prepare ASJCs --------------------------

      asjcs <- strsplit(Scopus_SourceList$ASJC, "[[:space:]]*;[[:space:]]*");
      asjcs <- lapply(asjcs, function(x)
         { x<- na.omit(as.integer(x));
           x<- x[x>=1000 & x<=9999];
           x;
         });

      stopifnot(length(asjcs) == nrow(Scopus_SourceList));




      if (verbose) window <- .gtk2.progressBar(0, n, info="Importing source list... ");


      ## ----- exec    queries ------------------------
      omitted <- numeric(0);
      k <- 0;
      dbExecQuery(conn, "PRAGMA journal_mode = MEMORY");
      dbBegin(conn);
      for (i in 1:n)
      {
         tryCatch(
         {
            res <- dbSendQuery(conn, queries[i]);
            dbClearResult(res);
   #        id <- as.numeric(dbGetQuery(conn, "SELECT last_insert_rowid()")[1,1]); ## == i

            m <- length(asjcs[[i]]);
            if (m > 0)
            {
               query2 <- sprintf("INSERT INTO Biblio_SourcesCategories('IdSource', 'IdCategory') %s",
                   paste(
                      sprintf("SELECT %.0f AS 'IdSource', %.0f AS 'IdCategory'", i, asjcs[[i]]),
                      collapse=" UNION "));

               dbExecQuery(conn, query2, TRUE);

               k <- k+m;
            } else warning(sprintf("No ASJC @ row=%g.", i));
         },
         error=function(err)
         {
   #        print(i);
   #        print(err);
            omitted <<- c(omitted, i);
         });

         if (verbose) .gtk2.progressBar(i, n, window=window);
      }
      dbCommit(conn);



      ## -----  now check which rows were added and report missing values ------
      if (verbose) cat(sprintf("Done, %g of %g records added; %g ASJC codes processed.\n", n-length(omitted), n, k));

      if (length(omitted) > 0 && verbose)
      {
         cat(sprintf("Note: %g records omitted @ rows=%s.\n",
            length(omitted), paste(omitted,collapse=","))
         );
      }
   }





   ## ---- Import categories and sources -----------------------------------

   .Scopus_ImportSources_Categories(conn, verbose)

   .Scopus_ImportSources_Sources(conn, verbose)


   ## --------- VACUUM -----------------------------------------------------
   dbExecQuery(conn, "VACUUM", FALSE);
   ## ----------------------------------------------------------------------

   return(TRUE);
}
