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



#' Cleans up a Local Bibliometric Storage
#' by removing all authors with no documents, fixing documents
#' with missing survey information, and executing the \code{VACUUM}
#' SQL command.
#'
#' @title Clean up a Local Bibliometric Storage
#' @param conn database connection object, see \code{\link{lbsConnect}}.
#' @param newSuveyDescription character; default survey description for documents with missing survey info.
#' @param newSuveyFilename character; default survey filename for documents with missing survey info.
#' @export
#' @return \code{TRUE} on success.
#' @seealso \code{\link{lbsConnect}}, \code{\link{lbsCreate}},
#' \code{\link{Scopus_ImportSources}},
#' \code{\link{lbsDeleteAllAuthorsDocuments}},
#' \code{\link{dbCommit}}, \code{\link{dbRollback}}
lbsTidy <- function(conn, newSuveyDescription="lbsTidy_Merged", newSuveyFilename="lbsTidy_Merged")
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection

   ## REMOVE ALL AUTHORS WITHOUT DOCUMENTS
   dbExecQuery(conn, "DELETE FROM Biblio_Authors WHERE IdAuthor IN (
      SELECT Biblio_Authors.IdAuthor
      FROM Biblio_Authors
      LEFT JOIN Biblio_AuthorsDocuments ON Biblio_Authors.IdAuthor=Biblio_AuthorsDocuments.IdAuthor
      WHERE IdDocument IS NULL
   )");
   chg <- dbGetQuery(conn, "SELECT changes()")[1,1];
   cat(sprintf("Deleted %g authors with no documents.\n", chg));



   ## FIX ALL DOCUMENTS WITH NO SURVEY INFO
   idDocNoSurvey <- dbGetQuery(conn, "SELECT Biblio_Documents.IdDocument FROM Biblio_Documents LEFT JOIN Biblio_DocumentsSurveys ON Biblio_Documents.IdDocument=Biblio_DocumentsSurveys.IdDocument WHERE Biblio_DocumentsSurveys.IdSurvey IS NULL");
   n <- nrow(idDocNoSurvey);

   if (n > 0)
   {
      query <- sprintf("INSERT INTO Biblio_Surveys('Description', 'FileName', 'Timestamp')
         VALUES(%s, %s, %s)",
         sqlStringOrNULL(newSuveyDescription),
         sqlStringOrNULL(newSuveyFilename),
         sqlStringOrNULL(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      );
      dbExecQuery(conn, query, TRUE);
      idSurvey <- (as.numeric(dbGetQuery(conn, "SELECT last_insert_rowid()")[1,1]));

#       dbBegin(conn);
      for (i in 1:n)
      {
         dbExecQuery(conn, sprintf("INSERT INTO Biblio_DocumentsSurveys('IdDocument', 'IdSurvey') VALUES (%s, %s);",
            sqlNumericOrNULL(idDocNoSurvey[i,1]),
            sqlNumericOrNULL(idSurvey)
         ));
      }
      dbCommit(conn);
   }

   cat(sprintf("Fixed %g documents with no survey information.\n", n));

   ## VACUUM
   dbExecQuery(conn, "VACUUM", FALSE);

   return(TRUE);
}

