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


#' Clears a Local Bibliometric Storage by dropping all tables
#' named \code{Biblio_*} and all views named \code{ViewBiblio_*}.
#'
#' For safety reasons, an SQL transaction opened at the beginning of the
#' removal process is not committed (closed) automatically.
#' You should do manually (or rollback it), see Examples below.
#'
#' @title Clear a Local Bibliometric Storage
#' @param conn database connection object, see \code{\link{lbsConnect}}.
#' @param verbose logical; \code{TRUE} to be more verbose.
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' lbsClear(conn);
#' dbCommit(conn);
#' lbsCreate(conn);
#' Scopus_ImportSources(conn);
#' ## ...
#' lbsDisconnect(conn);}
#'
#' @return \code{TRUE} on success.
#' @export
#' @importFrom RSQLite dbBegin
#'
#' @seealso \code{\link{lbsConnect}}, \code{\link{lbsCreate}},
#' \code{\link{Scopus_ImportSources}}, \code{\link{lbsDeleteAllAuthorsDocuments}}
#' \code{\link{dbCommit}}, \code{\link{dbRollback}}
lbsClear <- function(conn, verbose=TRUE)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection

   dbBegin(conn);

   objects <- dbListTables(conn);
   tables <- objects[substr(objects,1,7) == "Biblio_"];
   views  <- objects[substr(objects,1,11) == "ViewBiblio_"];


   if (length(tables) == 0 && length(views) == 0)
   {
      warning("Your Local Bibliometric Storage is already empty.");
      return(TRUE);
   }



   if (length(tables) != 0)
      for (i in 1:length(tables))
      {
         if (verbose) cat(sprintf("Dropping table '%s'... ", tables[i]));
         dbExecQuery(conn, sprintf("DROP TABLE %s;", tables[i]), TRUE);
         if (verbose) cat("DONE.\n");
      }


   if (length(views) != 0)
      for (i in 1:length(views))
      {
         if (verbose) cat(sprintf("Dropping view '%s'... ", views[i]));
         dbExecQuery(conn, sprintf("DROP VIEW %s;", views[i]), TRUE);
         if (verbose) cat("DONE.\n");
      }


   warning("Transaction has not been committed yet. Please use dbCommit(...) or dbRollback(...).");

   return(TRUE);
}





#' Deletes author, citation, document, and survey information from a Local Bibliometric
#' Storage.
#'
#' For safety reasons, an SQL transaction opened at the beginning of the
#' removal process is not committed (closed) automatically.
#' You should do manually (or rollback it), see Examples below.
#'
#' @title Delete all authors, documents and surveys from a Local Bibliometric Storage
#' @param conn database connection object, see \code{\link{lbsConnect}}.
#' @param verbose logical; \code{TRUE} to be more verbose.
#' @return \code{TRUE} on success.
#' @export
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db")
#' lbsDeleteAllAuthorsDocuments(conn)
#' dbCommit(conn)
#' ## ...
#' lbsDisconnect(conn)}
#' @seealso
#' \code{\link{lbsClear}},
#' \code{\link{dbCommit}},
#' \code{\link{dbRollback}}
lbsDeleteAllAuthorsDocuments <- function(conn, verbose=TRUE)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   if (verbose) cat(sprintf("Deleting all author and document information... "));


   dbBegin(conn);

   dbExecQuery(conn, "DELETE FROM Biblio_DocumentsSurveys", TRUE);
   dbExecQuery(conn, "DELETE FROM Biblio_AuthorsDocuments", TRUE);
   dbExecQuery(conn, "DELETE FROM Biblio_Surveys", TRUE);
   dbExecQuery(conn, "DELETE FROM Biblio_Authors", TRUE);
   dbExecQuery(conn, "DELETE FROM Biblio_Citations", TRUE);
   dbExecQuery(conn, "DELETE FROM Biblio_Documents", TRUE);

   if (verbose) cat(sprintf("DONE.\n"));

   warning("Transaction has not been committed yet. Please use dbCommit(...) or dbRollback(...).");

   return(TRUE);
}
