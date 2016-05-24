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

#' @title Connect to a Local Bibliometric Storage
#'
#' @description
#' Connects to a Local Bibliometric Storage handled by the SQLite engine
#' (see \pkg{RSQLite} package documentation).
#'
#' @details
#' Do not forget to close the connection (represented by the connection object returned)
#' with the \code{\link{lbsDisconnect}} function after use.
#'
#' Please note that the database may be also accessed by using
#' lower-level functions from the \pkg{DBI} package called on the
#' returned connection object. The table-view structure of a Local
#' Bibliometric Storage is presented in the man page of the
#' \code{\link{lbsCreate}} function.
#'
#' @export
#' @param dbfilename filename of an SQLite database.
#' @return An object of type \code{SQLiteConnection}, used to communicate with the SQLite engine.
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db")
#' ## ...
#' lbsDisconnect(conn)}
#' @seealso \code{\link{lbsCreate}},
#' \code{\link{lbsDisconnect}}
lbsConnect <- function(dbfilename)
{
   if (length(dbfilename)!=1 || !is.character(dbfilename))
      stop("incorrect 'dbfilename' given");

   drv <- dbDriver("SQLite");
   conn <- dbConnect(drv, dbname = dbfilename);

   objects <- dbListTables(conn);
   tables <- objects[substr(objects,1,7) == "Biblio_"];
   views  <- objects[substr(objects,1,11) == "ViewBiblio_"];


   if (length(tables) == 0 && length(views) == 0)
   {
      warning("Your Local Bibliometric Storage is empty. Use lbsCreate(...) to establish one.");
   } else if (any((tables == "Biblio_Countries") | (tables == "Biblio_Languages")))
   {
      warning("Your Local Bibliometric Storage seems to be created with an older version of CITAN. Please re-create the database.");
   }

   return(conn);
}


#' Disconnects from a Local Bibliometric Storage.
#'
#'
#' @title Disconnect from a Local Bibliometric Storage
#' @param conn database connection object, see \code{\link{lbsConnect}}.
#' @export
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' lbsDisconnect(conn);}
#' @seealso \code{\link{lbsConnect}}
lbsDisconnect <- function(conn)
{
   .lbsCheckConnection(conn); # will stop on an invalid/dead connection

   dbDisconnect(conn);
}
