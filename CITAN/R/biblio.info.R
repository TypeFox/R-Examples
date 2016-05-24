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


#' Retrieves basic information on given authors.
#'
#' @title Retrieve author information
#' @param conn a connection object as produced by \code{\link{lbsConnect}}.
#' @param idAuthors a numeric or integer vector with author identifiers (see column \code{IdAuthor} in the table \code{Biblio_Authors}).
#' @return
#' A list of \code{authorinfo} objects, that is lists with the following components:
#' \itemize{
#' \item \code{IdAuthor} --- numeric; author's identifier in the table \code{Biblio_Authors},
#' \item \code{Name} --- character; author's name.
#' \item \code{AuthorGroup} --- character; author group (used to merge author records).
#' }
#' @examples
#' \dontrun{
#' conn <- dbBiblioConnect("Bibliometrics.db");
#' ## ...
#' id <- lbsSearchAuthors(conn, c("Smith\%", "Knuth D.E.", "V_n \%"));
#' lbsGetInfoAuthors(conn, id);
#' ## ...}
#' @export
#' @seealso \code{\link{lbsSearchAuthors}}, \code{\link{lbsSearchDocuments}},
#' \code{\link{lbsGetInfoDocuments}},\cr
#' \code{\link{as.character.authorinfo}}, \code{\link{print.authorinfo}},
lbsGetInfoAuthors <- function(conn, idAuthors)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   if (is.null(idAuthors) || (class(idAuthors) != "numeric" && class(idAuthors) != "integer"))
      stop("'idAuthors' must be a nonempty numeric vector.");

   query <- sprintf("SELECT IdAuthor, Name, AuthorGroup FROM Biblio_Authors WHERE IdAuthor IN (%s)",
         paste(idAuthors, collapse=","));

   res <- dbGetQuery(conn, query);
   n <- nrow(res);
   out <- list();
   length(out) <- n;



   if (n < 1) return(NULL);

   for (i in 1:n)
   {
      out[[i]] <- list(IdAuthor=res[i,1], Name=res[i,2], AuthorGroup=res[i,3]);
      class(out[[i]]) <- "authorinfo";
   }

   return(out);
}


#' Retrieves information on given documents.
#'
#' @title Retrieve document information
#' @param conn a connection object as produced by \code{\link{lbsConnect}}.
#' @param idDocuments a numeric or integer vector with document identifiers (see column \code{IdDocument} in the table \code{Biblio_Documents}).
#' @return
#' A list of \code{docinfo} objects, that is lists with the following components:
#' \itemize{
#' \item \code{IdDocument} --- numeric; document identifier in the table \code{Biblio_Documents},
#' \item \code{Authors} --- list of \code{authorinfo} objects (see e.g. \code{\link{as.character.authorinfo}}).
#' \item \code{Title} --- title of the document,
#' \item \code{BibEntry} --- bibliographic entry,
#' \item \code{AlternativeId} --- unique character identifier,
#' \item \code{Pages} --- number of pages,
#' \item \code{Citations} --- number of citations,
#' \item \code{Year} --- publication year,
#' \item \code{Type} --- document type, e.g. \code{Article} or \code{Conference Paper}.
#' }
#' @examples
#' \dontrun{
#' conn <- dbBiblioConnect("Bibliometrics.db");
#' ## ...
#' id <- lbsSearchDocuments(conn,
#' 	idAuthors=lbsSearchAuthors(conn, "Knuth\%"));
#' lbsGetInfoDocuments(conn, id);
#' ## ...}
#' @export
#' @seealso \code{\link{print.docinfo}}, \code{\link{lbsSearchDocuments}},
#' \code{\link{lbsGetInfoAuthors}},\cr
#' \code{\link{as.character.authorinfo}}, \code{\link{as.character.docinfo}}
lbsGetInfoDocuments <- function(conn, idDocuments)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   if (is.null(idDocuments) || (class(idDocuments) != "numeric" && class(idDocuments) != "integer"))
      stop("'idDocuments' must be a nonempty numeric vector.");

   idDocuments_str <- paste(idDocuments, collapse=",");
   query <- sprintf("SELECT IdDocument, Title, BibEntry, AlternativeId, Pages, Citations, Year, Type
         FROM Biblio_Documents WHERE IdDocument IN (%s) ORDER BY IdDocument",
         idDocuments_str);

   res <- dbGetQuery(conn, query);
   n <- nrow(res);
   out <- list();
   length(out) <- n;


   if (n < 1) return(NULL);

   query2 <- sprintf("SELECT DISTINCT IdDocument, Biblio_Authors.IdAuthor, Biblio_Authors.Name, Biblio_Authors.AuthorGroup FROM
      Biblio_AuthorsDocuments
      JOIN Biblio_Authors ON Biblio_Authors.IdAuthor=Biblio_AuthorsDocuments.IdAuthor
      WHERE IdDocument IN (%s) ORDER BY IdDocument",
      idDocuments_str);
   res2 <- dbGetQuery(conn, query2);
   stopifnot(nrow(res2) >= 1);


   i <- 1;
   k <- 1;
   m <- nrow(res2);

   while (i <= m)
   {
      stopifnot(res2$IdDocument[i] == res$IdDocument[k]);

      j <- i+1;
      while (j<=m && res2$IdDocument[i] == res2$IdDocument[j])
         j <- j+1;

      authors <- list();
      length(authors) <- j-i;
      for (u in i:(j-1))
      {
         authors[[u-i+1]] <- list(IdAuthor=res2[u,2], Name=res2[u,3], AuthorGroup=res2[u,4]);
         class(authors[[u-i+1]]) <- "authorinfo";
      }

      doc <- list(IdDocument=res[k,1], Authors=authors, Title=res[k,2], BibEntry=res[k,3],
         AlternativeId=res[k,4], Pages=res[k,5], Citations=res[k,6],
         Year=res[k,7], Type=.lbs_DocumentType_ShortToFull(res[k,8]));

      class(doc) <- "docinfo";

      out[[k]] <- doc;


      i <- j;
      k <- k+1;
   }
   stopifnot(k-1 == n);


   return(out);
}


#' Converts an object of class \code{docinfo} to a character string.
#' Such an object is  returned by e.g. \code{\link{lbsGetInfoDocuments}}.
#'
#' A \code{docinfo} object is a list with the following components:
#' \itemize{
#' \item \code{IdDocument} --- numeric; document identifier in the table \code{Biblio_Documents},
#' \item \code{Authors} --- list of \code{authorinfo} objects (see e.g. \code{\link{as.character.authorinfo}}).
#' \item \code{Title} --- title of the document,
#' \item \code{BibEntry} --- bibliographic entry,
#' \item \code{AlternativeId} --- unique character identifier,
#' \item \code{Pages} --- number of pages,
#' \item \code{Citations} --- number of citations,
#' \item \code{Year} --- publication year,
#' \item \code{Type} --- type of document, see \code{\link{lbsCreate}}.
#' }
#'
#' @title Coerce a docinfo object to character string
#' @param x a single object of class \code{docinfo}.
#' @param ... unused.
#' @return A character string
#' @export
#' @method as.character docinfo
#' @seealso \code{\link{lbsSearchDocuments}},
#' \code{\link{as.character.authorinfo}}, \code{\link{print.docinfo}},\cr
#' \code{\link{lbsGetInfoDocuments}}
as.character.docinfo <- function(x, ...)
{
   ret <-            sprintf("IdDocument:    %g", x$IdDocument);
   ret <- paste(ret, sprintf("AlternativeId: %s", x$AlternativeId), sep="\n");
   ret <- paste(ret, sprintf("Title:         %s", x$Title), sep="\n");
   ret <- paste(ret, sprintf("BibEntry:      %s", x$BibEntry), sep="\n");
   ret <- paste(ret, sprintf("Year:          %s", x$Year), sep="\n");
   ret <- paste(ret, sprintf("Type:          %s", x$Type), sep="\n");
   ret <- paste(ret, sprintf("Citations:     %s", x$Citations), sep="\n");
   ret <- paste(ret, sprintf("Authors:       %s\n",
      paste(
         sapply(x$Authors, function(y) paste(y$Name, y$IdAuthor, y$AuthorGroup, sep="/")),
         collapse=", ")
      ), sep="\n");
   return(ret);
}

#' Prints out an object of class \code{docinfo}. Such an object is returned by e.g. \code{\link{lbsGetInfoDocuments}}.
#'
#' For more information see man page for \code{\link{as.character.docinfo}}.
#'
#' @title Print a docinfo object
#' @param x an object of class \code{docinfo}.
#' @param ... unused.
#' @export
#' @method print docinfo
#' @seealso \code{\link{as.character.docinfo}}, \code{\link{lbsSearchDocuments}}, \code{\link{lbsGetInfoDocuments}}
print.docinfo <- function(x, ...)
{
   cat(as.character(x));
}




#' Converts an object of class \code{authorinfo} to a character string.
#' Such an object is returned by e.g. \code{\link{lbsGetInfoAuthors}}.
#'
#' An \code{authorinfo} object  is a list with the following components:
#' \itemize{
#' \item \code{IdAuthor} --- numeric; author's identifier in the table \code{Biblio_Authors},
#' \item \code{Name} --- character; author's name.
#' }
#'
#' @title Coerce an authorinfo object to character string
#' @param x a single object of class \code{authorinfo}.
#' @param ... unused.
#' @return A character string
#' @export
#' @method as.character authorinfo
#' @seealso \code{\link{print.authorinfo}}, \code{\link{lbsSearchAuthors}}, \code{\link{lbsGetInfoAuthors}}
as.character.authorinfo <- function(x, ...)
{
   ret <-            sprintf("IdAuthor:    %g",   x$IdAuthor);
   ret <- paste(ret, sprintf("Name:        %s", x$Name), sep="\n");
   ret <- paste(ret, sprintf("AuthorGroup: %s\n", x$AuthorGroup), sep="\n");
   return(ret);
}

#' Prints out an object of class \code{authorinfo}. Such an object is returned by e.g. \code{\link{lbsGetInfoAuthors}}.
#'
#' For more information see man page for \code{\link{as.character.authorinfo}}.
#'
#' @title Print an authorinfo object
#' @param x an object of class \code{authorinfo}.
#' @param ... unused.
#' @method print authorinfo
#' @export
#' @seealso \code{\link{as.character.authorinfo}}, \code{\link{lbsSearchAuthors}}, \code{\link{lbsGetInfoAuthors}}
print.authorinfo <- function(x, ...)
{
   cat(as.character(x));
}



