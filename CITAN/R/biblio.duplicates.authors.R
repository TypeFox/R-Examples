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


#' Indicates, by finding similarities between authors' names,
#' groups of authors that possibly should be merged.
#'
#' The function uses a heuristic **EXPERIMENTAL** algorithm. Its behavior
#' is controlled by the \code{aggressiveness} parameter.
#'
#' Search results are presented in a convenient-to-use graphical dialog box.
#' Note that the calculation often takes a few minutes!
#'
#' The \code{names.like} parameter determines search patterns in an SQL \code{LIKE} format,
#' i.e. an underscore \code{_} matches a single character and a percent sign
#' \code{\%} matches any set of characters. The search is case-insensitive.
#'
#' @title Find groups of authors to be merged (**EXPERIMENTAL**)
#' @param conn connection object, see \code{\link{lbsConnect}}.
#' @param names.like character vector of SQL-LIKE patterns that allow for restricting
#' the search procedure to only given authors' names.
#' @param ignoreWords character vector; words to be ignored.
#' @param minWordLength numeric; minimal word length to be considered.
#' @param orderResultsBy determines results' presentation order; one of
#'    \code{citations}, \code{ndocuments} \code{name}.
#' @param aggressiveness nonnegative integer; controls the search depth.
#' @return
#' List of authors' identifiers to be merged.
#' The first element of each vector is the one marked by the user as \emph{Parent},
#' and the rest are the \emph{Children}.
#' @export
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' listauth <- lbsFindDuplicateAuthors(conn,
#'    ignoreWords=c("van", "von", "der", "no", "author", "name", "available"),
#'    minWordLength=4,
#'    orderResultsBy=c("citations"),
#'    aggressiveness=1);
#' lbsMergeAuthors(conn, listauth);
#' dbCommit(conn);
#' ## ...}
#' @seealso
#' \code{\link{lbsMergeAuthors}},
#' \code{\link{lbsFindDuplicateTitles}},
#' \code{\link{lbsGetInfoAuthors}}
lbsFindDuplicateAuthors <- function(conn,
   names.like=NULL,
   ignoreWords=c("van", "von", "der", "no", "author", "name", "available"),
   minWordLength=4,
   orderResultsBy=c("citations", "ndocuments", "name"),
   aggressiveness=0
)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection


   ## -------- Auxiliary function ---------------------------------------------
   .gtk2.selectAuthors <- function(conn, idAuthors, title="Select authors to merge")
   {
      stopifnot(length(idAuthors) > 1 && is.numeric(idAuthors));
      n <- length(idAuthors);

      if (n>2) {
         window <- .gtk2.progressBar(0, n, info="Preparing data...");
      } else window <- NULL;

      info <- list();
      length(info) <- n;
      for (i in 1:n)
      {
         info[[i]] <- as.list(dbGetQuery(conn, sprintf(
            "SELECT Biblio_Authors.IdAuthor AS IdAuthor, Name, AuthorGroup,
               COUNT(IdDocument) As Documents,
               SUM(Citations) As Citations
            FROM Biblio_Authors
            JOIN (
               SELECT Biblio_AuthorsDocuments.IdAuthor, Biblio_AuthorsDocuments.IdDocument, Biblio_Documents.Citations
               FROM Biblio_AuthorsDocuments
               JOIN Biblio_Documents ON Biblio_AuthorsDocuments.IdDocument=Biblio_Documents.IdDocument
            ) AS DocInfo ON Biblio_Authors.IdAuthor=DocInfo.IdAuthor
            WHERE Biblio_Authors.IdAuthor=%g
            GROUP BY Biblio_Authors.IdAuthor",
         idAuthors[i])));

         if (!is.null(window)) .gtk2.progressBar(i, n, window=window);
      }

      info <- info[order(sapply(info, function(x) x$Name))];


      dialog <- gtkDialogNewWithButtons(NULL, NULL, 0,
         "Merge selected", 1,
         "Do nothing", 0,
         "gtk-cancel", GtkResponseType["reject"], show=FALSE);
      dialog$setDefaultSize(700, 400);
      dialog$setTitle(title);


      data <- data.frame(
         "Parent"     =FALSE,
         "Child"      =FALSE,
         "Id"         =sapply(info, function(x) x$IdAuthor),
         "Name"       =sapply(info, function(x) x$Name),
         "Documents"  =sapply(info, function(x) x$Documents),
         "Citations"  =sapply(info, function(x) x$Citations),
         "AuthorGroup" = sapply(info, function(x) if (is.null(x$AuthorGroup)) "NULL" else x$AuthorGroup)
      );

      data$Id <- as.numeric(data$Id);
      model <- rGtkDataFrame(data);
      tree_view <- gtkTreeView(model);

      sapply(1:ncol(model), function(j) {
         if (j == 1)
         {
            renderer <- gtkCellRendererToggle();
            renderer[["radio"]] <- TRUE;
            column <- gtkTreeViewColumn(colnames(model)[j], renderer, active = j-1)
            tree_view$appendColumn(column)
            gSignalConnect(renderer, "toggled", function(widget, path)
            {
               if (model[as.integer(path)+1][1]$Parent == FALSE)
               {
                  for (o in 1:n)
                     model[o][1]$Parent <<- FALSE;
               }
               model[as.integer(path)+1][1]$Parent <<- TRUE;
               model[as.integer(path)+1][2]$Child <<- FALSE;
            })
         } else if (j == 2)
         {
            renderer <- gtkCellRendererToggle();
            renderer[["radio"]] <- FALSE;
            column <- gtkTreeViewColumn(colnames(model)[j], renderer, active = j-1)
            tree_view$appendColumn(column)
            gSignalConnect(renderer, "toggled", function(widget, path)
            {
               if (!model[as.integer(path)+1][1]$Parent)
                  model[as.integer(path)+1][2]$Child <<- !model[as.integer(path)+1][2]$Child;
            })
         } else {
            renderer <- gtkCellRendererText();
            renderer[["wrap-width"]] <- 350;
            renderer[["wrap-mode"]] <- GtkWrapMode["word"];
            column <- gtkTreeViewColumn(colnames(model)[j], renderer, text = j-1)
            tree_view$appendColumn(column)
         }
      })

      renderer <- gtkCellRendererToggle();
      renderer[["radio"]] <- TRUE;
      column <- gtkTreeViewColumn("List doc.", renderer)
      tree_view$insertColumn(column,4)
      gSignalConnect(renderer, "toggled", function(widget, path)
      {
         .gtk2.selectDocuments(conn,
            lbsSearchDocuments(conn, idAuthors=model[as.integer(path)+1][3]$Id),
            parent=dialog, remove=FALSE,
            title=sprintf("Documents by %s (IdAuthor=%g)", model[as.integer(path)+1][4]$Name, model[as.integer(path)+1][3]$Id));
      })

      if (is.null(gtkCheckVersion(2, 10, 0))) tree_view$setGridLines("both");



      swin <- gtkScrolledWindow()
      swin$add(tree_view);
      dialog[["vbox"]]$add(swin);

      doneVal <- dialog$run();

      dialog$destroy();

      data <- as.data.frame(model);


      if (doneVal==1) {
         x <- c(data$Id[data$Parent], data$Id[data$Child]);
         if (length(x) < 2)
         {
            warning("you should select one parent and at least one child.");
            return(integer(0));
         }
         return(x);
      } else if (doneVal<0) {
         return(NULL);
      } else return(integer(0));
   }
   ## --------------------------------------------------------------------------


   ## -------- Auxiliary function ----------------------------------------------
   .lbsFindDuplicateAuthors_split2Words <- function(what, ignoreWords, minWordLength)
   {
      n <- length(what);

      window <- .gtk2.progressBar(0, n, info="Preprocessing data...");

      mtch <- gregexpr("[^[:space:]]+", what);
      words <- list();
      length(words) <- n;

      for (i in 1:n)
      {
         match_i <- mtch[[i]];
         m <- length(match_i);
         words[[i]] <- character(0);
         if (m > 0)
         {
            for (j in 1:m)
            {
               wrd <- substr(what[i], match_i[j], match_i[j]+attr(match_i,"match.length")[j]-1);

               if (any(wrd == ignoreWords)) next;
               if (substr(wrd, nchar(wrd), nchar(wrd)) == ".") next; # last letter == "."
               if (nchar(wrd) < minWordLength) next;

               words[[i]] <- c(words[[i]], wrd);
            }
         }
         .gtk2.progressBar(i, n, window=window);
      }

      return(words);
   }
   ## --------------------------------------------------------------------------


   ## -------- Auxiliary function ----------------------------------------------
   .lbsFindDuplicateAuthors_getDupsGraph <- function(conn, names.like, ignoreWords, minWordLength, aggressiveness, orderResultsBy)
   {
      query <- sprintf(
         "SELECT Biblio_Authors.IdAuthor AS IdAuthor, Name,
            COUNT(IdDocument) As Documents,
            SUM(Citations) As Citations
         FROM Biblio_Authors
         JOIN (
            SELECT Biblio_AuthorsDocuments.IdAuthor, Biblio_AuthorsDocuments.IdDocument, Biblio_Documents.Citations
            FROM Biblio_AuthorsDocuments
            JOIN Biblio_Documents ON Biblio_AuthorsDocuments.IdDocument=Biblio_Documents.IdDocument
         ) AS DocInfo ON Biblio_Authors.IdAuthor=DocInfo.IdAuthor
         WHERE %s
         GROUP BY Biblio_Authors.IdAuthor
         ORDER BY %s DESC",
         ifelse(is.null(names.like), "1", sprintf("Name LIKE '%s'", names.like)),
         sqlSwitchOrNULL(orderResultsBy, c("citations", "ndocuments", "name"), c("Citations", "Documents", "Name"))
      );

      authors <- dbGetQuery(conn, query);
      n <- nrow(authors);

      if (n <= 1) return(NULL);

      authors$Group <- NA;
      authors$Name  <- toupper(authors$Name);
      authors$Row   <- 1:n;


      words <- .lbsFindDuplicateAuthors_split2Words(authors$Name, ignoreWords, minWordLength);

      window <- .gtk2.progressBar(0, n, info=sprintf("Generating dependency graph for %g author names... ",n));

      rel <- list();
      length(rel) <- n;
      for (i in 1:n)
      {
         rel[[i]] <- i;

         words_i <- words[[i]];
         nw <- length(words_i);

         if (nw > 0)
         {
            for (j in 1:nw)
            {
               if (aggressiveness == 0)
               {
                  idx <- grep(words_i[j], authors$Name, fixed=TRUE);
               } else
               {
                  idx <- agrep(words_i[j], authors$Name);
                  if (aggressiveness == 1)
                  {
                     if (length(idx) > 10)
                     {
                        idx2 <- grep(words_i[j], authors$Name, fixed=TRUE);
                        if (length(idx)/length(idx2)>1.5)
                           idx <- idx2;
                     }
                  }
               }

               if (length(idx) > 0)
               {
                  rel[[i]] <- c(rel[[i]],idx);
               }
            }
            rel[[i]] <- unique(rel[[i]]);
         }

         .gtk2.progressBar(i,n,window=window);
      }

      names(rel) <- authors$IdAuthor;
      return(rel);
   }
   ## --------------------------------------------------------------------------


   ## -------- Auxiliary function ----------------------------------------------
   .lbsFindDuplicateAuthors_getSimClusters <- function(x, y)
   {
      xy <- c(x, y);
      nxcapb <- length(x)+length(y)-length(unique(xy));

      return(max(nxcapb/length(x),nxcapb/length(y)));
   }
   ## --------------------------------------------------------------------------


   ## -------- Auxiliary function ----------------------------------------------
   .lbsFindDuplicateAuthors_getDupsFromGraph <- function(graph)
   {
      dups <- list();
      k <- 1;
      n <- length(graph);

      window <- .gtk2.progressBar(0, n, info=sprintf("Trying to group %g author names... ",n));
      selected <- rep(FALSE, n);

      for (i in 1:n)
      {
         graph_i <- graph[[i]];
         n_i <- length(graph_i);

         if ((n_i > 1) && (!selected[i]))
         {
            stopifnot(graph_i[1] == i);
            dups[[k]] <- graph_i;

            for (j in 2:n_i)
            {
               if (length(unique(c(dups[[k]], graph[[ graph_i[j] ]]))) > 15) next;
               # try to merge the clusters
               if (.lbsFindDuplicateAuthors_getSimClusters(dups[[k]], graph[[ graph_i[j] ]]) > 0.6)
               {
                  dups[[k]] <- unique(c(dups[[k]], graph[[ graph_i[j] ]]));
                  selected[j] <- TRUE;
               }
            }
            k <- k+1;
         }
         selected[i] <- TRUE;



         .gtk2.progressBar(i, n, window=window);
      }

      x <- sapply(dups, function(x) length(x));

      if (k-1 == 0) return(list());

      for (i in 1:(k-1))
      {
         dups[[i]] <- as.numeric(names(graph)[dups[[i]]]);
      }

      return(dups);
   }
   ## --------------------------------------------------------------------------




   orderResultsBy <- match.arg(orderResultsBy);

   if (!is.numeric(aggressiveness) || length(aggressiveness)!=1 || aggressiveness<0)
      stop("incorrect 'aggressiveness'.");

   if (!is.null(ignoreWords) && !is.character(ignoreWords))
      stop("incorrect 'ignoreWords'.");
   ignoreWords <- toupper(ignoreWords);

   if (!is.null(names.like) && (!is.character(names.like) || length(names.like)!=1))
      stop("incorrect 'names.like'.");



   graph <- .lbsFindDuplicateAuthors_getDupsGraph(conn, names.like, ignoreWords, minWordLength, aggressiveness, orderResultsBy);
   if (is.null(graph)) return(NULL);


   dups <- .lbsFindDuplicateAuthors_getDupsFromGraph(graph);



   n <- length(dups);

   if(n == 0) return(NULL);

   merged <- list();
   k <- 0;
   for (i in 1:n)
   {
      ret <- .gtk2.selectAuthors(conn, dups[[i]], sprintf("Select authors to merge (stage %g of %g)", i, n));

      if (is.null(ret))
         return(merged);

      if (length(ret) > 0)
      {
         k <- k+1;
         merged[[k]] <- ret;
      }
   }

   return(merged);
}




#' Deletes given documents from a Local Bibliometric Storage.
#'
#' For safety reasons, an SQL transaction  opened at the beginning of the
#' removal process is not committed (closed) automatically.
#' You should do it on your own (or rollback it), see Examples below.
#'
#' @title Delete given documents
#' @param conn a connection object as produced by \code{\link{lbsConnect}}.
#' @param idDocuments a list of numeric vectors or a numeric vector;
#' document identifiers (see \code{IdDocument} in the table \code{Biblio_Documents})
#' to be deleted.
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' listdoc <- lbsFindDuplicateTitles(conn,
#'    ignoreTitles.like=c("In this issue\%", "\%Editorial", "\%Introduction",
#'    "\%In this issue", "Letter to \%", "\%Preface"),
#'    aggressiveness=2);
#' lbsDeleteDocuments(conn, listdoc);
#' dbCommit(conn);
#' ## ...}
#' @export
#' @return
#' \code{TRUE} on success.
#' @seealso
#' \code{\link{lbsGetInfoDocuments}},
#' \code{\link{lbsFindDuplicateTitles}}
lbsDeleteDocuments <- function(conn, idDocuments)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection

   idDocuments <- unlist(idDocuments);
   if (!is.numeric(idDocuments) || length(idDocuments) == 0)
      stop("incorrect 'idDocuments'.");

   dbBegin(conn);

   for (i in 1:length(idDocuments))
   {
      dbExecQuery(conn, sprintf("DELETE FROM Biblio_DocumentsSurveys WHERE IdDocument=%g", idDocuments[i]));
      dbExecQuery(conn, sprintf("DELETE FROM Biblio_AuthorsDocuments WHERE IdDocument=%g", idDocuments[i]));
      dbExecQuery(conn, sprintf("DELETE FROM Biblio_Documents WHERE IdDocument=%g", idDocuments[i]));
   }

   dbExecQuery(conn,
      "DELETE FROM Biblio_Authors
      WHERE IdAuthor IN (
         SELECT DISTINCT Biblio_Authors.IdAuthor FROM Biblio_Authors
         LEFT JOIN Biblio_AuthorsDocuments ON Biblio_Authors.IdAuthor=Biblio_AuthorsDocuments.IdAuthor
         WHERE Biblio_AuthorsDocuments.IdDocument IS NULL
      )");

   warning("Transaction has not been committed yet. Do-it-yourself with dbCommit(...).");

   return(TRUE);
}



#' Merges given sets of authors. For each group, the function
#' maps all the related documents to a
#' distinguished \emph{parent} author (the first in a list) and
#' removes the other, unused from then on, records (\emph{children}).
#'
#' This function is useful when one author is represented by many
#' records in a Local Bibliometric Storage (a typical situation in case of
#' data gathered from on-line bibliographic databases),
#' e.g. prof. John Thomas Smith
#' appears as 'Smith J.' and 'Smith J.T.'. Some merge procedures
#' are often absolutely necessary if we would like to assess the
#' impact of authors reliably.
#'
#' Note that you may use \code{\link{lbsFindDuplicateAuthors}}
#' to generate input to this function. It will try to suggest which
#' records should be merged (see Examples below).
#'
#' For safety reasons, an SQL transaction  opened at the beginning of the
#' removal process is not committed (closed) automatically.
#' You should do it on your own (or rollback it), see Examples below.
#'
#' @title Merge given authors
#' @param conn a connection object as produced by \code{\link{lbsConnect}}.
#' @param idAuthors list of numeric vectors, each consisting of at least 2 authors' identifiers
#' (see \code{IdAuthor} in the table \code{Biblio_Authors});
#' every first element of a vector becomes a \emph{parent} to which other
#' records are merged.
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' listauth <- lbsFindDuplicateAuthors(conn,
#'    ignoreWords=c("van", "von", "der", "no", "author", "name", "available"),
#'    minWordLength=4,
#'    orderResultsBy=c("citations"),
#'    aggressiveness=1);
#' lbsMergeAuthors(conn, listauth);
#' dbCommit(conn);
#' ## ...}
#' @return
#' \code{TRUE} on success.
#' @export
#' @seealso
#' \code{\link{lbsFindDuplicateAuthors}},
#' \code{\link{lbsGetInfoAuthors}},
#' \code{\link{lbsAssess}}
lbsMergeAuthors <- function(conn, idAuthors)
{
   .lbsCheckConnection(conn); # will stop on invalid/dead connection

   if (!is.list(idAuthors) || length(idAuthors) == 0 || !is.numeric(unlist(idAuthors)))
      stop("incorrect 'idAuthors'.");

   dbBegin(conn);

   for (i in 1:length(idAuthors))
   {
      stopifnot(length(idAuthors[[i]]) >= 2);
      for (j in 2:length(idAuthors[[i]]))
      {
         dbExecQuery(conn, sprintf("UPDATE Biblio_AuthorsDocuments SET IdAuthor=%g WHERE IdAuthor=%g", idAuthors[[i]][1], idAuthors[[i]][j]));
         dbExecQuery(conn, sprintf("DELETE FROM Biblio_Authors WHERE IdAuthor=%g", idAuthors[[i]][j]));
      }
   }

   warning("Transaction has not been committed yet. Do-it-yourself with dbCommit(...).");

   return(TRUE);
}
