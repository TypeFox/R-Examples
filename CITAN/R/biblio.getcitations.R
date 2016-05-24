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


#' Creates ordered citation sequences of authors in a Local Bibliometric Storage.
#'
#' A citation sequence is a numeric vector
#' consisting of citation counts of all the documents mapped to
#' selected  authors.
#' However, the function may take into account only the documents
#' from a given Survey (using \code{surveyDescription}
#' parameter) or of chosen types (\code{documentTypes}).
#'
#' @title Fetch authors' citation sequences
#' @param conn a connection object as produced by \code{\link{lbsConnect}}.
#' @param documentTypes character vector or \code{NULL}; specifies document types to restrict to;
#'    a combination of \code{Article}, \code{Article in Press}, \code{Book}, \code{Conference Paper},
#'    \code{Editorial}, \code{Erratum}, \code{Letter}, \code{Note}, \code{Report}, \code{Review},
#'    \code{Short Survey}. \code{NULL} means no restriction.
#' @param surveyDescription single character string or \code{NULL}; survey to restrict to or \code{NULL} for no restriction.
#' @param idAuthors numeric vector of authors' identifiers for which the sequences are to be created or \code{NULL} for all authors in the database.
#' @param verbose logical; \code{TRUE} to inform about the progress of the process.
#' @return A list of non-increasingly ordered numeric vectors is returned. Each element of the list corresponds
#' to a citation sequence of some author. List \code{names} attribute are
#' set to authors' names. Moreover, each vector has a set \code{IdAuthor}
#' attribute, which uniquely identifies the corresponding record in the table \code{Biblio_Authors}.
#' Citation counts come together with \code{IdDocument}s (vector elements are named).
#'
#' The list of citation sequences may then be used to calculate
#' authors' impact using \code{\link{lbsAssess}} (see Examples below).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' conn <- lbsConnect("Bibliometrics.db");
#' ## ...
#' citseq <- lbsGetCitations(conn,
#' 	surveyDescription="Scientometrics", documentTypes="Article",
#' 	idAuthors=c(39264,39265,39266));
#' print(citseq);
#' ## $`Liu X.`                                # Author name
#' ## 40116 34128 39122 29672 32343 32775      # IdDocument
#' ##    11     4     1     0     0     0      # Citation count
#' ## attr(,"IdAuthor")
#' ## [1] 39264                                # IdAuthor
#' ##
#' ## $`Xu Y.`
#' ## 38680 38605 40035 40030 40124 39829 39745 29672
#' ##    30    14     8     6     6     5     3     0
#' ## attr(,"IdAuthor")
#' ## [1] 39265
#' ##
#' ## $`Wang Y.`
#' ## 29992 29672 29777 32906 33858 33864 34704
#' ##     1     0     0     0     0     0     0
#' ## attr(,"IdAuthor")
#' ## [1] 39266
#' print(lbsAssess(citseq,
#'    f=list(length, sum, index.h, index.g, function(x) index.rp(x,1),
#'        function(x) sqrt(prod(index.lp(x,1))),
#'        function(x) sqrt(prod(index.lp(x,Inf)))),
#'    captions=c("length", "sum", "index.h", "index.g", "index.w",
#'    "index.lp1", "index.lpInf")));
#' ##      Name length sum index.h index.g index.w index.lp1 index.lpInf
#' ## 3   Xu Y.      8  72       5       8       7  8.573214    5.477226
#' ## 2 Wang Y.      7   1       1       1       1  1.000000    1.000000
#' ## 1  Liu X.      6  16       2       4       3  4.157609    3.316625
#' ## ...
#' dbDisconnect(conn);}
#' @seealso
#' \code{\link{lbsConnect}},
#' \code{\link{lbsAssess}}
lbsGetCitations <- function(conn,
	documentTypes=NULL,
	surveyDescription=NULL,
	idAuthors=NULL,
	verbose=TRUE
)
{
	.lbsCheckConnection(conn); # will stop on invalid/dead connection

	# -----------------------------------------------------
	# Data set restrictions & subset stats

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


	if (verbose)
	{
		cat("Data set restrictions:\n");
		cat(sprintf("\tSurvey:         %s.\n", ifelse(is.null(surveyDescription), "<ALL>", surveyDescription)));
		cat(sprintf("\tDocument types: %s.\n", ifelse(is.null(documentTypesShort), "<ALL>", paste(documentTypesShort, collapse=", "))));
		cat("\n");
	}

	# ---------------------------------------------------------------------


	if (!is.null(idAuthors) && (!is.numeric(idAuthors) || any(!is.finite(idAuthors))))
		stop("incorrect 'idAuthors' given");

	if (length(idAuthors) == 0)
	{
		query <- sprintf("
		SELECT IdAuthor
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

		idAuthors <- dbGetQuery(conn, query)[,1];

		if (length(idAuthors) == 0) return(list());
	}




	i <- 1L;
	k <- 0L;
	n <- as.integer(length(idAuthors));
	citseq <- list();
	length(citseq) <- n;

	if (verbose)
	{
		cat("Creating citation sequences... ");
		window <- .gtk2.progressBar(0, n, info=sprintf("Creating %g citation sequences... ",n));
	}

	while (i <= n)
	{
		query <- sprintf("
			SELECT DISTINCT Biblio_Authors2.Name AS Name, Biblio_Documents.IdDocument AS IdDocument, Biblio_Documents.Citations AS Citations
			FROM
			(
				SELECT Name, IdAuthor FROM Biblio_Authors WHERE IdAuthor=%s
			) AS Biblio_Authors2
			JOIN Biblio_AuthorsDocuments ON (Biblio_Authors2.IdAuthor=Biblio_AuthorsDocuments.IdAuthor)
			JOIN ViewBiblio_DocumentsSurveys ON (Biblio_AuthorsDocuments.IdDocument=ViewBiblio_DocumentsSurveys.IdDocument)
			JOIN Biblio_Documents ON (Biblio_AuthorsDocuments.IdDocument=Biblio_Documents.IdDocument)
			WHERE %s
			ORDER BY Biblio_Documents.Citations DESC",
			sqlNumericOrNULL(idAuthors[i]),
			subQueryWhere
		);

		AuthorInfo <- dbGetQuery(conn, query);

		if (nrow(AuthorInfo) > 0)
		{
			names(citseq)[i] <- AuthorInfo$Name[1];
			citseq[[i]] <- as.numeric(AuthorInfo$Citations);
			names(citseq[[i]]) <- as.numeric(AuthorInfo$IdDocument);
			attr(citseq[[i]], "IdAuthor") <- idAuthors[i];
			k <- k+1L;
		}

		if (verbose) .gtk2.progressBar(i,n,window=window);

		i <- i+1L;
	}


	if (verbose) cat(sprintf("OK, %g of %g records read.\n", k, n));

	return(citseq);
}
