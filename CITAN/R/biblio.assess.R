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


#' Given a list of authors' citation sequences, the function calculates
#' values of many impact functions at a time.
#'
#' @title Calculate impact of given authors
#' @param citseq list of numeric vectors, e.g. the output of \code{\link{lbsGetCitations}}.
#' @param f a list of \eqn{n} functions which compute the impact of an author.
#'        The functions must calculate their values using numeric
#'        vectors passed as their first arguments.
#' @param captions a list of \eqn{n} descriptive captions for the functions in \code{f}.
#' @param orderByColumn column to sort the results on. \code{1} for author
#'        names, \code{2} for the first function in \code{f}, \code{3}
#'        for the second, and so on.
#' @param bestRanks if not \code{NULL}, only a given number of authors
#'        with the greatest impact (for each function in \code{f}) will be included in the output.
#' @param verbose logical; \code{TRUE} to inform about the progress of the process.
#' @return A data frame in which each row corresponds to the assessment
#' results of some citation sequence.
#' The first column stands for the authors' names (taken from \code{names(citseq)},
#' the second for the valuation of \code{f[[1]]}, the third for \code{f[[2]]}, and so on.
#' See Examples below.
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
#' library("agop")
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
#' @export
#'
#' @importFrom agop index_h
#'
#' @seealso \code{\link{lbsConnect}}, \code{\link{lbsGetCitations}}
lbsAssess <- function(citseq, f=list(length, index_h), captions=c("length", "index_h"), orderByColumn=2,
	bestRanks=20, verbose=T)
{
	if (!class(citseq)=="list") stop("incorrect 'citseq'");
	if (!is.null(bestRanks) && (!is.numeric(bestRanks) || bestRanks <= 0))
		stop("incorrect 'bestRanks'");

	if (class(f) != "list") stop("incorrect 'f'");
	if (class(captions) != "character") stop("incorrect 'captions'");
	if (length(f) != length(captions)) stop("'f' and 'captions' must be of the same length");

	result <- data.frame(Name=names(citseq));


	for (i in 1:length(f))
	{
		if (verbose) cat(sprintf("Calculating %s... ", captions[i]));
		result <- cbind(result, sapply(citseq,f[[i]]));
		if (verbose) cat("DONE.\n");
	}
	names(result)[2:(length(f)+1)] <- captions;
	result <- result[order(result[,orderByColumn]),]
	rownames(result) <- NULL;


	if (!is.null(bestRanks) && is.finite(bestRanks))
	{
		wh <- numeric(0);
		for (i in 1:length(f))
		{
			wh <- c(wh, which(rank(-result[,i+1],ties.method="min")<=bestRanks));
		}
		return(result[unique(sort(wh,decreasing=TRUE)),]);
	} else return(result);
}

