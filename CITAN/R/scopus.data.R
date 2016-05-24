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


#' List of Elsevier's \emph{SciVerse Scopus} covered titles (journals, conference proceedings, book series, etc.)
#'
#' Last update: October 2011. The data file is based on the official and publicly available
#' (no permission needed as stated by Elsevier) Scopus list of covered titles,
#' see \url{http://www.info.sciverse.com/documents/files/scopus-training/resourcelibrary/xls/title_list.xls}.
#'
#' This data frame consists of 30794 records.
#' It has the following columns.
#' \tabular{ll}{
#'   \code{SourceId} \tab Unique source identifier in \emph{SciVerse Scopus} (integer). \cr
#'   \code{Title} \tab Title of the source. \cr
#'   \code{Status} \tab Status of the source, either \code{Active} or \code{Inactive}. \cr
#'   \code{SJR_2009} \tab SCImago Journal Rank 2009. \cr
#'   \code{SNIP_2009} \tab Source Normalized Impact per Paper 2009. \cr
#'   \code{SJR_2010} \tab SCImago Journal Rank 2010. \cr
#'   \code{SNIP_2010} \tab Source Normalized Impact per Paper 2010. \cr
#'   \code{SJR_2011} \tab SCImago Journal Rank 2011. \cr
#'   \code{SNIP_2011} \tab Source Normalized Impact per Paper 2011. \cr
#'   \code{OpenAccess} \tab Type of Open Access, see below. \cr
#'   \code{Type} \tab Type of the source, see below. \cr
#'   \code{ASJC} \tab A list of semicolon-separated ASJC classification codes, see \code{\link{Scopus_ASJC}}. \cr
#' }
#'
#' \code{OpenAccess} is one of \code{DOAJ}, \code{Not OA} (not Open Access source),
#' \code{OA but not registered}, \code{OA registered}.
#'
#' \code{Type} is one of \code{Book Series}, \code{Conference Proceedings}, \code{Journal}, \code{Trade Journal}
#'
#' The \code{data.frame} is sorted by \code{Status} (\code{Active} sources first) and then by \code{SJR_2011} (higher values first).
#'
#' @export
#'
#' @title Scopus covered source list
#' @name Scopus_SourceList
#' @docType data
#' @seealso \code{\link{Scopus_ASJC}}, \code{\link{Scopus_ReadCSV}}, \code{\link{Scopus_ImportSources}}
#' @references \url{http://www.info.sciverse.com/scopus/scopus-in-detail/facts/}\cr
#' \url{http://info.scopus.com/journalmetrics/sjr.html}\cr
#' \url{http://info.scopus.com/journalmetrics/snip.html}\cr
#' @keywords Scopus, ASJC, journal, conference, proceedings
Scopus_SourceList <- NULL # will be loaded later


#' List of Elsevier's \emph{SciVerse Scopus} ASJC (All Science. Journals Classification)
#' source classification codes.
#'
#' Last update: October 2011. The data file is based on the official and publicly available
#' (no permission needed as stated by Elsevier) Scopus list of covered titles,
#' see \url{http://www.info.sciverse.com/documents/files/scopus-training/resourcelibrary/xls/title_list.xls}.
#'
#' It consists of 334 ASJC 4-digit integer codes (column \code{ASJC})
#' together with their group identifiers (column \code{ASJC_Parent})
#' and descriptions (column \code{Description}).
#'
#' ASJC codes are used to classify Scopus sources (see \code{\link{Scopus_SourceList}}).
#'
#' @export
#'
#' @title Scopus ASJC (All Science. Journals Classification) classification codes
#' @name Scopus_ASJC
#' @docType data
#' @seealso \code{\link{Scopus_SourceList}}, \code{\link{Scopus_ReadCSV}}, \code{\link{Scopus_ImportSources}}
#' @references \url{http://www.info.sciverse.com/scopus/scopus-in-detail/facts/}
#' @keywords Scopus, ASJC, journal
Scopus_ASJC <- NULL # will be loaded later
