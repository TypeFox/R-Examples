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


#' \pkg{CITAN} is a library of functions useful in --- but not limited to ---
#' quantitative research in the field of scientometrics.
#' It contains various tools for preprocessing
#' bibliographic data retrieved from, e.g., Elsevier's \emph{SciVerse Scopus} and
#' computing bibliometric impact of individuals.
#' Moreover, some functions dealing with Pareto-Type II (GPD)
#' and Discretized Pareto-Type II statistical models
#' are included (e.g., Zhang-Stephens and MLE estimators,
#' goodness-of-fit and two-sample tests, confidence intervals
#' for the theoretical Hirsch index etc.).
#' They may be used to describe and analyze many phenomena encountered
#' in the social sciences.
#'
#'
#' Fair and objective assessment methods of individual scientists
#' had become the focus of scientometricians' attention since the
#' very beginning of their discipline. A quantitative expression
#' of some publication-citation process'
#' characteristics is assumed to be a predictor of broadly conceived
#' scientific competence. It may be used e.g. in building decision support
#' systems for scientific quality control.
#'
#' The \eqn{h}-index, proposed by J.E. Hirsch (2005)
#' is among the most popular scientific impact indicators.
#' An author who has published \eqn{n}
#' papers has the Hirsch index equal to \eqn{H}, if each of his \eqn{H}
#' publications were cited at least \eqn{H} times, and each of the
#' remaining \eqn{n-H} items were cited no more than \eqn{H} times. This simple
#' bibliometric tool quickly received much attention in the academic community
#' and started to be a subject of intensive
#' research. It was noted that, contrary to earlier approaches,
#' i.e. publication count, citation count, etc.,
#' this measure  concerns both productivity and impact of an
#' individual.
#' \cr
#'
#' In a broader perspective, this issue is a special case
#' of the so-called \dfn{Producer Assessment Problem}
#' (PAP; see Gagolewski, Grzegorzewski, 2010b).
#'
#' Consider a \emph{producer} (e.g. a writer, scientist, artist,
#' craftsman) and a nonempty set of his \emph{products} (e.g. books,
#' papers, works, goods). Suppose that each product is given a
#' \emph{rating} (of quality, popularity, etc.) which is a single
#' number in \eqn{I=[a,b]}, where \eqn{a} denotes the lowest admissible
#' valuation. We typically choose \eqn{I=[0,\infty]} (an interval
#' in the extended real line).
#' Some instances of the PAP are listed below.
#'
#' \tabular{cllll}{
#'   \tab \strong{Producer}    \tab \strong{Products}   \tab \strong{Rating method} \tab \strong{Discipline}\cr
#' A \tab Scientist            \tab Scientific articles \tab Number of citations    \tab Scientometrics\cr
#' B \tab Scientific institute \tab Scientists          \tab The h-index            \tab Scientometrics\cr
#' C \tab Web server           \tab Web pages           \tab Number of in-links     \tab Webometrics\cr
#' D \tab Artist               \tab Paintings           \tab Auction price          \tab Auctions\cr
#' E \tab Billboard company    \tab Advertisements      \tab Sale results           \tab Marketing\cr
#' }
#'
#' Each possible state of producer's activity can therefore be represented by a point
#' \eqn{x\in I^n} for some \eqn{n}. Our aim is thus to construct
#' and analyze --- both theoretically and empirically ---
#' aggregation operators (cf. Grabisch et al, 2009) which can be used for rating
#' producers. A family of such functions should take  the two
#' following aspects of producer's quality into account:
#' \itemize{
#'     \item the ability to make highly-rated products,
#'     \item overall productivity, \eqn{n}.
#' }
#' For some more formal considerations please refer to (Gagolewski, Grzegorzewski, 2011).
#' \cr\cr
#'
#'
#'
#'
#' To \bold{preprocess and analyze bibliometric data} (cf. Gagolewski, 2011) retrieved
#' from e.g.  Elsevier's \emph{SciVerse Scopus}
#' we need the \pkg{RSQLite} package. It is an interface to the free
#' SQLite DataBase Management System (see \url{http://www.sqlite.org/}).
#' All data is stored in a so-called Local Bibliometric Storage (\acronym{LBS}),
#' created with the \code{\link{lbsCreate}} function.
#'
#' The data frames \code{\link{Scopus_ASJC}} and \code{\link{Scopus_SourceList}}
#' contain various information on current source coverage of SciVerse Scopus.
#' They may be needed during the creation of the LBS and \code{\link{lbsCreate}}
#' for more details.
#' \emph{License information: this data are publicly available
#'       and hence no special permission is needed to redistribute them
#'       (information from Elsevier).}
#'
#' \pkg{CITAN} is able to import publication data from Scopus CSV files
#' (saved with settings "Output: complete format" or "Output: Citations only",
#' see \code{\link{Scopus_ReadCSV}}). Note that the output limit in Scopus
#' is 2000 entries per file. Therefore, to perform
#' bibliometric research we often need to divide the query results
#' into many parts. \pkg{CITAN} is able to merge them back even if
#' records are repeated.
#'
#' The data may be accessed via functions from the \pkg{DBI} interface.
#' However, some typical tasks may be automated using
#' e.g. \code{\link{lbsDescriptiveStats}} (basic description of the whole sample
#' or its subsets, called \sQuote{Surveys}),
#' \code{\link{lbsGetCitations}} (gather citation sequences selected
#' authors), and \code{\link{lbsAssess}} (mass-compute impact functions'
#' values for given citation sequences).
#'
#' There are also some helpful functions (in **EXPERIMENTAL** stage) which use
#' the \pkg{RGtk2} library (see Lawrence, Lang, 2010)
#' to display some suggestions on which documents or authors should be
#' merged, see \code{\link{lbsFindDuplicateTitles}} and
#' \code{\link{lbsFindDuplicateAuthors}}.
#'
#'
#' For a complete list of functions, call \code{library(help="CITAN")}.
#' \cr\cr
#'
#' \bold{Keywords}: Hirsch's h-index, Egghe's g-index, L-statistics,
#' S-statistics, bibliometrics, scientometrics, informetrics,
#' webometrics, aggregation operators, arity-monotonicity,
#' impact functions, impact assessment.
#'
#' @name CITAN-package
#' @aliases CITAN
#' @docType package
#' @title CITation ANalysis toolpack
#' @author Marek Gagolewski
#' @references
#' GTK+ Project, \url{http://www.gtk.org}\cr
#' SQLite DBMS, \url{http://www.sqlite.org/}\cr
#' Dubois D., Prade H., Testemale C. (1988). Weighted fuzzy pattern matching,
#'    Fuzzy Sets and Systems 28, s. 313-331.\cr
#' Egghe L. (2006). Theory and practise of the g-index, Scientometrics 69(1),
#'    131-152.\cr
#' Gagolewski M., Grzegorzewski P. (2009). A geometric approach to the construction of
#'    scientific impact indices, Scientometrics 81(3), 617-634.\cr
#' Gagolewski M., Debski M., Nowakiewicz M. (2009). Efficient algorithms for computing
#'    ''geometric'' scientific impact indices, Research Report of Systems
#'    Research Institute, Polish Academy of Sciences RB/1/2009.\cr
#' Gagolewski M., Grzegorzewski P. (2010a). S-statistics and their basic properties,
#'    In: Borgelt C. et al (Eds.), Combining Soft Computing and Statistical
#'    Methods in Data Analysis, Springer-Verlag, 281-288.\cr
#' Gagolewski M., Grzegorzewski P. (2010b). Arity-monotonic extended aggregation
#'    operators, In: Hullermeier E., Kruse R., Hoffmann F. (Eds.),
#'    Information Processing and Management of Uncertainty in Knowledge-Based
#'    Systems, CCIS 80, Springer-Verlag, 693-702.\cr
#' Gagolewski M. (2011). Bibliometric Impact Assessment with R and the CITAN Package,
#' Journal of Informetrics 5(4), 678-692.\cr
#' Gagolewski M., Grzegorzewski P. (2011a). Axiomatic Characterizations of (quasi-)
#' L-statistics and S-statistics and the Producer Assessment Problem,
#; In: Galichet S., Montero J., Mauris G. (Eds.), Proc. 7th conf. European Society
#' for Fuzzy Logic and Technology (EUSFLAT/LFA 2011), Atlantic Press, 53-58.
#' Grabisch M., Pap E., Marichal J.-L., Mesiar R. (2009). Aggregation functions,
#'    Cambridge.\cr
#' Gagolewski M., Grzegorzewski P. (2011b). Possibilistic analysis of arity-monotonic
#'    aggregation operators and its relation to bibliometric impact assessment
#'    of individuals, International Journal of Approximate Reasoning 52(9), 1312-1324.\cr
#' Hirsch J.E. (2005). An index to quantify individual's scientific research output,
#'    Proceedings of the National Academy of Sciences 102(46),
#'    16569-16572.\cr
#' Kosmulski M. (2007). MAXPROD - A new index for assessment of the scientific output
#'    of an individual, and a comparison with the h-index, Cybermetrics 11(1).\cr
#' Lawrence M., Lang D.T. (2010). RGtk2: A graphical user interface toolkit for R,
#'    Journal of Statistical Software 37(8), 1-52.\cr
#' Woeginger G.J. (2008). An axiomatic characterization of the Hirsch-index,
#'    Mathematical Social Sciences 56(2), 224-232.\cr
#' Zhang J., Stevens M.A. (2009). A New and Efficient Estimation Method for the
#'    Generalized Pareto Distribution, Technometrics 51(3), 316-325.\cr
#' @importFrom stringi stri_trim_both
#' @importFrom stringi stri_replace_all_fixed
#' @importFrom hash hash
#' @importFrom hash clear
#' @importFrom RGtk2 gtkWindowNew
#' @importFrom RGtk2 gtkLabelNew
#' @importFrom RGtk2 gtkProgressBarNew
#' @importFrom RGtk2 gtkProgressBarSetText
#' @importFrom RGtk2 gtkProgressBarSetText
#' @importFrom RGtk2 gtkProgressBarSetFraction
#' @importFrom RGtk2 gObjectSetData
#' @importFrom RGtk2 gObjectGetData
#' @importFrom RGtk2 gtkWidgetQueueDraw
#' @importFrom RGtk2 gtkLabelSetText
#' @importFrom RGtk2 gtkWidgetDraw
#' @importFrom RGtk2 gtkDialogNewWithButtons
#' @importFrom RGtk2 GtkResponseType
#' @importFrom RGtk2 rGtkDataFrame
#' @importFrom RGtk2 gtkTreeView
#' @importFrom RGtk2 gtkVBox
#' @importFrom RGtk2 gtkTreeViewColumn
#' @importFrom RGtk2 gtkCellRendererToggle
#' @importFrom RGtk2 gSignalConnect
#' @importFrom RGtk2 gtkCellRendererText
#' @importFrom RGtk2 GtkWrapMode
#' @importFrom RGtk2 gtkCheckVersion
#' @importFrom RGtk2 gtkScrolledWindow
#' @importFrom RSQLite dbGetInfo
#' @importFrom RSQLite dbGetQuery
#' @importFrom RSQLite dbCommit
#' @importFrom DBI dbDriver
#' @importFrom RSQLite dbConnect
#' @importFrom DBI dbDisconnect
#' @importFrom RSQLite dbListTables
#' @importFrom grDevices as.graphicsAnnot
#' @importFrom grDevices dev.interactive
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics barplot
#' @importFrom graphics boxplot
#' @importFrom graphics mtext
#' @importFrom graphics par
#' @importFrom graphics pie
#' @importFrom stats na.omit
#' @importFrom utils read.csv
invisible(NULL)
