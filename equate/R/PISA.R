#' Programme for International Student Assessment 2009 USA Data
#' 
#' This dataset contains scored cognitive item response data from the 2009
#' administration of the Programme for International Student Assessment (PISA),
#' an international study education systems. The data, along with the license
#' under which they are released, are available online at
#' \url{http://www.oecd.org/pisa/}.
#' 
#' @name PISA
#' @docType data
#' @format \code{PISA} is a \code{list} containing four elements. The first,
#' \code{PISA$students}, is a \code{data.frame} containing 233 variables across
#' 5233 individuals, with one row per individual. All but one variable come
#' from the USA PISA data file "INT_COG09_S_DEC11.txt". The remaining variable,
#' language spoken at home, has been merged in from the student questionnaire
#' file "INT_STQ09_DEC11.txt". Variable names match those found in the original
#' files: \describe{ \item{list("stidstd")}{ Unique student ID (one for each of
#' the 5233 cases); } \item{list("schoolid")}{ School ID (there are 165
#' different schools); } \item{list("bookid")}{ ID for the test booklet given
#' to a particular student, of which there were 13; } \item{list("langn")}{
#' Student-reported language spoken at home, with 4466 students reporting
#' English (indicated by code 313), 484 students reporting Spanish (with code
#' 156) and 185 students reporting "another language" (code 859); }
#' \item{list("m033q01")}{ Scored item-response data across the 189 items
#' included in the general cognitive assessment, described below; and }\item{
#' to }{ Scored item-response data across the 189 items included in the general
#' cognitive assessment, described below; and }\item{list("s527q04t")}{ Scored
#' item-response data across the 189 items included in the general cognitive
#' assessment, described below; and } \item{list("pv1math")}{ PISA scale
#' scores, referred to in the PISA technical documentation as "plausible
#' values".  }\item{ to }{ PISA scale scores, referred to in the PISA technical
#' documentation as "plausible values".  }\item{list("pv5read5")}{ PISA scale
#' scores, referred to in the PISA technical documentation as "plausible
#' values".  } } Next, \code{PISA$booklets} is a \code{data.frame} containing 4
#' columns and 756 rows and describes the 13 general cognitive assessment
#' booklets. Variables include: \describe{ \item{list("bookid")}{ The test
#' booklet ID, as in \code{PISA$students}; } \item{list("clusterid")}{ ID for
#' the cluster or item subset in which an item was placed; items were fully
#' nested within clusters; however, each item cluster appeared in four
#' different test booklets; } \item{list("itemid")}{ Item ID, matching the
#' columns of \code{PISA$students}; each item appears in \code{PISA$booklets}
#' four times, once for each booklet; and } \item{list("order")}{ The order in
#' which the cluster was presented within a given booklet.  } }
#' \code{PISA$items} is a \code{data.frame} containing 4 columns and 189 rows,
#' with one row per item. Variables include: \describe{ \item{list("itemid")}{
#' Item ID, as in \code{PISA$booklets} } \item{list("clusterid")}{ Cluster ID,
#' as in \code{PISA$booklets} } \item{list("max")}{ Maximum possible score
#' value, either 1 or 2 points, with dichotomous scoring (max of 1) used for
#' the majority of items; and } \item{list("subject")}{ The subject of an item,
#' equivalent to the first character in \code{itemid} and \code{clusterid}.  }
#' \item{list("format")}{ Item format, abbreviated as \code{mc} for multiple
#' choice, \code{cmc} for complex multiple choice, \code{ocr} for open
#' constructed response, and \code{ccr} for closed constructed response.  }
#' \item{list("noptions")}{ Number of options, zero except for some multiple
#' choice items.  } } Finally, \code{PISA$totals} is a list of 13
#' \code{data.frame}s, one per booklet, where the columns correspond to total
#' scores for all students on each cluster for the corresponding booklet. These
#' total scores were calculated using \code{PISA$students} and
#' \code{PISA$booklets}. Elements within the \code{PISA$totals} list are named
#' by booklet, and the columns in the \code{data.frame} are named by cluster.
#' For example, \code{PISA$totals$b1$m1} contains the total scores on cluster
#' M1 for students taking booklet 1.
#' @source OECD (2012). PISA 2009 Technical Report, PISA, OECD Publishing.
#' http://dx.doi.org/10.1787/9789264167872-en
#' 
#' Addition information can be found at the PISA website:
#' \url{http://www.oecd.org/pisa/}
#' @keywords datasets
"PISA"