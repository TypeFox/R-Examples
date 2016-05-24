##' 19th-century international disputes
##' 
##' Dataset of militarized international disputes between 1816 and 1899.
##'
##' The dataset is taken from the Correlates of War project.  The unit of
##' observation is the dyad-year, and the variables are: \describe{
##' \item{\code{ccode1}}{Initiator's COW country code}
##' \item{\code{ccode2}}{Respondent's COW country code}
##' \item{\code{year}}{Year of dispute}
##' \item{\code{cap_1}}{Initiator's military capabilities (as percent of total
##' system capabilities)}
##' \item{\code{cap_2}}{Respondent's military capabilities (as percent of total
##' system capabilities)}
##' \item{\code{balanc}}{Balance of dyadic capabilities possessed by the
##' initiator (i.e., \code{cap_1 / (cap_1 + cap_2)})}
##' \item{\code{s_wt_re1}}{Dyadic S-score (see Signorino and Ritter 1998),
##' weighted by initiator's region}
##' \item{\code{s_wt_re2}}{Dyadic S-score, weighted by respondent's region}
##' \item{\code{dem1}}{Initiator's Polity score}
##' \item{\code{dem2}}{Respondent's Polity score}
##' \item{\code{distance}}{Distance (in miles) between initiator and respondent}
##' \item{\code{peaceyrs}}{Years since last dispute in this dyad}
##' \item{\code{midnum}}{Dispute's number in the MID data set}
##' \item{\code{revis1}}{Whether the initiator had "revisionist" aims}
##' \item{\code{revis2}}{Whether the respondent had "revisionist" aims}
##' \item{\code{sq}}{Indicator for status quo outcome}
##' \item{\code{capit}}{Indicator for capitulation outcome}
##' \item{\code{war}}{Indicator for war outcome}
##' \item{\code{esc}}{Indicator for escalation (i.e., either capitulation or war
##' occurs)}
##' \item{\code{regime1}}{Initiator's regime type (calculated from \code{dem1})}
##' \item{\code{regime2}}{Respondent's regime type (calculated from \code{dem2})}
##' }
##' @name war1800
##' @usage data(war1800)
##' @docType data
##' @references Daniel M. Jones, Stuart A. Bremer and J. David Singer.  1996.
##' "Militarized Interstate Disputes, 1816-1992: Rationale, Coding Rules, and
##' Empirical Patterns." \emph{Conflict Management and Peace Science} 15(2):
##' 163--213.
##' @seealso \code{\link{egame12}}
##' @keywords data
NULL

##' Currency attacks
##' 
##' Data on speculative currency attacks and devaluation decisions for 90
##' countries from 1985 to 1998.
##'
##' The dataset is taken from Leblang (2003).  The unit of
##' observation is the country-month, and the variables are: \describe{
##' \item{\code{outcome}}{Whether the country faced no speculative attack,
##' defended its currency against an attack, or devalued in response to an
##' attack in the given month}
##' \item{\code{preelec}}{Indicator for being in the three months prior to an
##' election}
##' \item{\code{postelec}}{Indicator for being in the three months following an
##' election}
##' \item{\code{rightgov}}{Indicator for a right-wing party being in power}
##' \item{\code{unifgov}}{Indicator for unified government: in presidential
##' systems, the same party controlling the presidency and the lower house of
##' the legislature; in parliamentary systems, one party/coalition having a
##' majority of seats}
##' \item{\code{lreserves}}{Logged ratio of currency reserves to base money in
##' the previous month}
##' \item{\code{realinterest}}{Domestic real interest rate in the previous
##' month}
##' \item{\code{lexports}}{Logged ratio of exports to GDP in the previous month}
##' \item{\code{capcont}}{Indicator for capital controls in the previous year}
##' \item{\code{overval}}{Overvaluation of the real exchange rate}
##' \item{\code{creditgrow}}{Domestic credit growth in the previous month}
##' \item{\code{service}}{External debt service (as percentage of exports) paid
##' in previous month}
##' \item{\code{USinterest}}{U.S. domestic interest rates in the previous month}
##' \item{\code{contagion}}{Number of other countries experiencing speculative
##' attacks in the same month}
##' \item{\code{prioratt}}{Number of prior speculative attacks experienced by the
##' country}
##' \item{\code{nation}}{Country name}
##' \item{\code{month}}{Month of observation}
##' \item{\code{year}}{Year of observation}
##' }
##' All of the non-binary variables other than \code{nation}, \code{month},
##' and \code{year} are standardized to have mean 0 and unit variance.
##'
##' We are grateful to David Leblang for allowing us to redistribute his data.
##' The original replication file is available in Stata format at
##' \url{https://sites.google.com/site/davidaleblang/data-1} (as of
##' 2015-02-22).
##' @name leblang2003
##' @usage data(leblang2003)
##' @docType data
##' @references David Leblang.  2003.  "To Defend or Devalue: The Political
##' Economy of Exchange Rate Policy."  \emph{International Studies Quarterly}
##' 47: 533--559.
##' @seealso \code{\link{egame12}}
##' @keywords data
##' @example inst/examples/leblang2003.r
NULL

##' Simulated egame122 data
##' 
##' Simulated data for illustrating \code{\link{egame122}}.
##'
##' The variables are: \describe{
##' \item{\code{f1}, \code{f2}}{Factors with levels "a", "b", "c"}
##' \item{\code{x1}--\code{x5}}{Numeric variables entering Player 1's utilities}
##' \item{\code{z1}--\code{z3}}{Numeric variables entering Player 2's utilities}
##' \item{\code{a1}}{Indicator for Player 1's move (L or R)}
##' \item{\code{a2}}{Indicator for Player 2's move (L or R)}
##' \item{\code{y}}{Factor containing outcome}
##' }
##' @name data_122
##' @usage data(data_122)
##' @docType data
##' @seealso \code{\link{egame122}}
##' @keywords data
NULL

##' Simulated egame123 data
##' 
##' Simulated data for illustrating \code{\link{egame123}}.
##'
##' The variables are: \describe{
##' \item{\code{x1}--\code{x8}}{Regressors}
##' \item{\code{a1}}{Indicator for Player 1's move (L or R)}
##' \item{\code{a2}}{Indicator for Player 2's move (L or R)}
##' \item{\code{a3}}{Indicator for Player 3's move (L or R)}
##' \item{\code{y}}{Numeric variable containing outcome number: 1, 3, 5, or 6,
##' corresponding to labels in the game tree in the \code{\link{egame123}}
##' documentation.}
##' }
##' @name data_123
##' @usage data(data_123)
##' @docType data
##' @seealso \code{\link{egame123}}
##' @keywords data
NULL

##' Simulated ultimatum data
##' 
##' Simulated data for illustrating \code{\link{ultimatum}}.
##'
##' The variables are: \describe{
##' \item{\code{offer}}{The offer made by Player 1}
##' \item{\code{accept}}{Whether Player 2 accepted the offer (0 for rejection, 1
##' for acceptance)}
##' \item{\code{w1}, \code{w2}}{Variables entering both players' reservation values}
##' \item{\code{x1}--\code{x4}}{Variables entering Player 1's reservation value}
##' \item{\code{z1}--\code{z4}}{Variables entering Player 2's reservation value}
##' }
##' The maximum offer size is 15.
##' @name data_ult
##' @usage data(data_ult)
##' @docType data
##' @seealso \code{\link{ultimatum}}
##' @keywords data
NULL

##' Data from students playing the ultimatum game
##' 
##' Data from a trial of the ultimatum game with graduate students.
##'
##' The variables are: \describe{
##' \item{\code{offer}}{The offer made by Player 1}
##' \item{\code{accept}}{Whether Player 2 accepted the offer (0 for rejection, 1
##' for acceptance)}
##' \item{\code{gender1}}{Whether Player 1 is female}
##' \item{\code{gender2}}{Whether Player 2 is female}
##' }
##' The maximum offer size is 100.
##' @name student_offers
##' @usage data(student_offers)
##' @docType data
##' @seealso \code{\link{ultimatum}}
##' @keywords data
NULL
