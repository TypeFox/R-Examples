#' Responses from 1983 American National Election Survey Pilot
#'
#' \code{ANES} contains political ideology survey data from the 1983 American National Election Survey Pilot. Respondents were read a value statement and asked 
#' to report their level of agreement: `strongly agree', `agree', `can't decide', `disagree', and strongly disagree'.
#' The 279 complete responses to the 19 statements selected by Feldman (1988)
#' and reanalyzed by Gross and Manrique-Vallier (2014) are included. `Strongly disagree' and `disagree' as well 
#' as `strongly agree' and `agree` have been collapsed into single categories. 
#' 0 indicates `agree', 1 indicates `can't decide', and 2 
#' indicates `disagree'. The statements have been grouped into 3 overarching themes as 
#' indicated in the variable names- Equality (EQ), Economic Individualism (IND), and Free Enterprise (ENT).
#' Variable ID's from the original ANES questionnaire are included in parentheses.
#' An example analysis of the data is included in the 
#' "Fitting Mixed Membership Models Using \code{mixedMem}" vignette. 
#'
#' @name ANES
#' @docType data
#' @usage data(ANES)
#' @format A data frame with 279 individuals and 19 variables:
#' \describe{
#'   \item{EQ1}{If people were treated more equally in this country,
#'    we would have many fewer problems (V832169)}
#'   \item{EQ2}{We should give up on the goal of equality, since people 
#'   are so different to begin with (V832172)}
#'   \item{EQ3}{Our society should do whatever is necessary to make sure 
#'   that everyone has an equal opportunity to succeed (V832175)}
#'   \item{EQ4}{Some people are just better cut out than others for important
#'    positions in society (V832178)}
#'   \item{EQ5}{Some people are better at running things and should be allowed 
#'    to do so (V832250)}
#'   \item{EQ6}{All kinds of people should have an equal say in running
#'    this country, not just those who are successful (V832253)}
#'   \item{EQ7}{One of the big problems in this country is that we don`t
#'    give everyone an equal chance (V832256)}
#'   \item{IND1}{Any person who is willing to work hard has a good chance
#'    of succeeding (V832170)}
#'   \item{IND2}{Hard work offers little guarantee of success (V832173)}
#'   \item{IND3}{Most people who don`t get ahead should not blame 
#'   the system; they really have only themselves to blame (V832176)}
#'   \item{IND4}{Even if people are ambitious, they often cannot succeed (V832251)}
#'   \item{IND5}{If people work hard, they almost always get what they want (V832254)}
#'   \item{IND6}{Even if people try hard, they often cannot reach their goals (V832257)}
#'   \item{ENT1}{The less government gets involved with business and the 
#'   economy, the better off this country will be (V832171)}
#'   \item{ENT2}{There are many goods and services that would never be
#'    available to ordinary people without governmental intervention (V832174)}
#'   \item{ENT3}{There should be no government interference with business
#'    and trade (V832177)}
#'   \item{ENT4}{Putting government regulations on business does not
#'    endanger personal freedom (V832252)}
#'   \item{ENT5}{Government intervention leads to too much red tape and 
#'   too many problems (V832255)}
#'   \item{ENT6}{Contrary to what some people think, a free enterprise
#'    system is not necessary for our form of government to survive (V832258)}
#' }
#' @references
#' Feldman, Stanley. "Structure and consistency in public opinion: The role of core beliefs and values." American Journal of Political Science (1988): 416-440.
#' 
#' Gross, J.H. and Manrique-Vallier, D. "A Mixed-Membership Approach to the Assessment of Political Ideology from Survey Responses." 
#' in Airoldi, E. M., Blei, D. M., Erosheva, E. A., and Fienberg, S. E.
#' Handbook of Mixed Membership Models and Its Applications. Chapman & Hall/CRC, 2014
#' 
#' National Election Studies, 1983 Pilot Election Study. Ann Arbor, MI: University of Michigan, Center for Political Studies, 1999
NULL

#' Point estimates from Gross and Manrique-Vallier 2014 
#'
#' \code{gmv_theta} contains the point estimates for the sub-populaton response probabilities, \eqn{\theta} from the original Gross and
#' Manrique-Vallier 2014 analysis. It is a 3 by 19 by 3 dimensional array with
#' the first dimension corresponding to each latent sub-population, the second dimension corresponding to
#' each variable, and the 3rd dimension corresponding to response choices (Agree, Can't Decide,
#' Disagree). Gross and Manrique-Vallier do not report results for group 3 for reasons discussed in the vignette
#' so all values for group 3 in theta[1, , ] are 0.
#'   
#' @name gmv_theta
#' @docType data
#' @usage data(gmv_theta)
#' 
#' @references
#' Gross, J.H.; Manrique-Vallier, D. "A Mixed-Membership Approach to the Assessment of Political Ideology from Survey Responses." 
#' in Airoldi, E. M., Blei, D. M., Erosheva, E. A., and Fienberg, S. E.
#' Handbook of Mixed Membership Models and Its Applications. Chapman & Hall/CRC, 2014
NULL