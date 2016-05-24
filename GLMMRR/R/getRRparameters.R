#' Compute Randomized Response parameters
#'
#' Represent Randomize Response models with two parameters, c and d. Reference:
#' Fox, J-P, Klotzke, K. and Veen, D. (2016). \emph{Generalized Linear Mixed Models for Randomized
#' Responses.} Manuscript submitted for publication.
#'
#'
#' @param vec.RRmodel
#' a character vector of Randomized Response models.
#' @param vec.p1
#' a numeric vector of p1 values.
#' @param vec.p2
#' a numeric vector of p2 values.
#'
#' @return
#' a list with c and d values.
#' @export
getRRparameters <- function(vec.RRmodel, vec.p1, vec.p2)
{
  getCandD <- function(RRmodel, p1, p2)
  {
    if (RRmodel == "DQ" | RRmodel == "Warner" | RRmodel == "Forced" | RRmodel == "UQM" | RRmodel == "Crosswise" |
        RRmodel == "Triangular" | RRmodel == "Kuk") {

      #%% Direct questioning, no RR influence on the data
      if (RRmodel == "DQ") {
        # Specifying c and d values for direct questioning
        c <- 0
        d <- 1
      }

      #%% Warner, p1 is meant to be the proportion that answers the "normal" question
      if (RRmodel == "Warner") {
        # Specifying c and d values for the warner model
        c <- (1 - p1)
        d <- (2 * p1 - 1)

        # Restrictions
        if (!all(p1 >= 0.5 & p1 <= 1)) {
          stop("P1 value should be between 0.5 and 1.
               p1 is meant to be the proportion that answers question that has an affirmative answer as sensitive characteristic.")
        }
      }

      #%% Forced, p1 is meant to be the proportion that answers the sensitive question
      #% p2 is meant to be the propotion of affirmative forced responses
      if (RRmodel == "Forced") {
        # Specifying c and d for Forced model
        d <- p1
        c <- (1 - p1) * p2

        # Restrictions
        if (!all(p1 >= 0.5 & p1 <= 1 & c >= 0 & c <= 1)) {
          stop("P1 value should be between 0.5 and 1.
               p1 is meant to be the proportion that answers the sensitive question.
               p2 is meant to be the propotion of affirmative forced responses.")
        }

        }

      ##%% Unrelated Question Method, p1 is meant to be the proportion that answers the sensitive question
      #% p2 is meant to be the proportion that answers the unsensitive question affirmative
      if (RRmodel == "UQM") {
        # Specifying c and d for UQM
        d <- p1
        c <- (1 - p1) * p2

        # Restrictions
        if (!all(p1 >= 0.5 & p1 <= 1 & d >= 0 & d <= 1)) {
          stop("P1 value should be between 0.5 and 1.
               p1 is meant to be the proportion that answers the sensitive question.
               p2 is meant to be the propotion of affirmative responses to the non sensitive question.")
        }
        }

      #%% Crosswise model, p1 is meant to be the propotion that answers the non sensitive question affirmative
      if (RRmodel == "Crosswise"){
        # Specifying c and d for Crosswise model
        c <- 1 - p1
        d <- (2 * p1 - 1)

        # Restrictions
        if (!all(p1 >= 0 & p1 <= 1)){
          stop("p1 should be between 0 and 1.
               p1 is meant to be the proportion that posesses non sensitive characteristic")
        }
        }

      #%% Triangular model, p1 is meant to be the proportion that answers the non sensitive question affirmative
      if (RRmodel == "Triangular"){
        # Specifying c and d for Triangular model
        c <- p1
        d <- (1 - p1)

        # Restrictions
        if (!all(p1 >= 0 & p1 <= 1)){
          stop("p1 should be between 0 and 1.
               p1 is meant to be the proportion that posesses non sensitive characteristic")
        }
        }

      if (RRmodel == "Kuk") {
        # Specifying c and d for Kuk model
        c <- p2
        d <- (p1 - p2)

        # Restrictions
        if (!all(p1 >= 0 & p1 <= 1 & p2 >= 0 & p2 <= 1)) {
          stop("p1 and p2 should be values between 0 and 1
               p1 is meant to be the binomial proportion to get an affirmitive answer when the respondent posesses the sensitive characteristic
               p2 is meant to be the binomial proportion to get an affirmative answer when the respondent does not posess the sensitive characteristic")
        }
      }

      return(list("c"= c, "d" = d))

    }
    else {
          stop("RRmodel should be one of the following; DQ, Warner, Forced, UQM, Crosswise, Triangular or Kuk.")
    }
  }

  c <- rep(-1, length(vec.RRmodel))
  d <- rep(-1, length(vec.RRmodel))

  for(ii in 1:length(vec.RRmodel))
  {
    parameters <- getCandD(vec.RRmodel[ii], vec.p1[ii], vec.p2[ii])
    c[ii] <- parameters$c
    d[ii] <- parameters$d
  }

  return(list("c"= c, "d" = d))

}
