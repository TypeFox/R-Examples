#' Simulate personality scale(s) and an outcome
#'
#' @description
#' Simulates a personality scale which correlates to an outcome. The function can specify
#'  the number of indicators (i.e. indicators) truly relating to the outcome. Also, the
#'   function can create a secondary scale, for instance mimicing informant-report
#' @param n Number of participants
#' @param to_n Number of indicators in a Trait relating to Outcome
#' @param tn_n Number of indicators in a Trait Not relating to outcome.
#' @param indicators2 if TRUE, a secondary set of indicators is created, e.g. to mimic informant-report. Defaults to FALSE
#' @param cor_to_tn Correlation between to and tn. Defaults to 0.3
#' @param cor_to_outcome correlation between to and outcome. Defaults to 0.4
#' @param to_min minimum factor loading for to_n. Defaults to 0.4
#' @param tn_min minimum factor loading for tn_n. Defaults to 0.4
#' @param to_max maximum factor loading for to_n. Defaults to 0.7
#' @param tn_max maximum factor loading for tn_n. Defaults to 0.7
#' @param n.cat number of response options. when you go larger than 5, update the
#' standard deviation as well. Defaults to 5
#' @param sdev standard deviation. Defaults to 0.8
#' @return A list object, first object is indicators' matrix and second object is outcome
#'  vector. If indicators2=TRUE, then a third object is added, which is the  secondary
#'   indicators' matrix.
#' @encoding utf-8

#' @examples
#' ## Create a scale-outcome set that violates ION. Only 2 indicators out of 8 relate
#' ## to the outcome, the others just relate to the 2 indicators This setting is similar
#' ## to the N5: Impulsiveness - BMI association in Vainik et al (2015) EJP paper.
#' set.seed(466)
#' a<-scale_sim(n=2500, to_n=2, tn_n=6)
#' # Last 2 indicators have considerably higher correlation with the outcome
#' cor(a[[1]],a[[2]])
#'
#' ## Create a scale-outcome set that has ION, all 8 indicators relate to the outcome
#' set.seed(466)
#' b<-scale_sim(n=2500, to_n=8)
#' # All indicators correlate largely on the same level with the outcome.
#' cor(b[[1]],b[[2]])
#'
#' ## Create a scale-outcome set that violates ION - only 1 indicator relates to
#' ##the outcome. Include other-report.
#' set.seed(466)
#' c<-scale_sim(n=2500, to_n=1, tn_n=7, indicators2=TRUE)
#' # Last 2 indicators have considerably higher correlation with the outcome
#' cor(c[[1]],c[[2]])
#' cor(c[[3]],c[[2]])
#' @export
#'


scale_sim <- function(n, to_n, tn_n = 0, indicators2 = FALSE, cor_to_tn = 0.3, cor_to_outcome = 0.4, to_min = 0.4,
    to_max = 0.7, tn_min = 0.4, tn_max = 0.7, n.cat = 5, sdev = 0.8) {

    # data generation function
    foo <- function(n, x, min, max, sdev) {
        r <- x * runif(1, min, max) + rnorm(n, sd = sdev)
        return(r)
    }

    n.indicators <- to_n + tn_n

    # generate traito, and outcome related to traito.
    traito <- rnorm(n)

    outcome <- traito * cor_to_outcome + rnorm(n)

    # generate traitn, if needed
    if (tn_n != 0)
        traitn <- traito * cor_to_tn + rnorm(n)

    # preallocate and create the scales

    indicators <- matrix(0, ncol = n.indicators, nrow = n)
    # the indexing will set the out-come related indicators towards the end of the trait. this is helpful when tn
    # exists
    indicators[, (n.indicators - to_n + 1):n.indicators] <- replicate(to_n, foo(n, traito, to_min, to_max, sdev))

    # traitn (related to traito) will only be created if, tn_n!=0
    if (tn_n != 0) {
        indicators[, 1:tn_n] <- replicate(tn_n, foo(n, traitn, tn_min, tn_max, sdev))
    }

    # create an ordinal version - round the values to full numbers, rescore extreme values.

    indicators <- round(indicators, 0)
    lwr <- 0 - (floor(n.cat/2))
    upr <- 0 + (floor(n.cat/2))


    indicators[indicators < lwr] <- lwr
    indicators[indicators > upr] <- upr

    # add another set of indicators, if asked

    if (indicators2 == T) {
        # preallocate and create the scales
        indicators2 <- matrix(0, ncol = n.indicators, nrow = n)
        # the indexing will set the out-come related indicators towards the end of the trait. this is helpful when tn
        # exists
        indicators2[, (n.indicators - to_n + 1):n.indicators] <- replicate(to_n, foo(n, traito, to_min, to_max, sdev))

        # traitn (related to traito) will only be created if, tn_n!=0
        if (tn_n != 0) {
            indicators2[, 1:tn_n] <- replicate(tn_n, foo(n, traitn, tn_min, tn_max, sdev))
        }

        # make them ordinal

        indicators2 <- round(indicators2)

        indicators2[indicators2 < lwr] <- lwr
        indicators2[indicators2 > upr] <- upr

        # indicators2=chopper(indicators2,des_n.cat=n.cat)
        out <- list(indicators, outcome, indicators2)

    } else {
        out <- list(indicators, outcome)
    }
    return(out)
}
