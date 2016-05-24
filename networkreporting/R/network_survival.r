#####################################################
## network_survival.r
##
## the network survival method for estimating the
## size and characteristics of a hidden population
##

#####################################################
##' network survival estimator
##'
##' use an aggregate multiplicity estimator
##' and the respondents' own network size estimates
##' to estimate hidden population sizes 
##'
##' This function takes two sources of data as input: first, it requires
##' a long-form dataframe with the attributes of the reported members of
##' the hidden population. For example, if we are asking about emigres and
##' we collect the age and sex of each reported emigrant, then the long form
##' dataset might look like:
##' \tabular{ccc}{
##'   age \tab sex \tab weight\cr
##'    15 \tab m   \tab 2.10  \cr
##'    58 \tab f   \tab 1.15  \cr
##'    33 \tab m   \tab 3.67  \cr
##' }
##'  
##'  The second source of data we need is the known population responses for the
##'  respondents, along with the *same* attributes for each respondent. For
##'  example, in the situation above, we would also require a dataset like this
##'  to be passed in
##'  \tabular{cccccc}{
##'   age  \tab sex  \tab  weight \tab hm.teachers \tab  hm.nurses  \tab ...  \cr
##'   20   \tab   f  \tab  2.10   \tab 4           \tab  0          \tab ...  \cr
##'   44   \tab   m  \tab  1.65   \tab 0           \tab  2          \tab ...  \cr
##'   60   \tab   m  \tab  2.75   \tab 1           \tab  1          \tab ...  \cr
##' }
##'
##' @section Technical note:
##' This function assumes that the sampling weights are standard analysis weights
##' and *not* relative weights. Standard analysis weights should provide an estimate
##' for the size of the frame population when added up; relative weights, on the other hand,
##' will sum to the number of respondents in the sample.  Demographic and Health surveys
##' typically have relative weights, which must be converted into standard sampling weights
##' before using this function.
##'
##' @section TODO:
##' \itemize{
##' \item{ handle missing values}
##' \item{ think about whether or not this is the best way to handle N.F}
##' \item{ write more general agg mult est fn and call that }
##' \item{ make unit tests }
##' }
##' 
##' @param resp.data the dataframe that has a row for each respondent, with reported
##'                  connections to the groups of known size, as well as the attributes.
##'                  Note that the column names of the attributes should match
##'                  their names in \code{attribute.data}
##' @param attribute.data A dataframe with the reported attributes of hidden population 
##'                       members reported by survey respondents. There should be one 
##'                       row for each time a respondent reports a hidden population member.
##'                       For example, to estimate death rates, there should be one row for
##'                       each report of a death.
##' @param attribute.names the names of the columns of attribute.data
##'                        and resp.data that contain the attribute information.
##' @param known.populations the names of the columns in \code{resp.data} that
##'          have responses to the known population questions
##' @param total.kp.size the size of the probe alters, i.e., the sum of the known
##'         population sizes
##' @param weights the weights or weights column for the respondent data
##' @param attribute.weights the weights or weights column for the alter data
##' @param dropmiss see \code{\link{report.aggregator}}
##' @param verbose if TRUE, print information to screen
##' @return the network reporting estimate of the hidden population's size
##'         (as a prevalence) broken down by the categories defined by all combinations
##'         of \code{attribute.names}.
##'
##' @rdname network.survival.estimator
##' @export
network.survival.estimator_ <- function(resp.data,
                                        attribute.data,
                                        attribute.names,
                                        known.populations,
                                        total.kp.size=1,
                                        weights,
                                        attribute.weights,
                                        dropmiss=NULL,
                                        verbose=TRUE) {

    ## estimate the average personal network size of the respondents
    ## for each combination of attributes
    deg.by.att <- kp.estimator_(resp.data=resp.data,
                                known.populations=known.populations,
                                attribute.names=attribute.names,
                                weights=weights,
                                total.kp.size=total.kp.size)

    ## count the number of connections from respondents to members
    ## of the hidden population with each combination of attributes
    attribute.data$death <- 1
    deaths.by.att <- report.aggregator_(attribute.data,
                                        attribute.names,
                                        "death",
                                        attribute.weights,
                                        qoi.name="deaths",
                                        dropmiss)

    tog.df <- plyr::join(deg.by.att, deaths.by.att, by=attribute.names)

    ## NB: we're using the sum of the sampling weights as N.F
    N.F <- sum(tog.df$wgt.total.y.kp)

    surveybootstrap:::vcat(verbose,
         "Taking N.F value implied by weights: ", N.F, "\n")

    ## to placate R CMD CHECK
    sum.deaths <- NULL
    sum.y.kp.over.kptot <- NULL

    ## for each combination of attributes, divide the number of
    ## reported connections by the estimated network size to produce
    ## an estimate of the overall size. then, divide this by the
    ## group size to obtain a rate.
    ## (these are done implicitly here because some of the factors cancel --
    ##  see TODO vignette)
    tog.df <- tog.df %>%
              dplyr::mutate(asdr.hat = sum.deaths / (sum.y.kp.over.kptot * N.F))

    return(tog.df)

}

#####################################################
##' @rdname network.survival.estimator
network.survival.estimator <- function(resp.data,
                                       attribute.data,
                                       attribute.names,
                                       known.populations,
                                       total.kp.size=1,
                                       weights,
                                       attribute.weights,
                                       verbose=TRUE) {

    network.survival.estimator_(resp.data,
                                attribute.data,
                                lazy(attribute.names),
                                lazy(known.populations),
                                lazy(total.kp.size),
                                lazy(weights),
                                lazy(attribute.weights),
                                lazy(verbose))

}

