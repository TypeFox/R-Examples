#' @import methods
#' @title RankData Class
#' @description A S4 class to represent ranking data
#' 
#' It is well understood that the ranking representation and ordering representation of ranking data can easily be confused.
#' I thus use a S4 class to store all the information about the ranking data. This can avoid unnecessary confusion.
#' 
#' @slot nobj The number of ranked objects. If not provided, it will be inferred as the maximum ranking in the data set. As a result, it must be provided if the data is top-q ranking.
#' @slot nobs the number of observations. No need to be provided during initialization since it must be equal to the sum of slot \code{count}.
#' @slot ndistinct the number of distinct rankings. No need to be provided during initialization since it must be equal to the number of rows of slot \code{ranking}.
#' @slot ranking a matrix that stores the ranking representation of distinct rankings. Each row contains one ranking. For top-q ranking, all unobserved objects have ranking \code{q+1}.
#' @slot count the number of observations for each distinct ranking corresponding to each row of \code{ranking}.
#' @slot topq a numeric vector to store top-q ranking information. More information in details section.
#' @slot subobs a numeric vector to store number of observations for each chunk of top-q rankings.
#' @slot q_ind a numeric vector to store the beginning and ending of each chunk of top-q rankings. The last element has to be \code{ndistinct+1}.
#' @details 
#' It is possible to store both complete and top-q rankings in the same RankData object. Three slots \code{topq}, \code{subobs}, and
#' \code{q_ind} are introduced for this purpose. Note that there is generally no need to specify these slots if your data set only contains 
#' a single "q" level (for example all data are top-10 rankings). The "q" level for complete ranking should be \code{nobj-1}. 
#' Moreover, if the rankings are organized in chunks of increasing "q" levels (for
#'  example, top-2 rankings followed by top-3 rankings followed by top-5 rankings etc.), then slots \code{subobs}, and \code{q_ind} can also be inferred 
#'  correctly by the initializer. Therefore it is highly recommender that you organise the ranking matrix in this way and utilize the initializer.
#' @examples
#' # creating a data set with only complete rankings
#' rankmat <- replicate(10,sample(1:52,52), simplify = "array")
#' countvec <- sample(1:52,52,replace=TRUE)
#' rankdat <- new("RankData",ranking=rankmat,count=countvec)
#' # creating a data set with both complete and top-10 rankings
#' rankmat_in <- replicate(10,sample(1:52,52), simplify = "array")
#' rankmat_in[rankmat_in>11] <- 11
#' rankmat_total <- cbind(rankmat_in, rankmat)
#' countvec_total <- c(countvec,countvec)
#' rankdat2 <- new("RankData",ranking=rankmat_total,count=countvec_total, nobj=52, topq=c(10,51))
#' @aliases RankData RankData-class
#' @seealso \code{\link{RankInit}}, \code{\link{RankControl}}
#' @export
setClass( "RankData",
          representation = representation(
              nobj = "numeric",
              nobs = "numeric",
              ndistinct = "numeric",
              ranking = "matrix",
              count = "numeric",
              topq = "numeric", # a numeric vector: topq 
			        subobs = "numeric", # number of observations for each topq
			        q_ind = "numeric"  # starting point of each topq location
          ),
          prototype = prototype(
              topq = -1,
			        subobs = 1,
			        q_ind = 1
          )
)


#' @title RankInit Class
#' @description A S4 class to store initialization information of model fitting
#' 
#' The \code{RankInit} class is used to give initial values of model fitting procedures.
#' 
#' @slot param.init a list containing initial values of the positive parametrization of weights.
#' @slot modal_ranking.init a list containing starting points for the modal ranking search.
#' @slot clu an integer containing the number of clusters used in the model.
#' @slot p.init a numeric vector containing the initial values for cluster probabilities.
#' @examples c1init = new("RankInit",param.init=list(rep(1,4)),
#'      modal_ranking.init=list(c(2,3,4,1,5)),clu=1L)
#' c2init = new("RankInit",param.init=list(rep(0.1,4),rep(0.1,4)),
#'      modal_ranking.init = list(c(2,3,4,1,5),c(2,5,1,4,3)),clu=2L,p.init=c(0.5,0.5))
#' @aliases RankInit RankInit-class
#' @seealso \code{\link{RankData}}, \code{\link{RankControl}}
#' @export
setClass( "RankInit",
          representation = representation(
              param.init = "list",
              modal_ranking.init = "list",
              clu = "integer",
              p.init = "numeric"
          ),
          prototype = prototype(
              clu = 1L,
              p.init = 1
          )
)

#' @title RankControl Class
#' @description A virtual S4 class to store control parameters for model fitting. 
#' 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current best. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours
#' @details RankControl class must be extended to reflect what distance metric should be used. Possibles extensions are \code{\link{RankControlWeightedKendall}}, \code{\link{RankControlKendall}}, \code{\link{RankControlPhiComponent}}, 
#' \code{\link{RankControlWtau}}, \code{\link{RankControlSpearman}}, \code{\link{RankControlFootrule}}, \code{\link{RankControlHamming}}, and \code{\link{RankControlCayley}}.
#' 
#' The control parameters that start with prefix \code{EM_} are intended for the EM iteration. The ones with prefix \code{SeachPi0} control the behaviour of searching model ranking.
#' @section User-defined Criterion:
#' You can specify user-defined criterion to choose modal rankings. The function object SearchPi0_FUN takes a list as argument. The components in the list include the following. \code{obs}: the number of observations.
#' \code{w.est}: the estimated weights. \code{log_likelihood}: the estimated log_likelihood. With this information, most of the popular information criterion can be supported and customized criterion can also be defined.
#' A larger returned value indicates a better fit. Note that if you are fitting a mixture model the EM algorithm always tries to maximized the log likelihood. Thus the default value should be used in this case.
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}
#' @aliases RankControl RankControl-class
#' @export
setClass( "RankControl",
          representation = representation(
              EM_limit = "numeric",
              EM_epsilon = "numeric",
              SearchPi0_limit = "numeric",
              SearchPi0_FUN = "function",
              SearchPi0_fast_traversal = "logical",
              SearchPi0_show_message = "logical",
              SearchPi0_neighbour = "character",
              # optimx_control = "list",
              "VIRTUAL"
          ),
          prototype = prototype(
              EM_limit=500,
              EM_epsilon=1e-5,
              SearchPi0_limit=500,
              SearchPi0_FUN = function(x){x$log_likelihood},
              SearchPi0_fast_traversal=TRUE,
              SearchPi0_show_message=FALSE,
              SearchPi0_neighbour="Cayley"
              # optimx_control = list(maximize=FALSE,starttests=TRUE,trace=0,dowarn=FALSE)
          )
)


#' @title RankControlWeightedKendall Class
#' @description A S4 class to store control parameters for Weighted Kendall distance model fitting. It is derived from class \code{\link{RankControl-class}}. 
#' 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours
#' @slot optimx_control a list to be passed to \code{\link[optimx]{optimx}}. The list must not contain a component \code{maximize=TRUE} since internally the negation of the likelihood function is minimized.
#' @slot assumption A character string specifying which assumption to use when handling top-q rankings. Supported choices are "equal-probability" and "tied-rank".
#' @details \code{RankControlWeightedKendall} is derived from virtual class \code{\link{RankControl}}. All slots in \code{\link{RankControl}} are still valid.
#' This control class tells the solver to fit a model based on Weighted Kendall distance.  
#' The control parameters that start with prefix \code{EM_} are intended for the EM iteration. The ones with prefix \code{SeachPi0} control the behaviour of searching model ranking.
#' @examples # enabling  warnings
#' testctrl = new("RankControlWeightedKendall",optimx_control=list(dowarn=TRUE))
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlWeightedKendall RankControlWeightedKendall-class
#' @export
setClass( "RankControlWeightedKendall",
        contains = "RankControl",
        representation = representation(
            optimx_control = "list",
            assumption = "character"
        ),
        prototype = prototype(
            optimx_control = list(maximize=FALSE,starttests=TRUE,trace=0,dowarn=FALSE),
            assumption = "tied-rank" # equal-probability
        )
)

#' @title RankControlKendall Class
#' @description A S4 class to store control parameters for Kendall distance model fitting (Mallow's Phi Model).
#' It is derived from class \code{\link{RankControl-class}}. 
#' 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours
#' @details \code{RankControlKendall} is derived from virtual class \code{\link{RankControl}}.
#' This control class tells the solver to fit a model based on Kendall distance.  
#' The control parameters that start with prefix \code{EM_} are intended for the EM iteration. The ones with prefix \code{SeachPi0} control the behaviour of searching model ranking.
#' @examples # enabling messages
#' testctrl = new("RankControlKendall",SearchPi0_show_message=TRUE)
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlKendall RankControlKendall-class
#' @export
setClass( "RankControlKendall",
        contains = "RankControl"
)

#' @title RankControlPhiComponent Class
#' @description A S4 class to store control parameters for Phi component model fitting.
#' It is derived from class \code{\link{RankControl-class}}. 
#' 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours
#' @details \code{RankControlKendall} is derived from virtual class \code{\link{RankControl}}.
#' This control class tells the solver to fit a model based on a stage-wise generalization of Kendall distance.  
#' The control parameters that start with prefix \code{EM_} are intended for the EM iteration. The ones with prefix \code{SeachPi0} control the behaviour of searching model ranking.
#' @examples # enabling messages
#' testctrl = new("RankControlPhiComponent",SearchPi0_show_message=TRUE)
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlPhiComponent RankControlPhiComponent-class
#' @export
setClass( "RankControlPhiComponent",
        contains = "RankControl"
)

#' @title RankControlWtau Class
#' @description A S4 class for the Weighted tau model fitting.
#' It is derived from class \code{\link{RankControl-class}}. 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours

#' @slot optimx_control a list to be passed to \code{\link[optimx]{optimx}}. The list must not contain a component \code{maximize=TRUE} since internally the negation of the likelihood function is minimized.
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlWtau RankControlWtau-class
#' @export
setClass( "RankControlWtau",
          contains = "RankControl",
          representation = representation(
              optimx_control = "list"
          ),
          prototype = prototype(
              optimx_control = list(maximize=FALSE,starttests=TRUE,trace=0,dowarn=FALSE)
          )
)

#' @title RankControlSpearman Class
#' @description A S4 class for the Spearman distance model fitting.
#' It is derived from class \code{\link{RankControl-class}}. 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours

#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlSpearman RankControlSpearman-class
#' @export
setClass("RankControlSpearman",
         contains = "RankControl"
)

#' @title RankControlFootrule Class
#' @description A S4 class for the Footrule distance model fitting.
#' It is derived from class \code{\link{RankControl-class}}. 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours

#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlFootrule RankControlFootrule-class
#' @export
setClass("RankControlFootrule",
         contains = "RankControl"
)

#' @title RankControlHamming Class
#' @description A S4 class for the Hamming distance model fitting.
#' It is derived from class \code{\link{RankControl-class}}. 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours

#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlHamming RankControlHamming-class
#' @export
setClass("RankControlHamming",
         contains = "RankControl"
)

#' @title RankControlCayley Class
#' @description A S4 class for the Cayley distance model fitting.
#' It is derived from class \code{\link{RankControl-class}}. 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local search of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours

#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlCayley RankControlCayley-class
#' @export
setClass("RankControlCayley",
         contains = "RankControl"
)


