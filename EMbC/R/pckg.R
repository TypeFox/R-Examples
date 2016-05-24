#' @title Expectation-Maximization binary Clustering package.
#'   
#' @description
#' 
#' The Expectation-maximization binary clustering (EMbC) is a general purpose, 
#' unsupervised, multi-variate, clustering algorithm, driven by two main 
#' motivations: (i) it looks for a good compromise between statistical soundness
#' and ease and generality of use - by minimizing prior assumptions and favouring
#' the semantic interpretation of the final clustering - and, (ii) it allows 
#' taking into account the uncertainty in the data. These two fetaures make it
#' specially suitable for the behavioural annotation of animal's movement 
#' trajectories.
#'
#' @details
#'
#' The method is a variant of the well sounded Expectation-Maximization
#' Clustering (EMC) algorithm, - i.e. under the assumption of an underlying
#' Gaussian Mixture Model (GMM) describing the distribution of the data-set - but 
#' constrained to generate a binary partition of the input space. This is achieved by
#' means of the *delimiters*, a set of parameters that discretizes the input features
#' into high and low values and define the binary regions of the input space.
#' As a result, each final cluster includes a unique  combination of either
#' low or high values of the input variables. Splitting the input features
#' into low and high values is what favours the semantic interpretation of the final
#' clustering.
#' 
#' The  initial assumptions implemented in the EMbC algorithm aim at minimizing biases
#' and sensitivity to initial conditions: (i) each data point is assigned a uniform 
#' probability of belonging to each cluster, (ii) the  prior mixture distribution  is
#' uniform (each  cluster starts with the  same number of  data points), (iii) the
#' starting partition, (*i.e.* initial delimiters position),  is selected  based
#' on  a global maximum  variance criterion,  thus conveying  the minimum  information
#' possible.
#'
#' The number of output clusters is $2^m$ determined by the number of input features
#' $m$. This number is only an upper bound as some of the clusters can be merged along
#' the likelihood optimization process. The EMbC algorithm is intented to be used with
#' not more than 5 or 6 input features, yielding a maximum of 32 or 64 clusters. This
#' limitation in the number of clusters is consistent with the main motivation of the
#' algorithm of favouring the semantic interpretation of the results.
#'
#' The algorithm deals very intuitively with data reliability: the larger the
#' uncertainty associated with a data point, the  smaller the leverage of that data
#' point in the clustering.
#'
#' Compared to close related methods like EMC and Hidden Markov Models (HMM),
#' the EMbC is specially useful when: (i) we can expect bi-modality, to some extent,
#' in the conditional distribution of the input features or, at least, we can assume
#' that a binary partition of the input space can provide useful information, and
#' (ii) a first order temporal dependence assumption, a necessary condition in HMM,
#' can not be guaranteed.
#'
#' The EMbC R-package is mainly intended for the behavioural annotation of animals'
#' movement trajectories where an easy interpretation of the final clustering and the
#' reliability of the data constitute two key issues, and the conditions of bi-modality
#' and unfair temporal dependence usually hold. In particular, the temporal dependence
#' condition is easily violated in animal's movement trajectories because of the
#' heterogeneity in empirical time series due to large gaps, or prefixed sampling
#' scheduling.
#'
#' Input movement trajectories are given either as a *data.frame* or a *Move* object
#' from the **move** R-package.  The pakage deals also with stacks of trajectories for
#' population level analysis. Segmentation is based on local estimates of velocity and
#' turning angle, eventually including a solar position covariate as a daytime
#' indicator.
#'
#' The core clustering method is complemented with a set of functions to easily
#' visualize and analyze the output:
#'  
#'  * clustering statistics,
#'  * clustering scatterplot (2D and 3D)
#'  * temporal labeling profile (ethogram),
#'  * plotting of intermediate variables,
#'  * confusion matrix (numerical validation with respect to an expert's labelling),
#'  * visual validation with external information (e.g. environmental data),
#'  * generation of kml or webmap docs for detailed inspection of the output.
#'  
#' Also, some functions are provided to further refine the output, either by
#' pre-processing (smoothing) the input data or by post-processing (smoothing,
#' relabelling, merging) the output labelling.
#'
#' The results obtained for different empirical datasets suggest that the EMbC
#' algorithm  behaves reasonably well for a wide range  of tracking technologies,
#' species, and ecological contexts (e.g. migration, foraging).
#'
#'
#' @rdname EMbC_pckg
#' @name EMbC-package
#' @aliases EMbC
#' @docType package
#' @title Expectation-Maximization binary Clustering.
#' @author Joan Garriga \email{jgarriga@@ceab.csic.es}
#'
#' @keywords package
#' @seealso \code{\link{move}}
NULL

#' binClst Instance definition
#' 
#' Unless otherwise specified, a \code{binClst} instance refers to any of the binary clustering objects defined in the package, either a \code{binClst} object itself, or any of its child classes, a \code{binClstPath} or a \code{binClstMove} instance.  The latter inherit all slots and functionality defined for the former.
#' 
#' @rdname binClst_Instance
#' @name binClst_instance
NULL

#' binClstPath Instance definition
#' 
#' Unless otherwise specified, a \code{binClstPath} instance refers to a \code{binClstPath} object itself, as well as its child class \code{binClstMove}. The latter inherits all slots and functionality defined for the former.
#' 
#' @rdname binClstPath_Instance
#' @name binClstPath_instance
NULL

#' Synthetic 2D object used in the examples
#' 
#' An ad-hoc object with a set of bivariate data points synthetically generated by sampling from a four component GMM and their corresponding labels indicating which component of the mixture generated each data point.
#' 
#' @rdname x2d
#' @name x2d
#' @docType data
#' @format See parameter \code{X} of the \link{embc} constructor.
#'
#' @keywords datasets
NULL

#' Synthetic path used in the examples
#' 
#' A data.frame with a synthetically generated trajectory with colum values (timeStamps, longituds, latitudes, labels) and column headers ('dTm','lon','lat','lbl'). The order of the columns is important. Column headers can be whatever but are expected to be there. The only exception is the header for the labels column: if headed as 'lbl' it will be used automatically by any methods that can make use of it.
#' 
#' @rdname expth
#' @name expth
#' @docType data
#' @format See parameter \code{pth} of the \link{stbc} constructor.
#'
#' @keywords datasets
NULL
