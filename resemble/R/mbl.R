#' @title A function for memory-based learning (mbl)
#' @description
#' This function is implemented for memory-based learning (a.k.a. instance-based learning or local regression) which is a non-linear lazy learning approach 
#' for predicting a given response variable from a set of (spectral) predictor variables. For each sample in an prediction set  a specific local 
#' regression is carried out based on a subset of similar samples (nearest neighbours) selected from a reference set. The local model is 
#' then used to predict the response value of the target (prediction) sample. Therefore this function does not yield a global regression model. 
#' @usage
#' mbl(Yr, Xr, Yu = NULL, Xu,
#'     mblCtrl = mblControl(), 
#'     dissimilarityM,
#'     group = NULL,
#'     dissUsage = "predictors", 
#'     k, k.diss, k.range,
#'     method, 
#'     pls.c, pls.max.iter = 1, pls.tol = 1e-6,
#'     noise.v = 0.001,
#'     ...)
#' @param Xr input \code{matrix} (or \code{data.frame}) of predictor variables of the reference data (observations in rows and variables in columns). 
#' @param Yr a numeric \code{vector} containing the values of the response variable corresponding to the reference data
#' @param Xu input \code{matrix} (or \code{data.frame}) of predictor variables of the data to be predicted (observations in rows and variables in columns). 
#' @param Yu an optional numeric \code{vector} containing the values of the response variable corresponding to the data to be predicted
#' @param mblCtrl a list (created with the \code{\link{mblControl}} function) which contains some parameters that control the some aspects of the \code{mbl} function. See the \code{\link{mblControl}} function for more details. 
#' @param dissimilarityM (optional) a dissimilarity matrix. This argument can be used in case a user-defined dissimilarity matrix is preferred over the automatic dissimilarity matrix computation specified in the \code{sm} argument of the  \code{\link{mblControl}} function. When \code{dissUsage = "predictors"}, \code{dissimilarityM} must be a square symmetric dissimilarity matrix (derived from a matrix of the form \code{rbind(Xr, Xu)}) for which the diagonal values are zeros (since the dissimilarity between an object and itself must be 0). 
#'        On the other hand if \code{dissUsage} is set to either \code{"weights"} or \code{"none"}, \code{dissimilarityM} must be a \code{matrix} representing the dissimilarity of each element in \code{Xu} to each element in \code{Xr}. The number of columns of the object correspondent to \code{dissimilarityM} must be equal to the number of rows in \code{Xu} and the number of rows equal to the number of rows in \code{Xr}.
#'        If both \code{dissimilarityM} and  \code{sm} are specified, only the \code{dissimilarityM} argument will be taken into account.
#' @param group an optional \code{factor} (or \code{vector} that can be coerced to a \code{factor} by \code{as.factor}) to be taken into account for internal validations. The length of the \code{vector} must be equal to \code{nrow(Xr)}, giving the identifier of related observations (e.g. spectra collected from the same batch of measurements, from the same sample, from samples with very similar origin, etc). When one observation is selected for cross validation, all observations of the same group are removed together and assigned to validation. See details.
#' @param dissUsage specifies how the dissimilarity information shall be used. The possible options are: \code{"predictors"}, \code{"weights"} and \code{"none"} (see details below).
#'        Default is "predictors".
#' @param k a numeric (integer) \code{vector} containing the sequence of k nearest neighbours to be tested. Either \code{k} or \code{k.diss} must be specified. Numbers with decimal values will be coerced to their next higher integer values. This \code{vector} will be automatically sorted into ascending order.
#' @param k.diss a \code{vector} containing the sequence of dissimilarity thresholds to be tested. When the dissimilarity between a sample in \code{Xr} and a sample \code{Xu} is below the given threshold, the sample in sample in \code{Xr} is treated as a neighbour of the sample in \code{Xu}, otherwise it is ignored.  These thresholds depend
#'        on the corresponding dissimilarity measure specified in \code{sm}. Either \code{k} or \code{k.diss} must be specified.
#' @param k.range a \code{vector} of length 2 which specifies the minimum (first value) and the maximum (second value) number of neighbours allowed when the \code{k.diss} argument is used. 
#' @param method a character string indicating the method to be used at each local multivariate regression. Options are: \code{"pls"}, \code{"wapls1"} and \code{"gpr"} (see details below). Note: \code{"wapls2"} from the previos version of the package is no longer available/supported.
#' @param pls.c the number of pls components to be used in the regressions if either \code{"pls"} or \code{"wapls1"} is used.
#'        When \code{"pls"} is used, this argument must be a single numerical value. When \code{"wapls1"} is used, this argument must be
#'        a \code{vector} of length 2 indicating the minimum (first value) and the maximum (second value) number of pls components
#'        used for the regressions (see details below).
#' @param pls.max.iter maximum number of iterations for the partial least squares methods.
#' @param pls.tol limit for convergence in the partial orthogonal scores partial least squares regressions using the nipals algorithm. Default is 1e-6.
#' @param noise.v a value indicating the variance of the noise for Gaussian process regression. Default is 0.001. 
#' @param ... additional arguments to be passed to other functions. 
#' @details
#' By using the \code{group} argument one can specify observations (spectra) groups of samples that have something in common e.g. spectra collected from the same batch of measurements, from the same sample, from samples with very similar origin, etc) which could produce biased cross-validation results due to pseudo-replication. This argument allows to select calibration points that are independent from the validation ones in cross-validation. In this regard, when  \code{valMethod = "loc_crossval"} (used in \code{\link{mblControl}} function), then the \code{p} argument refer to the percentage of groups of samples (rather than single samples) to be retained in each resampling iteration at each local segment. 
#' The \code{dissUsage} argument is used to specifiy whether the dissimilarity information must be used within the local regressions and (if so), how. When \code{dissUsage = "predictors"}
#' the local (square symmetric) dissimilarity matrix corresponding the selected neighbourhood is used as source
#' of additional predictors (i.e the columns of this local matrix are treated as predictor variables). In some cases this may result in an improvement 
#' of the prediction performance (Ramirez-Lopez et al., 2013a). If \code{dissUsage = "weights"}, the neighbours of the query point (\eqn{xu_{j}}) are weighted according to their dissimilarity (e.g. distance) to \eqn{xu_{j}} prior carrying out each local regression. 
#' The following tricubic function (Cleveland and Delvin, 1988; Naes et al., 1990) is used for computing the final weights based on the measured dissimilarities:
#' \deqn{W_{j}  =  (1 - v^{3})^{3}}
#' where if \eqn{xr_{i} \in} neighbours of \eqn{xu_{j}}:
#' \deqn{v_{j}(xu_{j})  =  d(xr_{i}, xu_{j})}
#' otherwise:
#' \deqn{v_{j}(xu_{j})  =  0}
#' In the above formulas \eqn{d(xr_{i}, xu_{j})} represents the dissimilarity between the query point and each object in \eqn{Xr}. When \code{dissUsage = "none"} is chosen
#' the dissimilarity information is not used. 
#' The possible options for performing regressions at each local segment implemented in the \code{mbl} function are described as follows: \itemize{
#'  \item{Partial least squares (\code{"pls"}):}{ It uses the orthogonal scores (non-linear iterative partial least squares, nipals) algorithm. The only parameter which needs to be optimized is the number of pls components. This can be done by cross-validation at each local segment.}
#'  \item{Weighted average pls (\code{"wapls1"}):}{ It uses multiple models generated by multiple pls components (i.e. between a minimum and a maximum number of pls components). At each local partition the final predicted value is a weighted average of all the predicted values generated by the multiple pls models. The weight for each component is calculated as follows: 
#'  \deqn{
#'        w_{j}  =  \frac{1}{s_{1:j}\times g_{j}}
#'        }
#'  where \eqn{s_{1:j}} is the root mean square of the spectral residuals of the unknown (or target) sample when a total of \eqn{j} pls components are used and \eqn{g_{j}} is the root mean square of the regression coefficients corresponding to the \eqn{j}th pls component (see Shenk et al., 1997 for more details). \code{"wapls1"} is not compatible with \code{valMethod = "loc_crossval"} since the weights are computed based on the sample to be predicted at each local iteration.}
#'  \item{Gaussian process with dot product covariance (\code{"gpr"}):}{ Gaussian process regression is a probabilistic and non-parametric Bayesian approach. It is commonly described as a collection of random variables which have a joint Gaussian distribution and it is characterized by both a mean and a covariance function (Williams and Rasmussen, 1996). The covariance function used in the implemented method is the dot product, which inplies that there are no parameters to be optimized for the computation of the covariance.
#'  Here, the process for predicting the response variable of a new sample (\eqn{y_{new}}) from its predictor variables (\eqn{x_{new}}) is carried out first by computing a prediction vector (\eqn{A}). It is derived from a set of reference spectra (\eqn{X}) and their respective response vector (\eqn{Y}) as follows: 
#'   \deqn{
#'       A = (X  X^\textup{T} + \sigma^2 I)^{-1} Y
#'    }
#'   where  \eqn{\sigma^2} denotes the variance of the noise and \eqn{I} the identity matrix (with dimensions equal to the number of observations in \eqn{X}). The prediction of \eqn{y_{new}} is then carried out by:
#'   \deqn{
#'     y_{new} = (x_{new}x_{new}^\textup{T}) A
#'     }
#'     }
#'  }
#' The loop used to iterate over the \code{Xu} samples in \code{mbl} uses the \code{\%dopar\%} operator of the \code{\link{foreach}} package, which can be used to parallelize this internal loop. The last example given in the \code{\link{mbl}} function ilustrates how to parallelize the \code{\link{mbl}} function.
#' Note that the computational cost depends largely on the way on which the arguments of the function are set. For big datasets, it is recommended to carefully select the values of the parameters to test (e.g. validation method, regression method, dissimilarity information usage, dissimilarity method).
#' @return a \code{list} of class \code{mbl} with the following components (sorted by either \code{k} or \code{k.diss} according to the case):
#' \itemize{
#'  \item{\code{call}:}{ the call used.}
#'  \item{\code{cntrlParam}:}{ the list with the control parameters used. If one or more control parameters were reset automatically, then a list containing a list with the initial control parameters specified and a list with the parameters which were finally used.}
#'  \item{\code{dissimilarities}:}{ a list with the method used to obtain the dissimilarity matrices and the dissimilarity matrices corresponding to \eqn{D(Xr, Xu)} and \eqn{D(Xr,Xr)} if \code{dissUsage = "predictors"}. This object is returned only if the \code{returnDiss} argument in the \code{mblCtrl} list was set to \code{TRUE} in the the call used.}
#'  \item{\code{totalSamplesPredicted}}{ the total number of samples predicted.}
#'  \item{\code{pcAnalysis}:}{ a list containing the results of the principal component analysis. The first two objects (\code{scores_Xr} and \code{scores_Xu}) are the scores of the \code{Xr} and \code{Xu} matrices. It also contains the number of principal components used (\code{n.componentsUsed}) and another object which is  a \code{vector} containing the standardized Mahalanobis dissimilarities (also called GH, Global H distance) between each sample in \code{Xu} and the centre of \code{Xr}.}
#'  \item{\code{components}:}{ a list containing either the number of  principal components  or partial least squares components used for the computation of the orthogonal dissimilarities. This object is only returned if the dissimilarity meausre specified in \code{mblCtrl} is any of the following options: \code{'pc'}, \code{'loc.pc'}, \code{"pls"}, \code{'loc.pls'}. If any of the local orthogonal dissimilarities was used (\code{'loc.pc'} or \code{"pls"})
#'                            a \code{data.frame} is also returned in his list. This object is equivalent to the \code{loc.n.components} object returned by the \code{\link{orthoDiss}} function. It specifies the number of local components (either principal components or partial least squares components) used for computing the dissimilarity between each query sample and its neighbour samples, as returned by the \code{\link{orthoDiss}} function.
#'                            }
#'  \item{\code{nnValStats}:}{ a data frame containing the statistics of the nearest neighbour cross-validation for each either \code{k} or \code{k.diss} depending on the arguments specified in the call. It is returned only if \code{'NNv'} or \code{'both'} were selected as validation method}
#'  \item{\code{localCrossValStats}:}{ a data frame containing the statistics of the local leave-group-out cross validation for each either \code{k} or \code{k.diss} depending on the arguments specified in the call. It is returned only if \code{'local_crossval'} or \code{'both'} were selected as validation method}
#'  \item{\code{YuPredictionStats}:}{ a data frame containing the statistics of the cross-validation of the prediction of  \code{Yu} for each either \code{k} or \code{k.diss} depending on the arguments specified in the call. It is returned only if \code{Yu} was provided.}
#'  \item{\code{results}:}{ a list of data frames which contains the results of the predictions for each either \code{k} or \code{k.diss}. Each \code{data.frame} contains the following columns:}
#'  \itemize{
#'  \item{\code{o.index}:}{ The index of the sample predicted in the input \code{matrix}}
#'  \item{\code{k.diss}:}{ This column is only ouput if the \code{k.diss} argument is used. It indicates the corresponding dissimilarity threshold for selecting the neighbors used to predict a given sample.}
#'  \item{\code{distance}:}{ This column is only ouput if the \code{k.diss} argument is used. It is a logical that indicates whether the neighbors selected by the given dissimilarity threshold were outside the boundaries specified in the \code{k.range} argument. In that case the number of neighbors used is coerced to on of the boundaries.}
#'  \item{\code{k.org}:}{ This column is only ouput if the \code{k.diss} argument is used. It indicates the number of neighbors that are retained when the given dissimilarity threshold is used.}
#'  \item{\code{pls.comp}:}{ This column is only ouput if \code{pls} regression was used. It indicates the final number of pls components used. If no optimization was set, it retrieves the original pls components specified in the \code{pls.c} argument.}
#'  \item{\code{min.pls}:}{ This column is only ouput if \code{wapls1} regression was used. It indicates the final number of minimum pls components used. If no optimization was set, it retrieves the original minimum pls components specified in the \code{pls.c} argument.}
#'  \item{\code{max.pls}:}{ This column is only ouput if \code{wapls1} regression was used. It indicates the final number of maximum pls components used. If no optimization was set, it retrieves the original maximum pls components specified in the \code{pls.c} argument.}
#'  \item{\code{yu.obs}:}{ This column is only ouput if the \code{Yu} argument is used. It indicates the input values given in \code{Yu} (the response variable corresponding to the data to be predicted).}     
#'  \item{\code{pred}:}{ The predicted values}     
#'  \item{\code{yr.min.obs}:}{ The minimum reference value (of the response variable) in the neighborhood.}
#'  \item{\code{yr.max.obs}:}{ The maximum reference value (of the response variable) in the neighborhood.}
#'  \item{\code{index.nearest.in.ref}}{ The index in \code{Xr} of the nearest neighbor.}
#'  \item{\code{y.nearest}:}{ The reference value (of the response variable) of the nearest neighbor in \code{Xr}.}
#'  \item{\code{y.nearest.pred}:}{ This column is only ouput if the validation method (selected with the \code{\link{mblControl}} function) is equal to \code{'NNv'}. It represents the predicted value of the nearest neighbor sample in \code{Xr} using the neighborhood of the predicted sample in \code{Xu}.}
#'  \item{\code{loc.rmse.cv}:}{ This column is only ouput if the validation method (selected with the \code{\link{mblControl}} function) is equal to \code{'loc_crossval'}. It represents the cross validation RMSE value computed in for the neighborhood of the sample of the predicted sample in \code{Xu}.}
#'  \item{\code{loc.st.rmse.cv}:}{ This column is only ouput if the validation method (selected with the \code{\link{mblControl}} function) is equal to \code{'loc_crossval'}. It represents the cross validation standardized RMSE value computed in for the neighborhood of the sample of the predicted sample in \code{Xu}.}
#'  \item{\code{dist.nearest}:}{ The distance to the nearest neighbor.}
#'  \item{\code{dist.k.farthest}:}{ The distance to the farthest neighbor selected.}
#'  }
#'  }
#' When the \code{k.diss} argument is used, the printed results show a table with a column named 'p.bounded'. It represents the percentage of samples
#' for which the neighbors selected by the given dissimilarity threshold were outside the boundaries specified in the \code{k.range} argument.   
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @references 
#' Cleveland, W. S., and Devlin, S. J. 1988. Locally weighted regression: an approach to regression analysis by local fitting. Journal of the American Statistical Association, 83, 596-610.
#' 
#' Fernandez Pierna, J.A., Dardenne, P. 2008. Soil parameter quantification by NIRS as a Chemometric challenge at "Chimiomitrie 2006". Chemometrics and Intelligent Laboratory Systems 91, 94-98
#' 
#' Naes, T., Isaksson, T., Kowalski, B. 1990. Locally weighted regression and scatter correction for near-infrared reflectance data. Analytical Chemistry 62, 664-673.  
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Stevens, A., Dematte, J.A.M., Scholten, T. 2013a. The spectrum-based learner: A new local approach for modeling soil vis-NIR spectra of complex datasets. Geoderma 195-196, 268-279.
#' 
#' Ramirez-Lopez, L., Behrens, T., Schmidt, K., Viscarra Rossel, R., Dematte, J. A. M.,  Scholten, T. 2013b. Distance and similarity-search metrics for use with soil vis-NIR spectra. Geoderma 199, 43-53.
#'
#' Rasmussen, C.E., Williams, C.K. Gaussian Processes for Machine Learning. Massachusetts Institute of Technology: MIT-Press, 2006.
#'   
#' Shenk, J., Westerhaus, M., and Berzaghi, P. 1997. Investigation of a LOCAL calibration procedure for near infrared instruments. Journal of Near Infrared Spectroscopy, 5, 223-232.
#' 
#' @seealso \code{\link{mblControl}}, \code{\link{fDiss}}, \code{\link{corDiss}}, \code{\link{sid}}, \code{\link{orthoDiss}}, \code{\link{neigCleaning}}
#' @examples
#' \dontrun{
#' require(prospectr)
#' 
#' data(NIRsoil)
#' 
#' # Filter the data using the Savitzky and Golay smoothing filter with 
#' # a window size of 11 spectral variables and a polynomial order of 3 
#' # (no differentiation).
#' sg <- savitzkyGolay(NIRsoil$spc, p = 3, w = 11, m = 0) 
#'
#' # Replace the original spectra with the filtered ones
#' NIRsoil$spc <- sg
#' 
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' 
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' Xu <- Xu[!is.na(Yu),]
#' Xr <- Xr[!is.na(Yr),]
#' 
#' Yu <- Yu[!is.na(Yu)]
#' Yr <- Yr[!is.na(Yr)]
#' 
#' # Example 1
#' # A mbl implemented in Ramirez-Lopez et al. (2013, 
#' # the spectrum-based learner)
#' # Example 1.1
#' # An exmaple where Yu is supposed to be unknown, but the Xu 
#' # (spectral variables) are known 
#' ctrl1 <- mblControl(sm = "pc", pcSelection = list("opc", 40), 
#'                     valMethod = "NNv", 
#'                     scaled = FALSE, center = TRUE)
#' 
#' sbl.u <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
#'              mblCtrl = ctrl1, 
#'              dissUsage = "predictors", 
#'              k = seq(40, 150, by = 10), 
#'              method = "gpr")
#' sbl.u
#' plot(sbl.u)
#' 
#' 
#'  
#' # Example 1.2
#' # If Yu is actually known... 
#' sbl.u2 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'               mblCtrl = ctrl1, 
#'               dissUsage = "predictors", 
#'               k = seq(40, 150, by = 10), 
#'               method = "gpr")
#' sbl.u2
#'
#' # Example 1.3
#' # A variation of the spectrum-based learner implemented in 
#' # Ramirez-Lopez et al. (2013) where the dissimilarity matrices are 
#' # recomputed based on partial least squares scores
#' ctrl_1.3 <- mblControl(sm = "pls", pcSelection = list("opc", 40), 
#'                        valMethod = "NNv", 
#'                        scaled = FALSE, center = TRUE)
#'                           
#' sbl_1.3 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'                mblCtrl = ctrl_1.3,
#'                dissUsage = "predictors",
#'                k = seq(40, 150, by = 10), 
#'                method = "gpr")
#' sbl_1.3
#' 
#' # Example 2
#' # A mbl similar to the ones implemented in 
#' # Ramirez-Lopez et al. (2013) 
#' # and Fernandez Pierna and Dardenne (2008)
#' ctrl.mbl <- mblControl(sm = "cor", 
#'                        pcSelection = list("cumvar", 0.999), 
#'                        valMethod = "NNv", 
#'                        scaled = FALSE, center = TRUE)
#'                           
#' local.mbl <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'                  mblCtrl = ctrl.mbl,
#'                  dissUsage = "none",
#'                  k = seq(40, 150, by = 10), 
#'                  pls.c = c(5, 15),
#'                  method = "wapls1")
#' local.mbl
#'
#' # Example 3
#' # A variation of the previous example (using the optimized pc 
#' # dissmilarity matrix) using the control list of the example 1
#'                          
#' local.mbl2 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'                   mblCtrl = ctrl1,
#'                   dissUsage = "none",
#'                   k = seq(40, 150, by = 10), 
#'                   pls.c = c(5, 15),
#'                   method = "wapls1")
#' local.mbl2
#' 
#' # Example 4
#' # Using the function with user-defined dissimilarities:
#' # Examples 4.1 - 4.2: Compute a square symetric matrix of 
#' # dissimilarities between 
#' # all the elements in Xr and Xu (dissimilarities will be used as 
#' # additional predictor variables later in the mbl function)
#' 
#' # Examples 4.3 - 4.4: Derive a dissimilarity value of each element 
#' # in Xu to each element in Xr (in this case dissimilarities will 
#' # not be used as additional predictor variables later in the 
#' # mbl function)
#' 
#' # Example 4.1
#' # the manhattan distance 
#' manhattanD <- dist(rbind(Xr, Xu), method = "manhattan") 
#' manhattanD <- as.matrix(manhattanD)
#'
#' ctrl.udd <- mblControl(sm = "none", 
#'                        pcSelection = list("cumvar", 0.999), 
#'                        valMethod = c("NNv", "loc_crossval"), 
#'                        resampling = 10, p = 0.75,
#'                        scaled = FALSE, center = TRUE)
#'
#' mbl.udd1 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'                 mblCtrl = ctrl.udd, 
#'                 dissimilarityM = manhattanD,
#'                 dissUsage = "predictors",
#'                 k = seq(40, 150, by = 10), 
#'                 method = "gpr")
#' mbl.udd1
#' 
#' #Example 4.2
#' # first derivative spectra
#' Xr.der.sp <- t(diff(t(rbind(Xr, Xu)), lag = 7, differences = 1)) 
#' Xu.der.sp <- t(diff(t(Xu), lag = 7, differences = 1)) 
#' 
#' # The principal components dissimilarity on the derivative spectra 
#' der.ortho <- orthoDiss(Xr = Xr.der.sp, X2 = Xu.der.sp,
#'                        Yr = Yr,
#'                        pcSelection = list("opc", 40),
#'                        method = "pls",
#'                        center = FALSE, scale = FALSE) 
#' 
#' der.ortho.diss <- der.ortho$dissimilarity
#' 
#' # mbl applied to the absorbance spectra
#' mbl.udd2 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'                 mblCtrl = ctrl.udd, 
#'                 dissimilarityM = der.ortho.diss,
#'                 dissUsage = "none",
#'                 k = seq(40, 150, by = 10), 
#'                 method = "gpr")
#'                                 
#' #Example 4.3
#' # first derivative spectra
#' der.Xr <- t(diff(t(Xr), lag = 1, differences = 1)) 
#' der.Xu <- t(diff(t(Xu), lag = 1, differences = 1))
#' # the sid on the derivative spectra
#' der.sid <- sid(Xr = der.Xr, X2 = der.Xu, mode = "density", 
#'                center = TRUE, scaled = FALSE) 
#' der.sid <- der.sid$sid
#' 
#' mbl.udd3 <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'                 mblCtrl = ctrl.udd, 
#'                 dissimilarityM = der.sid,
#'                 dissUsage = "none",
#'                 k = seq(40, 150, by = 10), 
#'                 method = "gpr")
#' mbl.udd3
#' 
#' # Example 5
#' # For running the mbl function in parallel
#' n.cores <- detectCores() - 1
#' if(n.cores == 0) n.cores <- 1
#' 
#' # Set the number of cores according to the OS
#' if (.Platform$OS.type == "windows") {
#'   require(doParallel)
#'   clust <- makeCluster(n.cores)   
#'   registerDoParallel(clust)
#' }else{
#'   require(doSNOW)
#'   clust <- makeCluster(n.cores, type = "SOCK")
#'   registerDoSNOW(clust)
#'   ncores <- getDoParWorkers()
#' }
#' 
#' ctrl <- mblControl(sm = "pc", pcSelection = list("opc", 40), 
#'                    valMethod = "NNv",
#'                    scaled = FALSE, center = TRUE)
#'
#' mbl.p <- mbl(Yr = Yr, Xr = Xr, Yu = Yu, Xu = Xu,
#'              mblCtrl = ctrl, 
#'              dissUsage = "none",
#'              k = seq(40, 150, by = 10), 
#'              method = "gpr")
#' registerDoSEQ()
#' try(stopCluster(clust))
#' mbl.p
#' }
#' @export

######################################################################
# resemble
# Copyrigth (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
######################################################################

## History:
## 09.03.2014 Leo     Doc examples  were formated with a max. line width
## 13.03.2014 Antoine The explanation of the cores argument was modified   
## 23.04.2014 Leo     Added default variable names when they are missing 
## 08.09.2014 Leo     A bug related with the computations of the weights
##                    for wapls2 was fixed
## 23.09.2014 Leo     A bug that prevented the mbl function of using 
##                    the 'dissimilarityM' argument was fixed
## 03.10.2014 Antoine Fix bug when scale = T and add allowParallel argument
## 12.10.2014 Leo     noise.v was missing in the locFit function used for 
##                    the nearest neighbor validation
## 16.11.2015 Leo     Now the scale argument for gaussian process is a 
##                    indicates that both x and y variables must be scaled 
##                    to zero mean and unit variance. Before it only the x 
##                    variables were scaled to unit variance
## 18.11.2015 Leo     The mbl examples were corrected (unnecessary arguments 
##                    were deleted)
## 01.12.2015 Leo     The wapls2 was removed from the options of methods
## 10.12.2015 Leo     The locFit function has been renamed to locFitnpred
##                    and now it always performs a prediction.
## 10.12.2015 Leo     Several redundant/repaeated sanity checks (ocurring 
##                    when combining functions) were deactivated. This does not
##                    impact the finaly sanity checks of the overall mbl fucntion.
## 11.12.2015 Leo     A bug when the k.diss argument was used was corrected. 
##                    The results about the percentage of samples that were bounded 
##                    by k.range was not not correct.  
## 11.12.2015 Leo     The explanation of the output variables in the results element of 
##                    the mbl objects was extended. The rep variable is not output anymore
##                    in the results element.
## 03.01.2016 Leo     Now it is possible to optimize the max and min pls components of wapls1 
## 04.02.2016 Leo     An extrange bug was fixed. The object predobs (in the parallel loop) had a
##                    variable named pls.c (predobs$pls.c). When when method = "gpr" was used, 
##                    and mbl was runing in parallel it retrieved and error saying that pls.c was missing!!!
##                    This was perhaps due to the fact that the pls.c  was variable (in predobs) and an argument.
## 16.02.2016 Leo     Bug fixed. It caused the mbl functio to return an error (sometimes) when the group argument was used together 
##                    with local cross validation. The errors ocurred when groups containing very few samples (e.g. 1 or 2) were used.
                   
mbl <- function(Yr, Xr, Yu = NULL, Xu,
                mblCtrl = mblControl(),
                dissimilarityM,
                group = NULL,
                dissUsage = "predictors",
                k, 
                k.diss,
                k.range,
                method,
                pls.c,
                pls.max.iter = 1,
                pls.tol = 1e-6,
                noise.v = 0.001,                
                ...){
  
  '%mydo%' <- if (mblCtrl$allowParallel) get('%dopar%') else get('%do%')    
  
  ini.cntrl <- mblCtrl
  
  # Sanity checks
  if(ncol(Xr) != ncol(Xu))
    stop("The number of predictor variables in Xr must be equal to the number of variables in Xu")
  
  if(ncol(Xr) < 4)
    stop("This function works only with matrices containing more than 3 predictor variables")
  
  if(length(Yr) != nrow(Xr))
    stop("length(Yr) must be equal to nrow(Xr)")
  
  if(any(is.na(Yr)))
    stop("The current version of the mbl function does not handle NAs in the response variable of the reference observations (Yr)")
  
  Xr <- as.matrix(Xr)
  Xu <- as.matrix(Xu)
  Yr <- as.matrix(Yr)
  
  rownames(Xr) <- 1:nrow(Xr)
  rownames(Xu) <- 1:nrow(Xu)
  
  if(is.null(colnames(Xr)))
    colnames(Xr) <- 1:ncol(Xr)
  
  if(is.null(colnames(Xu)))
    colnames(Xu) <- 1:ncol(Xu)
  
  if(sum(!colnames(Xu) ==  colnames(Xr)) != 0)
    stop("The names of the variables in Xr do not match the names of the variables in Xu")
  
  if(mblCtrl$sm=="none" & missing(dissimilarityM))
    stop("mblCtrl$sm is set to 'none' while 'dissimilarityM' is missing. Either similarity/disimilarity metric must be specified in mblCtrl$sm or a proper similarity/disimilarity matrix must be specified in the dissimilarityM argument")
  
  if(!missing(dissimilarityM)){
    if(!is.matrix(dissimilarityM))
      stop("'dissimilarityM' must be a matrix")
    if(mblCtrl$sm != "none"){
      warning(paste("Both 'dissimilarityM' and 'sm' ('mblCtrl$sm = ", mblCtrl$sm,"') were specified, only the 'dissimilarityM' argument will be taken into account and mblCtrl$sm will be set to 'none'", sep = ""))
      mblCtrl$sm <- "none"
    } 
    if(dissUsage == "predictors")
      if(sum(dim(dissimilarityM) - (nrow(Xr) + nrow(Xu))) != 0 & sum(diag(dissimilarityM) == 0) != 0)
        stop("When dissUsage = 'predictors', 'dissimilarityM' must be a square symmetric matrix of dissimilarities (derived from rbind(Xr, Xu)) for which the diagonal values are zeros")
    if(dissUsage %in% c("weights", "none"))
      if(nrow(dissimilarityM) != nrow(Xr) & ncol(dissimilarityM) != nrow(Xu))
        stop("When the 'dissUsage' argument is set to either 'weights' or 'none', 'dissimilarityM' must be a matrix representing the dissimilarity of each element in 'Xu' to each element in 'Xr'. The number of columns in 'dissimilarityM' must be equal to the number of rows of 'Xu' and the number of rows equal to the number of rows of 'Xr'")
  } else{   
    if(mblCtrl$sm == "movcor"){
      if(mblCtrl$ws < 3 | mblCtrl$ws > (ncol(Xr) - 1) | length(mblCtrl$ws) != 1) 
        mblCtrl$ws <- round(ncol(Xr) * 0.10)
      if(!mblCtrl$ws %% 2){
        mblCtrl$ws <- mblCtrl$ws - 1
        warning(paste("In this case the ws specified in mblCtrl$ws must be an unique odd value between 3 and ", (ncol(Xr) - 1), ". Therefore the window size was reset to ", mblCtrl$ws,"."))
      }
    }
  }
  
  if(!is.null(group))
  {
    if(length(group) != nrow(Xr))
      stop("The length of 'group' must be equal to the number of observations in 'Xr'")
  }
  
  pcSel <- match.arg(mblCtrl$pcSelection[[1]], c("opc", "var", "cumvar", "manual"))
  
  trsh <- mblCtrl$pcSelection$value
  
  if(pcSel == "opc" & mblCtrl$pcSelection$value > min(nrow(Xr) + nrow(Xu), ncol(Xr))){
    warning(paste("If 'mblCtrl$pcSelection$method = 'opc'' the value specified in 'mblCtrl$pcSelection$value' cannot be greater than  min(nrow(Xr) + nrow(Xu), ncol(Xr)) (i.e ", min(nrow(Xr) + nrow(Xu), ncol(Xr)),") in this case. Therefore the value was reset to ",min(nrow(Xr) + nrow(Xu), ncol(Xr)), sep=""))
    trsh <- min(nrow(Xr) + nrow(Xr), ncol(Xr))
  }
  
  if(pcSel == "manual" & mblCtrl$pcSelection$value > min(nrow(Xr) + nrow(Xu), ncol(Xr))){
    warning(paste("If 'mblCtrl$pcSelection$method = 'manual'' the value specified in 'mblCtrl$pcSelection$value' cannot be greater than  min(nrow(Xr) + nrow(Xu), ncol(Xr)) (i.e ", min(nrow(Xr) + nrow(Xu), ncol(Xr)),") in this case. Therefore the value was reset to ",min(nrow(Xr) + nrow(Xu), ncol(Xr)), sep=""))
    trsh <- min(nrow(Xr) + nrow(Xr), ncol(Xr))
  }
  
  match.arg(dissUsage, c("predictors", "weights", "none"))
  
  if(missing(k) & missing(k.diss)) 
    stop("Either k or k.diss must be specified")
  
  if(!missing(k) & !missing(k.diss)) 
    stop("Only one of k or k.diss can be specified")  
  
  if(mblCtrl$sm == "loc.pc")
    if(!missing(k))
      if(max(k) >= mblCtrl$k0)
        stop(paste("The maximum number of neighbours specified in the k vector (", max(k) ,") cannot be greather than the initial number of neighbours specified in 'mblCtrl$k0' (",mblCtrl$k0,")", sep = ""))
  
  if(!missing(k.diss)){
    dtc <- k.diss
    if(missing(k.range))
      stop("If the k.diss argument is used, k.range must be specified")
    if(length(k.range) != 2 | !is.numeric(k.range) | diff(k.range) < 0)
      stop("k.range must be a vector (of length 2) which specifies the minimum (first value) and the maximum (second value) number of neighbours") 
    k.min <- as.integer(k.range[1])
    k.max <- as.integer(k.range[2])
    if(k.min < 10)
      stop("Minimum number of nearest neighbours allowed is 10")
    if(k.max > nrow(Xr))
      stop("Maximum number of nearest neighbours cannot exceed the number of reference observations")
    if(mblCtrl$sm == "loc.pc")
      if(mblCtrl$k0 < max(k.range))
        stop(paste("The parameter 'k0' specified in 'mblCtrl$k0' cannot be lower than the maximum specified in the 'k.range' argument (", k.max,")"))
  }else
    dtc <- NULL
  
  ## Remove this line for the package version following 1.2
  if(method == "wapls2")
    stop("The method 'wapls2' is deprecated. Since it is no longer supported, choose instead one of the following options: 'pls', 'wapls1', 'gpr'")
  
  match.arg(method, c("pls", "wapls1", "gpr"))
  
  if(method %in% c("pls", "wapls1")){
    if(method %in% c("wapls1", "wapls2")){
      if(missing(pls.c))
        stop("When 'wapls1'is chosen, 'pls.c' must be specified")
      if(length(pls.c) != 2 | !is.numeric(pls.c) | missing(pls.c))
        stop("When either 'wapls1' or 'wapls2' are chosen, 'pls.c' must be a numerical vector of length 2 which specifies both the minimum (first value) and the maximum (second value) number of PLS components")
      if(diff(pls.c) < 0)
        stop("The number of minimum PLS components specified exedes the maximum. The first value in 'pls.c' must be the minimum number of PLS components")
    }else{
      if(missing(pls.c))
        stop("When 'pls' is chosen, 'pls.c' must be specified")
      if(length(pls.c) != 1 | !is.numeric(pls.c))
        stop("When 'pls'' is chosen, 'pls.c' must be a single numerical value specifiying the maximum number of PLS components to be evaluated")
    }
    plsF <- pls.c
  }else{
    plsF <- NULL
  }
  
  if(sum(c("loc_crossval", "NNv") %in% mblCtrl$valMethod) == 2)
    mblCtrl$valMethod <- "both"
  
  if("NNv" %in% mblCtrl$valMethod & nrow(Xu) < 3)
    stop("For nearest neighbour validation ('NNv') Xu must contain at least 3 observations")
  
  if(mblCtrl$valMethod %in% c("NNv", "both") & nrow(Xu) < 3)
    stop("For nearest neighbour validation (specified in mblCtrl$valMethod with 'NNv') Xu must contain at least 3 observations")
  
  if(!is.null(Yu)){
    Yu <- as.matrix(Yu)
    if(length(Yu) != nrow(Xu))
      stop("length(Yu) must be equal to nrow(Xu)")
    val <- data.frame(y = Yu) 
  } else {
    val <- data.frame(y = rep(NA, nrow(Xu)))
  }
  val$mat <- Xu
  
  if(!missing(k)){
    k <- as.integer(k)
    if(min(k) < 10)
      stop("Minimum number of nearest neighbours allowed is 10")
    if(max(k) > nrow(Xr))
      stop("The number of nearest neighbours cannot exceed the number of reference observations")
  }else
    k <- NULL
  
  call.f <-(match.call())    
  components <- NULL
  
  if(!missing(dissimilarityM)){
    if(dissUsage == "predictors"){
      d.cal.mat <- dissimilarityM[1:nrow(Xr), 1:nrow(Xr)]
      d.mat <- dissimilarityM[1:nrow(Xr),(1+nrow(Xr)):ncol(dissimilarityM)]
    }else {
      d.mat <- dissimilarityM
    }
    
  }else{
    # dissimilarity matrix computation between Xr and Xu     
    if(dissUsage != "predictors"){
      if(mblCtrl$sm == "euclid"){  
        f <- ncol(Xr)
        # fast dissimilarity function implemented in C++
        d.mat <- fDiss(Xr = Xr, X2 = Xu,  method = "euclid", center = mblCtrl$center, scaled = mblCtrl$scaled) # standardized distance
      } 
      
      if(mblCtrl$sm == "pc"){
        mtd <- ifelse(mblCtrl$pcMethod == "svd", "pca", "pca.nipals")  
        pcDistance <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, Yr = Yr, method = mtd, local = FALSE, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = FALSE, call. = FALSE, cores = mblCtrl$cores)
        d.mat <- pcDistance$dissimilarity[ , , drop = TRUE]
        npcs <- pcDistance$n.components
        rm(pcDistance)
        components <- list(n.components = npcs)
      } 
      
      if(mblCtrl$sm == "pls"){  
        pcDistance <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, method = "pls", Yr = Yr, local = FALSE, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = FALSE, call. = FALSE, cores = mblCtrl$cores)
        d.mat <- pcDistance$dissimilarity[ , , drop = TRUE]
        npcs <- pcDistance$n.components
        rm(pcDistance)
        components <- list(n.components = npcs)
      }       
      
      if(mblCtrl$sm == "loc.pc"){
        mtd <- ifelse(mblCtrl$pcMethod == "svd", "pca", "pca.nipals")  
        pcDissim <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, Yr = Yr, method = mtd, local = TRUE, k0 = mblCtrl$k0, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = FALSE, call.=FALSE)
        d.mat <- pcDissim$dissimilarity[ , , drop = TRUE]
        npcs <- pcDissim$n.components
        
        loc.r <- pcDissim$loc.n.components
        loc.r$sample.nm <- sub("X2", "Xu", loc.r$sample.nm)
        loc.npcs <- loc.r$loc.n.components
        
        rm(pcDissim)
        components <- list(n.components = npcs, loc.n.components = loc.npcs)
      } 
      
      if(mblCtrl$sm == "loc.pls"){
        pcDissim <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, Yr = Yr, method = "pls", local = TRUE, k0 = mblCtrl$k0, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = FALSE, call.=FALSE)
        d.mat <- pcDissim$dissimilarity[ , , drop = TRUE]
        
        loc.r <- pcDissim$loc.n.components
        loc.r$sample.nm <- sub("X2", "Xu", loc.r$sample.nm)
        loc.npcs <- loc.r$loc.n.components
        
        npcs <- pcDissim$n.components
        rm(pcDissim)
        components <- list(n.components = npcs, loc.n.components = loc.npcs)
      } 
      
      if(mblCtrl$sm == "movcor") {
        # correlation moving window approach as implemented in movcorDist.cpp   
        d.mat <-  corDiss(Xr = Xr, X2= Xu, ws = mblCtrl$ws, center = mblCtrl$center, scaled = mblCtrl$scaled)
      }
      
      if(mblCtrl$sm == "cor"){
        d.mat <-  corDiss(Xr = Xr, X2= Xu, ws = NULL, center = mblCtrl$center, scaled = mblCtrl$scaled)
      }      
      
      if(mblCtrl$sm == "sidF") {
        if(mblCtrl$center)
          warning("Centering was not applied in the computation of the sidF, because when the sid function is computed using mode = 'feature', it does not accept negative values (which are generated by centering).")
        d.mat <- sid(Xr = Xr, X2 = Xu, mode = "feature", center = FALSE, scaled = mblCtrl$scaled, reg = 10^-4)$sid
      }
      
      if(mblCtrl$sm == "sidD") {
        d.mat <- sid(Xr = Xr, X2 = Xu, mode = "density", center = mblCtrl$center, scaled = mblCtrl$scaled, reg = 10^-4)$sid
      }
      
      if(mblCtrl$sm == "cosine") {
        d.mat <-  fDiss(Xr = Xr, X2 = Xu, center = mblCtrl$center, scaled = mblCtrl$scaled, method = "cosine") 
      }
    } else{     
      # Distance matrix computation between all the objects resulting from rbind(Xr,Xu)  
      if(mblCtrl$sm == "euclid"){ 
        f <- ncol(Xr)
        # fast distance function implemented in C++
        d.mat <-   fDiss(Xr = rbind(Xr, Xu), X2 = NULL, center = mblCtrl$center, scaled = mblCtrl$scaled, method= "euclid") # standardized distance
        colnames(d.mat) <- paste("X2", 1:ncol(d.mat), sep = ".")
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
      } 
      
      if(mblCtrl$sm == "pc"){  
        mtd <- ifelse(mblCtrl$pcMethod == "svd", "pca", "pca.nipals")  
        pcDistance <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, Yr = Yr, method = mtd, local = FALSE, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = TRUE, call.=FALSE, cores = mblCtrl$cores)
        d.mat <- pcDistance$dissimilarity[ , , drop = TRUE]
        npcs <- pcDistance$n.components
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
        rm(pcDistance)
        components <- list(n.components = npcs)
      } 
      
      if(mblCtrl$sm == "pls"){  
        pcDistance <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, Yr = Yr, method = "pls", local = FALSE, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = TRUE, call.=FALSE, cores = mblCtrl$cores)
        d.mat <- pcDistance$dissimilarity[ , , drop = TRUE]
        npcs <- pcDistance$n.components
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
        rm(pcDistance)
        components <- list(n.components = npcs)
      } 
      
      if(mblCtrl$sm == "loc.pc"){
        mtd <- ifelse(mblCtrl$pcMethod == "svd", "pca", "pca.nipals") 
        pcDissim <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, Yr = Yr, method = mtd, local = TRUE, k0 = mblCtrl$k0, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = FALSE, call.=FALSE)
        d.mat <- pcDissim$dissimilarity[ , , drop = TRUE]
        npcs <- pcDissim$n.components
        loc.r <- pcDissim$loc.n.components
        loc.r$sample.nm <- sub("X2", "Xu", loc.r$sample.nm)
        loc.npcs <- loc.r$loc.n.components
        rm(pcDissim)
        components <- list(n.components = npcs, loc.n.components = loc.npcs)
      } 
      
      if(mblCtrl$sm == "loc.pls"){
        pcDissim <- orthoDiss(Xr = Xr, X2 = Xu, pcSelection = mblCtrl$pcSelection, Yr = Yr, method = "pls", local = TRUE, k0 = mblCtrl$k0, center = mblCtrl$center, scaled = mblCtrl$scaled, return.all = FALSE, call.=FALSE)
        d.mat <- pcDissim$dissimilarity[ , , drop = TRUE]
        npcs <- pcDissim$n.components    
        loc.r <- pcDissim$loc.n.components
        loc.r$sample.nm <- sub("X2", "Xu", loc.r$sample.nm)
        loc.npcs <- loc.r$loc.n.components
        rm(pcDissim)
        components <- list(n.components = npcs, loc.n.components = loc.npcs)
      } 
      
      if(mblCtrl$sm == "movcor") {
        # correlation moving window approach as implemented in movcorDist.cpp
        d.mat <- corDiss(Xr = rbind(Xr, Xu), X2 = NULL, ws = mblCtrl$ws, center = mblCtrl$center, scaled = mblCtrl$scaled)
        colnames(d.mat) <- paste("X2", 1:ncol(d.mat), sep = ".")
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
      }
      
      if(mblCtrl$sm == "cor") {
        d.mat <-  corDiss(Xr = rbind(Xr, Xu), X2 = NULL, ws = NULL, center = mblCtrl$center, scaled = mblCtrl$scaled)
        colnames(d.mat) <- paste("X2", 1:ncol(d.mat), sep = ".")
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
      }
      
      if(mblCtrl$sm == "sidF") {
        if(mblCtrl$center)
          warning("Centering was not applied in the computation of the sidF, because when the sid function is computed using mode = 'feature', it does not accept negative values (which are generated by centering).")
        d.mat <- sid(Xr = rbind(Xr, Xu), X2 = NULL, mode = "feature", center = FALSE, scaled = mblCtrl$scaled, ...)$sid
        colnames(d.mat) <- paste("X2", 1:ncol(d.mat), sep = ".")
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
      }
      
      if(mblCtrl$sm == "sidD") {
        d.mat <- sid(Xr = rbind(Xr, Xu), X2 = NULL, mode = "density", center = mblCtrl$center, scaled = mblCtrl$scaled, ...)$sid
        colnames(d.mat) <- paste("X2", 1:ncol(d.mat), sep = ".")
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
      }
      
      if(mblCtrl$sm == "cosine") {
        d.mat <-   fDiss(Xr = rbind(Xr, Xu), X2 = NULL, center = mblCtrl$center, scaled = mblCtrl$scaled, method = "cosine") 
        colnames(d.mat) <- paste("X2", 1:ncol(d.mat), sep = ".")
        d.cal.mat <- d.mat[1:nrow(Xr), 1:nrow(Xr)]
        d.mat <- d.mat[1:nrow(Xr),(1+nrow(Xr)):ncol(d.mat)]
      }
    }
  }
  dimnames(d.mat) <- NULL
  
  if(dissUsage!="predictors")
    d.cal.mat <- NULL
  
  # Global H
  mtd <- ifelse(mblCtrl$pcMethod == "svd", "pca", "pca.nipals")  
  pca <- pcProjection(Xr = rbind(Xr, Xu), X2 = NULL, pcSelection = mblCtrl$pcSelection, method = mtd, Yr = c(Yr, rep(NA, nrow(Xu))), center = mblCtrl$center, scaled = mblCtrl$scaled, call.=FALSE, cores = mblCtrl$cores)
  npcs <- pca$n.components 
  
  scores <- sweep(pca$scores, 2, pca$sc.sdv, "/")
  scores.cal <- scores[1:nrow(Xr), , drop=F]
  scores.val <- scores[(1+nrow(Xr)):nrow(scores), , drop=F]
  
  GH <- as.vector((fDiss(Xr = scores.val, X2 = t(colMeans(scores.cal)), center = FALSE, scaled = FALSE, method = "euclid")))
  rm(scores)
  
  # convert k.dis into k
  if(is.null(k)) {
    dtc <- sort(dtc)    
    k <- matrix(NA, nrow = nrow(Xu), ncol = length(dtc))
    for(j in 1:length(dtc)){
      k[,j] <- apply(d.mat,2,function(x) sum(x <= dtc[j], na.rm = TRUE))
    }      
    k.org <- k
    k_mat <- pmin(pmax(k.org, k.min), k.max) 
  } else {
    k.org <- k_mat <- matrix(k, nrow = nrow(Xu), ncol = length(k), byrow = TRUE)    
  }  
  
  
  if(mblCtrl$valMethod %in% c("loc_crossval", "both") & method %in% c("pls", "wapls1"))
  {
    if(((floor(min(min(k_mat), ncol(Xr)) * mblCtrl$p)) - 1) < min(plsF)){
      stop(paste("The number of pls components must be lower than ", mblCtrl$p*100, "% (i.e. mblCtrl$p) of the number of neighbours to be used"))
    }
  }
  
  if(method != "gpr"){
    if(floor(min(min(k_mat), ncol(Xr)) * mblCtrl$p) < min(plsF)){
      stop("The number of pls components must be lower than the number of neighbours to be used")
    }
  }
  
  
  
  it <- i <- kn.org <- NULL
  predobs <- foreach(tmp.val = iter(val, by = "row"), kn.org = iter(k.org, by = "row"), i = icount(),
                     it = isubset(x=Xr,Dx=d.cal.mat,D=d.mat,k=k_mat),
                     .combine = rbind, .inorder = FALSE,
                     .export = c("orthoDiss",  
                                 "orthoProjection", 
                                 "simEval", "corDiss", "fDiss", 
                                 "locFitnpred", "plsCv", "cSds", "wapls.weights")) %mydo%{            
                                   
                                   #    if(mblCtrl$valMethod %in% c("loc_crossval", "both") & method %in% c("pls", "wapls1"))
                                   #    {
                                   #      if(((floor(min(it$k, nrow(it$d)) * mblCtrl$p)) - 1) < min(plsF)){
                                   #        stop(paste("The number of pls components must be lower than ", mblCtrl$p*100, "% (i.e. mblCtrl$p) of the number of neighbours to be used"))
                                   #      }
                                   #    }
                                   #    
                                   #    if(method != "gpr"){
                                   #      if(((floor(min(it$k, nrow(it$d)))) - 1) < min(plsF)){
                                   #        stop("The number of pls components must be lower than the number of neighbours to be used")
                                   #      }
                                   #    }
                                   #    
                                   # k nearest neighbours loop. Must be within the sample loop to allow possibly for 'local' tuning
                                   
                                   predobs <- data.frame(o.index = i, 
                                                         k.org = c(kn.org), k = it$k,
                                                         pls.comp = NA,
                                                         min.pls = NA,
                                                         max.pls = NA,
                                                         yu.obs = tmp.val$y, pred = NA, 
                                                         yr.min.obs = NA, yr.max.obs = NA,
                                                         index.nearest.in.ref = it$index[1],
                                                         y.nearest = Yr[it$index[1]],
                                                         y.nearest.pred = NA,
                                                         loc.rmse.cv = NA,
                                                         loc.st.rmse.cv = NA,
                                                         dist.nearest = min(it$d),
                                                         dist.k.farthest = NA, rep = NA)
                                   
                                   
                                   if(mblCtrl$sm == "loc.pc" & dissUsage == "predictors")
                                   {
                                     pca.n <- plsProjection(Xr = it$x, Yr = Yr[it$index[1:max(it$k)]], Xu = tmp.val$mat, nPf = loc.npcs[i], scaled = mblCtrl$scaled)
                                     pca.i <- list(x = pca.n$scores, sdev = cSds(pca.n$scores))
                                     scores.i <- sweep(pca.i$x[,1:loc.npcs[i]], 2, pca.i$sdev[1:loc.npcs[i]], "/")
                                     locD.i <-   fDiss(Xr = scores.i, X2 = scores.i, center = FALSE, scaled = FALSE, method = "euclid")     
                                   }
                                   
                                   if(mblCtrl$progress)
                                   {
                                     cat ((paste("Predicting sample:",i)) ," ") 
                                     pb <- txtProgressBar(width = 10, char = "-")
                                   }
                                   
                                   for(kk in 1:length(it$k))
                                   {
                                     knkk <- it$k[kk]
                                     tmp.val.k <- tmp.val
                                     
                                     if(mblCtrl$progress)
                                       setTxtProgressBar(pb, kk/length(it$k))
                                     
                                     # If the sample has not been predicted before, then create a model and predict it
                                     if(knkk != ifelse(kk == 1, 0, it$k[kk-1])){
                                       
                                       predobs$rep[kk] <- 0
                                       
                                       ###############################################################################################
                                       dk <- it$index[1:knkk]
                                       
                                       # Select cal
                                       tmp.cal <- data.frame(y = Yr[dk])
                                       
                                       if(!is.null(group))
                                         tmp.group <- as.factor(as.character(group[dk]))
                                       else
                                         tmp.group <- NULL
                                       
                                       tmp.cal$mat  <- it$x[1:knkk,,drop=F]                           
                                       
                                       minmax <- range(tmp.cal$y, na.rm=T) 
                                       predobs$yr.min.obs[kk] <- minmax[1]
                                       predobs$yr.max.obs[kk] <- minmax[2]
                                       predobs$dist.k.farthest[kk] <- max(it$d[dk])
                                       
                                       if(dissUsage == "weights")
                                       {
                                         # Weights are defined according to a tricubic function 
                                         # as in Cleveland and Devlin (1988) and Naes and Isaksson (1990).
                                         stdd <- it$d[dk]/max(it$d[dk])
                                         i.wgts <- (1 - (stdd^3))^3
                                         i.wgts[which(i.wgts == 0)] <- 1e-04
                                       } else {i.wgts <- rep(1,knkk)}
                                       
                                       if(dissUsage == "predictors")
                                       {
                                         if(kk == 1) {tmp.val.k$mat <- tmp.val$mat}
                                         tmp.val$mat <- tmp.val.k$mat
                                         nms <- c(paste("k",(1:knkk), sep = ""), colnames(tmp.cal$mat))
                                         if(mblCtrl$sm == "loc.pc")
                                         {           
                                           tmp.cal$mat <- cbind(locD.i[1:knkk, 1:knkk] ,tmp.cal$mat)
                                         }else{
                                           tmp.cal$mat <- cbind(it$dx[1:knkk, 1:knkk] ,tmp.cal$mat)
                                         }
                                         tmp.val.k$mat <- cbind(t(it$d[dk]), tmp.val$mat)
                                         colnames(tmp.cal$mat) <- nms
                                         colnames(tmp.val.k$mat) <- nms
                                       }
                                       
                                       scale <- mblCtrl$scaled
                                       
                                       # local fit
                                       i.pred <- locFitnpred(x = tmp.cal$mat, y = tmp.cal$y, 
                                                             predMethod = method, 
                                                             scaled = scale, pls.c = plsF,
                                                             weights = i.wgts,
                                                             newdata = as.vector(tmp.val.k$mat), 
                                                             CV = mblCtrl$valMethod %in% c("loc_crossval", "both"),
                                                             .optimize = mblCtrl$localOptimization,
                                                             group = tmp.group,
                                                             p = mblCtrl$p, resampling = mblCtrl$resampling,
                                                             noise.v = noise.v, range.pred.lim = mblCtrl$range.pred.lim,
                                                             pls.max.iter = pls.max.iter, pls.tol = pls.tol)
                                        
                                        predobs$pred[kk] <- i.pred$prediction
                                       
                                        
                                        if(mblCtrl$valMethod %in% c("loc_crossval","both")){
                                          if(mblCtrl$localOptimization & method %in% c("pls", "wapls1")){
                                            if(method == "pls"){
                                              o <- i.pred$validation$bestpls.c
                                              predobs$pls.comp[kk] <- o
                                              oplsfs <- o
                                            }else{
                                              o <- which.min(i.pred$validation$cvResults$rmse.cv)
                                              oplsfs <- i.pred$validation$cvResults[o,c("minF", "maxF")]
                                              predobs$min.pls[kk] <- oplsfs[[1]]
                                              predobs$max.pls[kk] <- oplsfs[[2]]
                                              
                                            }
                                            predobs$loc.rmse.cv[kk] <- i.pred$validation$cvResults$rmse.cv[o]
                                            predobs$loc.st.rmse.cv[kk] <- i.pred$validation$cvResults$st.rmse.cv[o]
                                          }else{
                                            oplsfs <- plsF
                                            o <- ifelse(method == "pls", plsF, 1)
                                            predobs$pls.comp[kk] <- ifelse(method == "pls", plsF, NA)
                                            predobs$min.pls[kk] <- ifelse(method == "wapls1", plsF[[1]], NA)
                                            predobs$max.pls[kk] <- ifelse(method == "wapls1", plsF[[2]], NA)
                                            predobs$loc.rmse.cv[kk] <- i.pred$validation$cvResults$rmse.cv[o]
                                            predobs$loc.st.rmse.cv[kk] <- i.pred$validation$cvResults$st.rmse.cv[o]
                                          }
                                        }else{
                                          predobs$pls.comp[kk] <- ifelse(method == "pls", plsF, NA)
                                          predobs$min.pls[kk] <- ifelse(method == "wapls1", plsF[[1]], NA)
                                          predobs$max.pls[kk] <- ifelse(method == "wapls1", plsF[[2]], NA)
                                          oplsfs <- plsF
                                        }
                                       
                                       
                                       if(mblCtrl$valMethod %in% c("NNv","both"))
                                       {
                                         if(!is.null(group))
                                           out.g <- which(tmp.group == tmp.group[[1]])
                                         else
                                           out.g <- 1
                                         
                                         nearest.pred <- locFitnpred(x = tmp.cal$mat[-c(out.g),], y = tmp.cal$y[-c(out.g)], 
                                                                     predMethod = method, 
                                                                     scaled = scale, pls.c = oplsfs, 
                                                                     noise.v = noise.v,
                                                                     newdata = as.vector(tmp.cal$mat[1,]), 
                                                                     CV = FALSE,  .optimize = FALSE, range.pred.lim = mblCtrl$range.pred.lim, 
                                                                     pls.max.iter = pls.max.iter, pls.tol = pls.tol)$prediction
                                         
                                         predobs$y.nearest.pred[kk] <- nearest.pred/i.wgts[1]
                                       }
                                     } else {
                                       predobs[kk,] <- predobs[(kk-1),]
                                       predobs$rep[[kk]] <- 1
                                     }     
                                   }
                                   if(mblCtrl$progress)
                                     close(pb)
                                   
                                   predobs$k <- as.factor(predobs$k)
                                   predobs
                                 }
  
  out <-c(if(missing(Yu)){"yu.obs"},
          if(is.null(dtc)){"k.org"},
          if(!(mblCtrl$valMethod %in% c("NNv", "both"))){"y.nearest.pred"},
          if(method != "wapls1"){c("min.pls","max.pls")},
          if(method != "pls"){"pls.comp"},
          if(!(mblCtrl$valMethod %in% c("loc_crossval", "both"))){c("loc.rmse.cv", "loc.st.rmse.cv")},
          "rep")
  
  predobs <- predobs[,!(colnames(predobs) %in% out)]
  predobs <- predobs[order(predobs$o.index),]
  
  dfrn <- function(x,...){
    x <- data.frame(x)
    rownames(x) <- 1:nrow(x)
    return(x)
  }
  
  if(is.null(dtc)){  
    predobs <- by(predobs, as.numeric(paste(predobs$k, ".",predobs$rep, sep="")), dfrn)
    names(predobs) <- paste("Nearest_neighbours_", k, sep = "")    
    trh <- data.frame(k = k)
  } else{
    d.dist <- rep(dtc, nrow(Xu))
    distance <- (predobs$k.org - as.numeric(as.vector(predobs$k))) == 0
    nms <- colnames(predobs)
    predobs <- cbind(predobs, k.diss = d.dist, distance = distance)
    nms <- c(nms[1], "k.diss", "distance" ,nms[-1])
    predobs <- predobs[,nms]
    predobs <- by(predobs, predobs$k.diss, FUN = dfrn)
    names(predobs) <- paste("k distance:", dtc)
    bounded <- function(x, ...){
      return(sum(x$distance == FALSE))
    }
    trh <- data.frame(k.diss = dtc, p.bounded = paste(round(100*sapply(predobs, bounded)/nrow(Xu), 3),"%", sep =""))
  }
  
  if(mblCtrl$valMethod %in% c("loc_crossval", "both"))
  {
    mean.loc.res <- function(x,...){
      mean.loc.rmse <- mean(x$loc.rmse.cv)
      mean.loc.st.rmse <- mean(x$loc.st.rmse.cv)
      return(c(loc.rmse = mean.loc.rmse, loc.st.rmse = mean.loc.st.rmse))
    }
    loc.res <- as.data.frame(t(sapply(predobs, mean.loc.res)))
    loc.res <- cbind(trh, 
                     rmse = loc.res$loc.rmse, 
                     st.rmse = loc.res$loc.st.rmse)
  } else {loc.res <- NULL}
  
  if(mblCtrl$valMethod %in% c("NNv", "both"))
  {
    nn.stats <- function(x,...){
      nn.rmse <- (mean((x$y.nearest - x$y.nearest.pred)^2))^0.5
      nn.st.rmse <- nn.rmse / diff(range(x$y.nearest))  
      nn.rsq <- (cor(x$y.nearest, x$y.nearest.pred))^2
      return(c(nn.rmse = nn.rmse, nn.st.rmse = nn.st.rmse, nn.rsq = nn.rsq))
    }
    
    loc.nn.res <- as.data.frame(t(sapply(predobs, nn.stats)))
    loc.nn.res <- cbind(trh, 
                        rmse = loc.nn.res$nn.rmse, 
                        st.rmse = loc.nn.res$nn.st.rmse, 
                        r2 = loc.nn.res$nn.rsq)
  } else {loc.nn.res <- NULL}
  
  if(!is.null(Yu))
  {
    yu.stats <- function(x,...){
      yu.rmse <- (mean((x$yu.obs - x$pred)^2))^0.5
      yu.st.rmse <- yu.rmse / diff(range(x$yu.obs)) 
      yu.rsq <- (cor(x$yu.obs, x$pred))^2
      return(c(yu.rmse = yu.rmse, yu.st.rmse = yu.st.rmse, yu.rsq = yu.rsq))
    }
    pred.res <- as.data.frame(t(sapply(predobs, yu.stats)))
    pred.res <- cbind(trh, 
                      rmse = pred.res$yu.rmse, 
                      st.rmse = pred.res$yu.st.rmse, 
                      r2 = pred.res$yu.rsq)
  } else {pred.res <- NULL}
  
  if(mblCtrl$valMethod == "none"){loc.res <- NULL}
  
  mblCtrl$valMethod <- ini.cntrl$valMethod  
  
  if(!is.null(mblCtrl$sm)){
    if(mblCtrl$sm %in% c("movcor", "loc.pc", "loc.pls")){
      cntrlParam <- list(sm = mblCtrl$sm,
                         smParam = ifelse(mblCtrl$sm == "movcor", mblCtrl$ws, mblCtrl$k0),
                         returnDiss = mblCtrl$returnDiss,
                         pcSelection = list(method = pcSel, value = trsh),
                         pcMethod = mblCtrl$pcMethod,
                         center = mblCtrl$center,
                         scaled = mblCtrl$scaled,
                         valMethod = mblCtrl$valMethod,
                         resampling = mblCtrl$resampling, 
                         p = mblCtrl$p,
                         range.pred.lim = mblCtrl$range.pred.lim,
                         progress = mblCtrl$progress, 
                         cores = 1)
      nm <- ifelse(mblCtrl$sm == "movcor", "ws", "k0")
      names(cntrlParam)[names(cntrlParam) == "smParam"] <- nm
    } else{
      cntrlParam <- list(sm = mblCtrl$sm,
                         returnDiss = mblCtrl$returnDiss,
                         pcSelection = list(method = pcSel, value = trsh),
                         pcMethod = mblCtrl$pcMethod,
                         center = mblCtrl$center,
                         scaled = mblCtrl$scaled,
                         valMethod = mblCtrl$valMethod,
                         localOptimization = mblCtrl$localOptimization,
                         resampling = mblCtrl$resampling, 
                         p = mblCtrl$p,
                         range.pred.lim = mblCtrl$range.pred.lim,
                         progress = mblCtrl$progress, 
                         cores = mblCtrl$cores,
                         allowParallel = mblCtrl$cores)
    }
  }
  
  if(!sum(unlist(ini.cntrl) == unlist(cntrlParam))){
    cntrlParam <- list(initial = ini.cntrl, finallyUsed = cntrlParam)
  }
  
  if(ini.cntrl$returnDiss){
    s.meth <- ifelse(mblCtrl$sm == "none", "A matrix provided by the user through the 'dissimilarityM' argument", mblCtrl$sm)
    if(dissUsage == "predictors"){
      dissimilarities <- list(method = s.meth, Xr_Xu = d.mat, Xr_Xr = d.cal.mat)
    }else{
      dissimilarities <- list(method = s.meth, Xr_Xu = d.mat) 
    }
  }
  
  resultsList <- list(call = call.f,
                      cntrlParam = cntrlParam,
                      dissimilarities = if(ini.cntrl$returnDiss){dissimilarities}else{NULL}, 
                      totalSamplesPredicted = nrow(Xu),
                      pcAnalysis = list(scores_Xr = scores.cal, scores_Xu = scores.val, n.componentsUsed = npcs, GH = GH),
                      components = components,
                      localCrossValStats = loc.res, 
                      nnValStats = loc.nn.res, 
                      YuPredictionStats = pred.res, 
                      results = predobs)
  
  resultsList <- resultsList[!sapply(resultsList,is.null)]
  
  attr(resultsList, "call") <- call.f
  class(resultsList) <- c("mbl","list")  
  return(resultsList)
}

#' @title Gaussian process regression with dot product covariance
#' @description internal
#' @keywords internal
gpr.dp <- function(Xn, Yn, noise.v = 0.001, scaled = TRUE, center = FALSE, x.add = NULL){
  #   if(!is.null(x.add)){
  #     Xn <- rbind(x.add, Xn)
  #   }
  if(center)
  {
    center <- colMeans(Xn)
    Xn <- sweep(x = Xn, MARGIN = 2, STATS = center, FUN = "-")
  }else{
    center <- NULL
  }
  if(sum(scaled))
  {
    scale <- cSds(Xn)
    if(length(scaled)> 1){
      scale <- scale * scaled
      scale[scale == 0] <- 1
    }
    Xn <- sweep(x = Xn, MARGIN = 2, STATS = scale, FUN = "/")
  }else{
    scale <- NULL
  }
  
  #   if(!is.null(x.add)){
  #     Xn <- Xn[-1,]
  #   }
  #   
  K <- Xn %*% t(Xn)
  vrnc <- diag(rep(noise.v, nrow(Xn)))
  alpha <- solve((K + vrnc), Yn)
  return(list(Xn = Xn, alpha = alpha, scale = scale, center = center))
}

#' @title Prediction function for "gpr.dp" (Gaussian process regression with dot product covariance)
#' @description internal
#' @keywords internal
pred.gpr.dp <- function(object, newdata){
  if(!is.null(object$center)){
    newdata <- sweep(x = newdata, MARGIN = 2, STATS = object$center, FUN = "-")
  }
  if(!is.null(object$scale)){
    newdata <- sweep(x = newdata, MARGIN = 2, STATS = object$scale, FUN = "/")
  }
  predicted <- newdata %*% t(object$Xn) %*% object$alpha
}

#' @title Local multivariate regression
#' @description internal
#' @keywords internal
locFitnpred <- function(x, y, predMethod, scaled = TRUE, weights = NULL, newdata, 
                        pls.c = NULL, CV = FALSE, .optimize = FALSE, resampling = 10, p = 0.75, group = NULL, noise.v = 0.001, 
                        range.pred.lim = TRUE,
                        pls.max.iter = 1, pls.tol = 1e-6){
  
  ## Sanity checks are dactivated here for speeding-up the code. 
  ## This function is only used within the mbl and the checks here are already done either in that function or
  ## in the mblControl function. Since the present function goes into a main loop in the mbl funtion,
  ## this means that the checks will be done during every iteration
  #   if(nrow(x) != length(y))
  #     stop("length of the response variable (y) does not match with the number of observations in x")
  
  #   if(!is.vector(weights) & !is.null(weights))
  #     stop("if the argument weights is specified it must be a vector")
  #   if(!is.null(weights)){
  #     if(length(weights) != length(y))
  #       stop("Length of the vector of weights does not match with the number of observations in x and y")
  #   }
  
  #   if(predMethod != "gpr"){
  #     if(predMethod == "pls" & length(pls.c) != 1)
  #       stop("When predMethod = 'pls', pls.c must be a vector of length 1")
  #     if(predMethod %in% c("wapls1"))
  #       if(length(pls.c) != 2)
  #         stop("When 'wapls1' is used as regression method, pls.c must be a vector of length 2 indicating a minimum and a maximum number of PLS components")
  #   }
  
  if(is.null(weights))
    weights <- 1
  
  #   if(CV){
  #     if(!is.numeric(resampling) | length(resampling) != 1)
  #       stop("resampling must be a single numeric value")
  #     if(!is.numeric(p) | length(p) != 1 | p >= 1 | p <= 0)
  #       stop("p must be a single numeric value higher than 0 and lower than 1")
  #   }
  
  
  #   if(!is.vector(newdata))
  #     stop("'newdata' must be a vector")
  
  #   if(length(newdata) != ncol(x))
  #     stop("Length of vector newdata must be equal to the number of columns of x")
  
  newdata <- t(newdata)
  
  if(any(cSds(x)==0)){
    warning("One of the variables has zero variance. Data will not be scaled")
    scaled <- FALSE
  }
  
  results <- cvVal <- pred <- NULL
  if(predMethod=="gpr") {    
    if(CV)
    {
      cvVal <- gprCv(x = x, y = as.matrix(y), scaled = scaled, weights = weights, p = p, resampling = resampling, 
                     group = group, retrieve = "final.model")
      fit <- cvVal$model
    } else { 
      x <- sweep(x, 1, weights, "*")   ###MODIFIED
      y <- y * weights
      ##fit <- gpr.dp(Xn = x, Yn = y, noise.v = noise.v, scaled = scaled, x.add = newdata)
      fit <- gprdp(X = x, Y = as.matrix(y), noisev = noise.v, scale = scaled)
    }
    
    ##pred <- pred.gpr.dp(fit, newdata = newdata)
    pred <- predgprdp(Xz = fit$Xz, 
                      alpha = fit$alpha, 
                      newdata = newdata,
                      scale = fit$is.scaled,
                      Xcenter = fit$Xcenter,
                      Xscale = fit$Xscale,
                      Ycenter = fit$Ycenter,
                      Yscale = fit$Yscale)[[1]]
  } 
  
  if(predMethod == "pls"){
    
    if(CV)
    {
      cvVal <- plsCv(x = x, y = as.matrix(y), ncomp = pls.c, 
                     method = "pls",
                     center = TRUE, scaled = scaled, 
                     weights = weights, 
                     p = p, resampling = resampling, 
                     group = group, 
                     retrieve = TRUE,
                     .optimize = .optimize,
                     max.iter = pls.max.iter, tol = pls.tol)
      
      
      
      fit <- cvVal$models
      ncomp <- cvVal$bestpls.c
      
    } else {
      x <- sweep(x, 1, weights, "*")   ###MODIFIED
      y <- y * weights
      fit <- fopls(X = x, Y = as.matrix(y), scale = scaled, ncomp = pls.c, maxiter = pls.max.iter, tol = pls.tol)
      
      ncomp <- pls.c
    }
    pred <- predopls(bo = fit$bo,
                     b = fit$coefficients, 
                     ncomp = ncomp, 
                     newdata = newdata,
                     scale = ifelse(nrow(fit$transf$Xscale) == 1, TRUE, FALSE),
                     Xscale = fit$transf$Xscale)[,ncomp]
  }
  if(predMethod == "wapls1"){
    if(CV)
    {
      cvVal <- plsCv(x = x, y = as.matrix(y), 
                     ncomp = pls.c[[2]],
                     method = "wapls1",
                     center = TRUE, scaled = scaled, 
                     minF = pls.c[[1]],
                     newX = newdata,
                     weights = weights, 
                     p = p, resampling = resampling, 
                     group = group, 
                     retrieve = TRUE,
                     .optimize = .optimize,
                     max.iter = pls.max.iter, tol = pls.tol)
      
      fit <- cvVal$models
      #rstls <- w[pls.c[[1]]:pls.c[[2]]] * colMeans(cvVal$cvResults[pls.c[[1]]:pls.c[[2]],])
      if(!.optimize){
        cvVal$cvResults <- data.frame(minF = pls.c[[1]], 
                                      maxF = pls.c[[2]], 
                                      rmse.cv = cvVal$cvResults$rmse.cv,
                                      st.rmse.cv = cvVal$cvResults$st.rmse.cv,
                                      rmse.sd.cv = cvVal$cvResults$rmse.sd.cv,     
                                      rsq.cv = cvVal$cvResults$rsq.cv)
        
      }
      
      
      ## here all the weights are output from 1 to plsMax
      # w <- cvVal$compweights[!cvVal$compweights == 0]
      w <- cvVal$compweights
    } else {
      x <- sweep(x, 1, weights, "*")   ###MODIFIED
      y <- y * weights
      fit <- fopls(X = x, Y = as.matrix(y), scale = scaled, ncomp = pls.c[[2]], maxiter = pls.max.iter, tol = pls.tol)
      
      # compute weights for PLS components selected (from plsMin to plsMax)
      w <- rep(0, pls.c[[2]])
      w[pls.c[[1]]:pls.c[[2]]] <- wapls.weights(plsO = fit, orgX = x, type = "w1", newX = t(newdata), pls.c = pls.c)       
    }
    
    
    # compute the weighted average of the multiple PLS predictions
    #pred <- sum(c(predict(fit, newdata, ncomp = pls.c[[1]]:pls.c[[2]])) * w) # weighted preds
    
    predctn <- predopls(bo = fit$bo,
                        b = fit$coefficients, 
                        ncomp = pls.c[[2]], 
                        newdata = newdata,
                        scale = ifelse(nrow(fit$transf$Xscale) == 1, TRUE, FALSE),
                        Xscale = fit$transf$Xscale)
    pred <- sum(c(predctn) * w) # weighted preds
  }
  
  if(range.pred.lim)
  {
    if(pred > max(y))
      pred <- max(y)
    if(pred < min(y))
      pred <- min(y)
  }
  
  results$validation <- cvVal
  results$prediction <- as.vector(pred)
  
  return(results)
}

#' @title A simple ramdom sampling function
#' @description internal
#' @keywords internal
smpl <- function(x){
  if(length(x) == 1)
    x <- c(x, x)
  return(sample(x, size = 1))
}


#' @title Cross validation for Gaussian process regression
#' @description internal
#' @keywords internal
gprCv <- function(x, y, scaled, weights = NULL, p = 0.75, resampling = 10, group = NULL, noise.v = 0.001, retrieve = c("final.model", "none")){
  
  ## Create the resampling groups
  if(is.null(group)){
    nv <- floor((1-p)*nrow(x))
    y.str <- quantile(y, probs = seq(0, 1, length = (nv+1)))
    y.cut <- cut(y, unique(y.str), include.lowest = TRUE)
    levels(y.cut) <- 1:length(y)
    prob.g <- data.frame(orig.order = 1:length(y), prob.group = as.numeric (y.cut))
    
    rspi <- matrix(NA, nv, resampling)
    mrspi <- matrix(NA, nrow(x) - nv, resampling)
    colnames(rspi) <- colnames(mrspi) <- paste("Resample", seq(1:resampling), sep = "")
    rownames(rspi) <- paste("index", seq(1:nv), sep = "")
    rownames(mrspi) <- paste("index", seq(1:nrow(mrspi)), sep = "")
    
    idxs <- 1:nrow(x)
    for(jj in 1:resampling)
    {
      strs <- tapply(X = prob.g$orig.order, FUN = smpl, INDEX = prob.g$prob.group)
      if(length(strs) < nv)
      {
        adds <- sample((1:length(y))[-strs], nv - length(strs))
        strs <- c(strs, adds)
      }
      rspi[,jj]<- strs
      mrspi[,jj]<- idxs[!idxs %in% strs]
    }
  }else{
    ygr <- data.frame(y, grp = as.factor(as.character(group)), order.orig = 1:length(y))
    lvls <- levels(ygr$grp)
    nlv <- nlevels(ygr$grp)
    # Compute means per group and make the stratified sampling based on that 
    # vector of means
    yg <- tapply(X = ygr$y, INDEX = ygr$grp, FUN = mean)
    
    nv <- floor((1 - p) * nlv)
    y.str <- quantile(yg, probs = seq(0, 1, length = (nv+1)))
    y.cut <- cut(yg, unique(y.str), include.lowest = TRUE)
    levels(y.cut) <- paste("l_", 1:nlevels(y.cut), sep = "")
    prob.g <- data.frame(orig.order = 1:length(yg), prob.group = factor(y.cut))
    
    rspi <- matrix(0, nv, resampling)
    colnames(rspi) <- paste("Resample", seq(1:resampling), sep = "")
    rownames(rspi) <- paste("index", seq(1:nv), sep = "")
    
    listrspi <- mlistrspi <- NULL
    for(jj in 1:resampling)
    {
      strs <- tapply(X = prob.g$orig.order, FUN = smpl, INDEX = prob.g$prob.group)
      if(length(strs) < nv)
      {
        adds <- sample((1:nrow(prob.g))[-strs], nv - length(strs))
        strs <- c(strs, adds)
      }
      
      strs <- lvls[strs]  
      
      listrspi[[jj]] <- ygr$order.orig[ygr$grp %in% strs]
      mlistrspi[[jj]] <- ygr$order.orig[!ygr$grp %in% strs]
    }
    
    lensam <- lengths(listrspi)
    mlensam <- lengths(mlistrspi)
    rspi <- matrix(NA, max(lensam), resampling)
    mrspi <- matrix(NA, max(mlensam), resampling)
    colnames(rspi) <- colnames(mrspi) <- paste("Resample", seq(1:resampling), sep = "")
    rownames(rspi) <- paste("index", 1:nrow(rspi), sep = "")
    rownames(mrspi) <- paste("index", 1:nrow(mrspi), sep = "")
    idxs <- 1:nrow(x)
    for(jj in 1:resampling)
    {
      ## ramdomly select the replacements for the missing indexes
      ## replacements are necessary!
      rspi[,jj] <- c(sample(listrspi[[jj]], max(lensam) - lensam[jj]), listrspi[[jj]])
      mrspi[,jj] <- c(sample(mlistrspi[[jj]], max(mlensam) - mlensam[jj]), mlistrspi[[jj]])
    }
  }
  
  cvre <- pgpcv_cpp(X = x, 
                    Y = y, 
                    mindices = mrspi,
                    pindices = rspi,
                    noisev = noise.v,
                    scale = scaled)
  
  val <- NULL
  val$resamples <- rspi
  val$cvResults <- data.frame(rmse.cv = mean(cvre$rmse.seg), st.rmse.cv = mean(cvre$st.rmse.seg), rmse.sd.cv = sd(cvre$rmse.seg), rsq.cv = mean(cvre$rsq.seg))
  
  if(retrieve == "final.model")
  {
    if(is.null(weights)){
      x <- sweep(x, 1, weights, "*") #### MODIFIED
      y <- y * weights
    }
    val$model <-   gprdp(X = x, Y = as.matrix(y), noisev = noise.v, scale = scaled)
  }
  return(val)
}

# #' @title Standard deviation of columns
# #' @description internal
# #' @keywords internal
# cSds <-  function(x){
#   return(apply(x,2,sd))    
# }

#' @title Cross validation for PLS regression
#' @description for internal use only!
#' @keywords internal
plsCv <- function(x, y, ncomp, 
                  method = c("pls", "wapls1"),
                  center = TRUE, scaled, 
                  minF = 1,
                  newX = matrix(0,1,1),
                  weights = NULL, 
                  p = 0.75, 
                  resampling = 10, 
                  group = NULL,
                  retrieve = TRUE,
                  .optimize = TRUE,
                  max.iter = 1, tol = 1e-6)
{
  
  if(min(ncol(x), nrow(x)) < floor(ncomp * (1 + p))) {ncomp <- (floor(min(ncol(x), nrow(x)) * p))-1} 
  
  #   smpl <- function(x){
  #     if(length(x) == 1)
  #       x <- c(x, x)
  #     return(sample(x, size = 1))
  #   }
  
  
  if(is.null(group)){
    nv <- floor((1-p)*nrow(x))
    y.str <- quantile(y, probs = seq(0, 1, length = (nv+1)))
    y.cut <- cut(y, unique(y.str), include.lowest = TRUE)
    levels(y.cut) <- 1:length(y)
    prob.g <- data.frame(orig.order = 1:length(y), prob.group = as.numeric (y.cut))
    
    rspi <- matrix(NA, nv, resampling)
    mrspi <- matrix(NA, nrow(x) - nv, resampling)
    colnames(rspi) <- colnames(mrspi) <- paste("Resample", seq(1:resampling), sep = "")
    rownames(rspi) <- paste("index", seq(1:nv), sep = "")
    rownames(mrspi) <- paste("index", seq(1:nrow(mrspi)), sep = "")
    
    idxs <- 1:nrow(x)
    for(jj in 1:resampling)
    {
      strs <- tapply(X = prob.g$orig.order, FUN = smpl, INDEX = prob.g$prob.group)
      if(length(strs) < nv)
      {
        adds <- sample((1:length(y))[-strs], nv - length(strs))
        strs <- c(strs, adds)
      }
      rspi[,jj]<- strs
      mrspi[,jj]<- idxs[!idxs %in% strs]
      
    }
  }else{
    ygr <- data.frame(y, grp = as.factor(as.character(group)), order.orig = 1:length(y))
    lvls <- levels(ygr$grp)
    nlv <- nlevels(ygr$grp)
    # Compute means per group and make the stratified sampling based on that 
    # vector of means
    yg <- tapply(X = ygr$y, INDEX = ygr$grp, FUN = mean)
    
    nv <- floor((1 - p) * nlv)
    y.str <- quantile(yg, probs = seq(0, 1, length = (nv+1)))
    y.cut <- cut(yg, unique(y.str), include.lowest = TRUE)
    levels(y.cut) <- paste("l_", 1:nlevels(y.cut), sep = "")
    prob.g <- data.frame(orig.order = 1:length(yg), prob.group = factor(y.cut))
    
    rspi <- matrix(0, nv, resampling)
    colnames(rspi) <- paste("Resample", seq(1:resampling), sep = "")
    rownames(rspi) <- paste("index", seq(1:nv), sep = "")
    
    listrspi <- mlistrspi <- NULL
    for(jj in 1:resampling)
    {
      strs <- tapply(X = prob.g$orig.order, FUN = smpl, INDEX = prob.g$prob.group)
      if(length(strs) < nv)
      {
        adds <- sample((1:nrow(prob.g))[-strs], nv - length(strs))
        strs <- c(strs, adds)
      }
      
      strs <- lvls[strs]  
      
      listrspi[[jj]] <- ygr$order.orig[ygr$grp %in% strs]
      mlistrspi[[jj]] <- ygr$order.orig[!ygr$grp %in% strs]
    }
    lensam <- lengths(listrspi)
    mlensam <- lengths(mlistrspi)
    rspi <- matrix(NA, max(lensam), resampling)
    mrspi <- matrix(NA, max(mlensam), resampling)
    colnames(rspi) <- colnames(mrspi) <- paste("Resample", seq(1:resampling), sep = "")
    rownames(rspi) <- paste("index", 1:nrow(rspi), sep = "")
    rownames(mrspi) <- paste("index", 1:nrow(mrspi), sep = "")
    idxs <- 1:nrow(x)
    for(jj in 1:resampling)
    {
      ## ramdomly select the replacements for the missing indexes
      ## replacements are necessary!
      rspi[,jj] <- c(sample(listrspi[[jj]], max(lensam) - lensam[jj]), listrspi[[jj]])
      mrspi[,jj] <- c(sample(mlistrspi[[jj]], max(mlensam) - mlensam[jj]), mlistrspi[[jj]])
    }
    
  }
  
  if(method == "wapls1" & retrieve & .optimize)
  {
    pls.c <- c(minF, ncomp)
    seq.pls <- min(pls.c):max(pls.c)
    sgr <- expand.grid(minF = seq.pls, maxF = seq.pls)
    sgr <- as.matrix(sgr[sgr$minF < sgr$maxF,])
    rownames(sgr) <- 1:nrow(sgr)
    fmethod <- "completewapls1"
  }else{
    sgr <- matrix(0,0,0)
    fmethod <- method
  }
  
  cvre <- pplscv_cpp(X = x, 
                     Y = as.matrix(y),
                     scale = scaled,
                     method = fmethod,
                     mindices = mrspi,
                     pindices = rspi,
                     minF = minF,
                     ncomp = ncomp,
                     newX = newX,
                     maxiter = max.iter, 
                     tol = tol,
                     waplsgrid = sgr)
  
  val <- NULL
  val$resamples <- rspi
  if(method == "pls"){
    val$cvResults <- data.frame(nPLS = 1:ncomp, rmse.cv = rowMeans(cvre$rmse.seg), st.rmse.cv = rowMeans(cvre$st.rmse.seg), rmse.sd.cv = cSds(t(cvre$rmse.seg)), rsq.cv = rowMeans(cvre$rsq.seg))
    bestpls.c <- val$cvResults$nPLS[which(val$cvResults$rmse.cv == min(val$cvResults$rmse.cv))]
    val$bestpls.c <- bestpls.c
    
    if(retrieve)
    {
      if(!is.null(weights)){
        x <- sweep(x, 1, weights, "*")  ### MODIFIED
        y <- y * weights
      }
      if(.optimize)
      {
        val$models <- fopls(X = x, 
                            Y = y, 
                            scale = scaled, 
                            ncomp = bestpls.c,
                            maxiter = max.iter, 
                            tol = tol)
      } else {
        val$models <- fopls(X = x, 
                            Y = as.matrix(y), 
                            scale = scaled, 
                            ncomp = ncomp,
                            maxiter = max.iter, 
                            tol = tol)
        
      }
    }
  }
  
  if(method == "wapls1" & retrieve & !.optimize){
    val$compweights <- cvre$compweights
    val$cvResults <- data.frame(minF = minF, maxF = ncomp, rmse.cv = mean(cvre$rmse.seg), st.rmse.cv = mean(cvre$st.rmse.seg), rmse.sd.cv = sd(cvre$rmse.seg), rsq.cv = mean(cvre$rsq.seg))
    #bestpls.c <- val$cvResults$nPLS[which(val$cvResults$rmse.cv == min(val$cvResults$rmse.cv))]
    #val$bestpls.c <- bestpls.c
    
    if(!is.null(weights)){
      x <- sweep(x, 1, weights, "*")  ### MODIFIED
      y <- y * weights
    } 
    val$models <- fopls(X = x, 
                        Y = as.matrix(y), 
                        scale = scaled, 
                        ncomp = ncomp,
                        maxiter = max.iter, 
                        tol = tol)
  }
  
  if(method == "wapls1" & retrieve & .optimize){
    
    val$cvResults <- data.frame(sgr, rmse.cv = rowMeans(cvre$rmse.seg), st.rmse.cv = rowMeans(cvre$st.rmse.seg), rmse.sd.cv = cSds(t(cvre$rmse.seg)), rsq.cv = rowMeans(cvre$rsq.seg))
    opmls <- which.min(val$cvResults$rmse.cv)
    val$bestpls.c$minF <- val$cvResults$minF[opmls]
    val$bestpls.c$maxF <- val$cvResults$maxF[opmls]
    val$compweights <- cvre$compweights
    
    
    if(!is.null(weights)){
      x <- sweep(x, 1, weights, "*")  ### MODIFIED
      y <- y * weights
    }
    
    tmp.compweights <- val$compweights[val$bestpls.c$minF:val$bestpls.c$maxF]
    val$compweights[] <- 0
    tmp.compweights <- tmp.compweights/sum(tmp.compweights)
    val$compweights[val$bestpls.c$minF:val$bestpls.c$maxF] <- tmp.compweights
    
    val$models <- fopls(X = x, 
                        Y = as.matrix(y), 
                        scale = scaled, 
                        ncomp = ncomp,
                        maxiter = max.iter, 
                        tol = tol)
  }
  return(val)
}


#' @title Internal function for computing the weights of the PLS components necessary for weighted average PLS 
#' @description internal
#' @keywords internal
#' @param plsO either an object returned by the \code{plsCv} function or an object as returned by the \code{opls} function which contains a pls model.
#' @param orgX the original spectral \code{matrix} which was used for calibrating the pls model.
#' @param type type of weight to be computed. The only available option (for the moment) is \code{"w1"}. See details on the \code{mbl} function where it is explained how \code{"w1"} is computed whitin the \code{"wapls1"} regression.
#' @param newX a \code{vector} of a new spectral sample. When "w1" is selected, newX must be specified.
#' @param pls.c a \code{vector} of length 2 which contains both the minimum and maximum number of PLS components for which the weights must be computed.
#' @return \code{wapls.weights} returns a \code{vector} of weights for each PLS component specified
#' @author Leonardo Ramirez-Lopez and Antoine Stevens

wapls.weights <- function(plsO, orgX, type = "w1", newX = NULL, pls.c){
  
  if(!is.null(newX)){newX <- t(as.vector(newX))}
  
  if(!is.null(plsO$models))
  {
    rmse.cv <- plsO$cvResults$rmse.cv
    plsO <- plsO$models
    is_plsCv <- TRUE
  } else is_plsCv <-FALSE
  if(max(pls.c) > max(plsO$ncomp))
    stop("the maximum number of PLS components specified in pls.c exceeds the number of PLS components contained in plsO!")
  if(type == "w1" & is.null(newX))
    stop("newX is missng! When type = 'w2', a vector of a new spectral sample must be specified")
  if(type == "w1" & (length(newX) != ncol(orgX)))
    stop("the length of the vector newX does not match with the number of variables of orgX")
  if(length(pls.c) != 2)
    stop("'pls.c' must be a vector of length 2 which specifies both the minimum and the maximum number of PLS components")
  if(pls.c[[1]] > pls.c[[2]])
    stop("for the vector of 'pls.c' the first value must be the minimum number of PLS components")
  
  minF <- pls.c[[1]]
  maxF <- pls.c[[2]]
  
  
  whgt <- waplswCpp(projectionm = plsO$projectionM, 
                    xloadings = plsO$X.loadings,
                    coefficients = plsO$coefficients,
                    newX = newX,
                    minF = minF, 
                    maxF = maxF, 
                    scale = ifelse(nrow(plsO$transf$Xscale) == 1, TRUE, FALSE),
                    Xcenter = plsO$transf$Xcenter,
                    Xscale = plsO$transf$Xscale)[1,]
  
  return(whgt)
}


#' @title iterator that re-organize and subset a matrix based on dissimilarity matrix (used in mbl)
#' @description internal
#' @param x an input \code{matrix}
#' @param Dx a (symmetric) cross-dissimilarity \code{matrix} of x
#' @param D a distance \code{matrix} between x and another matrix
#' @param k \code{matrix} number of neighbours (rows)  to keep
#' @param by \code{'column'} or \code{'row'}
#' @return an object of \code{class} iterator giving a \code{list} with (a) a subset of the input \code{matrix} with rows re-arranged in the order 
#' of the distance \code{vector}, (b) the distance \code{vector}, (c) the order of the distance \code{vector}
#' @details isubset will look at the order of knn in each col of D and re-organize the rows of x accordingly
#' @author Antoine Stevens
#' @keywords internal
isubset <- function(x, Dx = NULL, D, k, by = "column") {
  
  it <- iter(D, by= by)
  k <- iter(k, by = "row")
  nextEl <- function() {
    d <- nextElem(it)       
    kk <- nextElem(k)       
    n <- order(d)
    n1 <- n[1:max(kk)]
    if(is.null(Dx))
      list(x=x[n1,],d=d,index=n,k=c(kk)) # resulting iterator
    else
      list(x=x[n1,],dx=Dx[n1,n1],d=d,index=n,k=c(kk))
  }
  
  obj <- list(nextElem=nextEl)
  class(obj) <- c("isubset", "abstractiter", "iter")
  obj
}

#' @title Standard deviation of columns
#' @description For internal use only!
#' @keywords internal
cSds <-  function(x){
  return(as.vector(colSds(x)))
}