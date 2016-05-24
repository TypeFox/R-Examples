#' @include GLCM.R GLRLM.R GLSZM.R MGLSZM.R
NULL

#' Calculate texture and first order statistics.
#'
#' \code{calc_features} Calculates features of given texture matrix. If a simple
#' matrix is given, will calculate first order features. If desired, user may input the
#' features they wish to calculate for a given matrix type by passing them as a vector
#' of strings to the \code{features} argument.
#' 
#' Lists of features available for each matrix type can be accessed through
#' \code{?first_order_features} \code{?glcm_features}, \code{?glrlm_features}, \code{?glszm_features}.
#' 
#' Matrices of class \code{mglszm} accept features belonging to the glszm.
#'
#' @param object An object of class "matrix", "glcm", "glrlm", "glszm", or "mglszm"
#' @param features A vector containing the features the user wishes to calculate for a 
#' given matrix type. 
#' @return A data frame with a single observation. The columns of the dataframe 
#'   correspond to the calculated features.
#'
#' @examples
#' \dontrun{
#' calc_features(glcm(hallbey))
#' calc_features(glrlm(psf, n_grey=10))
#' calc_features(glcm(hallbey), features=c("glcm_mean", "glcm_variance", "pickles"))
#' }
#' 
#' @references \url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0102107}
#' 
#' @seealso \code{\link{glcm}}
#'   \code{\link{glrlm}}
#'   \code{\link{glszm}}
#'   \code{\link{mglszm}}
#' @export
#' @importFrom grDevices heat.colors
#' @importFrom graphics axis
#' @importFrom stats median
#' @importFrom stats sd
#' @importFrom stats var
setGeneric("calc_features", function(object, features = c()) standardGeneric("calc_features"))

#' @describeIn calc_features Calculate first order features of a numeric matrix 
#' @export
setMethod("calc_features", 
          signature = "matrix",
          definition = function(object, features = c()){
            
            #Set up allowed features
            feature_list<- list("calc_energy",
                                "calc_entropy",
                                "calc_kurtosis",
                                "calc_meanDeviation",
                                "calc_skewness",
                                "calc_uniformity",
                                "calc_mean",
                                "calc_median",
                                "calc_max",
                                "calc_min",
                                "calc_variance",
                                "calc_RMS",
                                "calc_sd"
            )
            
            #If user entered features, choose those in list that apply,
            #give message for those that dont
            
            if(length(features != 0)){
              if(any(!features %in% feature_list)){
                message("Some requested features not available for glcm:")
                message(paste0("    ", features[!features %in% feature_list], "\n"))
                message("For list of available features enter ?glcm_features")
              } 
              feature_list <- feature_list[feature_list %in% features]
            }
            
            feature_df <- data.frame(lapply(feature_list, function(f) tryCatch(get(f)(object),
                                                                               error=function(cond) return(NA),
                                                                               warning=function(cond) return(NA))))
            colnames(feature_df) <- feature_list
            return(feature_df)
          }
          
)     




#' @describeIn calc_features Calculate texture features of a glcm matrix 
#' @export
setMethod("calc_features", 
          signature = "glcm",
          definition = function(object, features = c()){
            
            #Set up allowed features
            feature_list<- list(
              "glcm_mean", "glcm_variance", "glcm_autoCorrelation",
              "glcm_cProminence", "glcm_cShade", "glcm_cTendency",
              "glcm_contrast", "glcm_correlation", "glcm_differenceEntropy",
              "glcm_dissimilarity", "glcm_energy", "glcm_entropy", 
              "glcm_homogeneity1", "glcm_homogeneity2", "glcm_IDMN",
              "glcm_IDN", "glcm_inverseVariance", "glcm_maxProb", 
              "glcm_sumAverage", "glcm_sumEntropy", "glcm_sumVariance"
            )
            
            #If user entered features, choose those in list that apply,
            #give message for those that dont
            
            if(length(features != 0)){
              if(any(!features %in% feature_list)){
                message("Some requested features not available for glcm:")
                message(paste0("    ", features[!features %in% feature_list], "\n"))
                message("For list of available features enter ?glcm_features")
              } 
              feature_list <- feature_list[feature_list %in% features]
            }
            feature_df <- data.frame(lapply(feature_list, function(f) tryCatch(get(f)(object),
                                            error=function(cond) return(NA),
                                            warning=function(cond) return(NA))))
            colnames(feature_df) <- feature_list
            return(feature_df)
          }
          
)     

#' @describeIn calc_features Calculate texture features of a glrlm matrix 
#' @export 
setMethod("calc_features", 
          signature = "glrlm",
          definition = function(object, features = c()){
            
            #Set up allowed features
            feature_list<- list("glrlm_GLN", "glrlm_HGLRE", "glrlm_LRE", 
                                "glrlm_LRHGLE", "glrlm_LRLGLE",
                                "glrlm_LGLRE", "glrlm_RLN", "glrlm_RP",
                                "glrlm_SRE", "glrlm_SRHGLE", "glrlm_SRLGLE"
            )
            
            #If user entered features, choose those in list that apply,
            #give message for those that dont
            if(length(features != 0)){
              if(any(!features %in% feature_list)){
                message("Some requested features not available for glrlm:")
                message(paste0("    ", features[!features %in% feature_list], "\n"))
                message("For list of available features enter ?glrlm_features")
              } 
              feature_list <- feature_list[feature_list %in% features]
            }
            feature_df <- data.frame(lapply(feature_list, function(f) tryCatch(get(f)(object),
                                            error=function(cond) return(NA),
                                            warning=function(cond) return(NA))))
            colnames(feature_df) <- feature_list
            return(feature_df)
          }
          
)   


#' @describeIn calc_features Calculate texture features of a glszm matrix 
#' @export
setMethod("calc_features", 
          signature = "glszm",
          definition = function(object, features = c()){
            
            #Set up allowed features
            feature_list<- list("glszm_SAE", "glszm_LAE", "glszm_IV", 
                                "glszm_SZV", "glszm_ZP", "glszm_LIE",
                                "glszm_HIE", "glszm_LISAE", "glszm_HISAE", 
                                "glszm_LILAE", "glszm_HILAE"
            )
            
            #If user entered features, choose those in list that apply,
            #give message for those that dont
            if(length(features != 0)){
              if(any(!features %in% feature_list)){
                message("Some requested features not available for glszm:")
                message(paste0("    ", features[!features %in% feature_list], "\n"))
                message("For list of available features enter ?glszm_features")
              } 
              feature_list <- feature_list[feature_list %in% features]
            }
            feature_df <- data.frame(lapply(feature_list, function(f) tryCatch(get(f)(object),
                                                                               error=function(cond) return(NA),
                                                                               warning=function(cond) return(NA))))
            
            colnames(feature_df) <- feature_list
            return(feature_df)
          }
          
)  


#' @describeIn calc_features Calculate texture features of an mglszm matrix 
#' @export
setMethod("calc_features", 
          signature = "mglszm",
          definition = function(object, features = c()){
            
            #Set up allowed features
            feature_list<- list("glszm_SAE", "glszm_LAE", "glszm_IV", 
                                "glszm_SZV", "glszm_ZP", "glszm_LIE",
                                "glszm_HIE", "glszm_LISAE", "glszm_HISAE", 
                                "glszm_LILAE", "glszm_HILAE"
            )
            
            #If user entered features, choose those in list that apply,
            #give message for those that dont
            if(length(features != 0)){
              if(any(!features %in% feature_list)){
                message("Some requested features not available for glszm:")
                message(paste0("    ", features[!features %in% feature_list], "\n"))
                message("For list of available features enter ?glszm_features")
              } 
              feature_list <- feature_list[feature_list %in% features]
            }
            feature_df <- data.frame(lapply(feature_list, function(f) tryCatch(get(f)(object),
                                                                               error=function(cond) return(NA),
                                                                               warning=function(cond) return(NA))))
            
            
            colnames(feature_df) <- paste0("m",feature_list)
            return(feature_df)
          }
          
)  










