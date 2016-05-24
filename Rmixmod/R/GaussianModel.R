###################################################################################
##                          GaussianModel.R                                      ##
###################################################################################

###################################################################################
##' @include global.R
##' @include Model.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{GaussianModel}}] class
##'
##' This class defines a gaussian Model. Inherits the [\code{\linkS4class{Model}}] class.
##' 
##' \describe{
##'   \item{family}{character defining a family of models.}
##' }
##'
##' @examples
##'   new("GaussianModel")
##'   new("GaussianModel",family="general")
##'
##'   getSlots("GaussianModel")
##' 
##' @name GaussianModel-class
##' @rdname GaussianModel-class
##' @exportClass GaussianModel
##'
setClass(
    Class="GaussianModel",
    representation=representation(
        family = "character"
    ),
    contains=c("Model"),
    prototype=prototype(
        family = character(0)
    ),
    validity=function(object){
      # spherical models
      spherical.free=c("Gaussian_pk_L_I", "Gaussian_pk_Lk_I")
      spherical.equal=c("Gaussian_p_L_I", "Gaussian_p_Lk_I")
      # all spherical
      spherical=c(spherical.free,spherical.equal)
      # diagonal models
      diagonal.free=c("Gaussian_pk_L_B", "Gaussian_pk_Lk_B", "Gaussian_pk_L_Bk", "Gaussian_pk_Lk_Bk")
      diagonal.equal=c("Gaussian_p_L_B", "Gaussian_p_Lk_B", "Gaussian_p_L_Bk", "Gaussian_p_Lk_Bk")
      # all diagonal
      diagonal=c(diagonal.free,diagonal.equal)
      #  general models
      general.free=c("Gaussian_pk_L_C", "Gaussian_pk_Lk_C", "Gaussian_pk_L_D_Ak_D", "Gaussian_pk_Lk_D_Ak_D", "Gaussian_pk_L_Dk_A_Dk", "Gaussian_pk_Lk_Dk_A_Dk", "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck")
      general.equal=c( "Gaussian_p_L_C", "Gaussian_p_Lk_C", "Gaussian_p_L_D_Ak_D", "Gaussian_p_Lk_D_Ak_D", "Gaussian_p_L_Dk_A_Dk", "Gaussian_p_Lk_Dk_A_Dk", "Gaussian_p_L_Ck", "Gaussian_p_Lk_Ck")
      # all general
      general=c(general.free,general.equal)
      # all models
      all.free=c(spherical.free,diagonal.free,general.free)
      all.equal=c(spherical.equal,diagonal.equal,general.equal)
      all=c(spherical,diagonal,general)
          
      # check listModels validity
      if ( sum(object@listModels %in% all) != length(object@listModels) )
          stop("At least one model is not a valid model. See ?mixmodGaussianModel for the list of all gaussian models.")
          
      # check proportions parameters validity
      if ( !object@equal.proportions & !object@free.proportions )
        stop("equal.proportions and free.porportions cannot be both as FALSE !")
      
      # check whether some family names are not valid
      if ( sum( !(object@family %in% c("all","general","diagonal","spherical")) ) ){
        warning( object@family[which(!(object@family %in% c("all","general","diagonal","spherical")))], ": unknown family name !")
      }
      
      # check proportions
      if ( !object@free.proportions & (sum(object@listModels %in% all.free)>0) )
        stop("At least one model has a free proportions but free.proportions is set as FALSE. See ?mixmodGaussianModel for the list of models with equal proportions.")
      if ( !object@equal.proportions & (sum(object@listModels %in% all.equal)>0) )
        stop("At least one model has an equal proportions but equal.proportions is set as FALSE. See ?mixmodGaussianModel for the list of models with free proportions.")
      
      return(TRUE)
    }
)
###################################################################################

###################################################################################
##' Create an instance of the [\code{\linkS4class{GaussianModel}}] class using new/initialize.
##' 
##' Initialization method. Used internally in the `Rmixmod' package.
##' 
##' @seealso \code{\link{initialize}}
##'
##' @keywords internal
##'
##' @rdname initialize-methods
##'
setMethod(
  f="initialize",
  signature=c("GaussianModel"),
  definition=function(.Object, listModels, family, free.proportions, equal.proportions){
          
    # spherical models
    spherical.free=c("Gaussian_pk_L_I", "Gaussian_pk_Lk_I")
    spherical.equal=c("Gaussian_p_L_I", "Gaussian_p_Lk_I")
    # all spherical
    spherical=c(spherical.free,spherical.equal)
    
    # diagonal models
    diagonal.free=c("Gaussian_pk_L_B", "Gaussian_pk_Lk_B", "Gaussian_pk_L_Bk", "Gaussian_pk_Lk_Bk")
    diagonal.equal=c("Gaussian_p_L_B", "Gaussian_p_Lk_B", "Gaussian_p_L_Bk", "Gaussian_p_Lk_Bk")
    # all diagonal
    diagonal=c(diagonal.free,diagonal.equal)
    
    #  general models
    general.free=c("Gaussian_pk_L_C", "Gaussian_pk_Lk_C", "Gaussian_pk_L_D_Ak_D", "Gaussian_pk_Lk_D_Ak_D", "Gaussian_pk_L_Dk_A_Dk", "Gaussian_pk_Lk_Dk_A_Dk", "Gaussian_pk_L_Ck", "Gaussian_pk_Lk_Ck")
    general.equal=c( "Gaussian_p_L_C", "Gaussian_p_Lk_C", "Gaussian_p_L_D_Ak_D", "Gaussian_p_Lk_D_Ak_D", "Gaussian_p_L_Dk_A_Dk", "Gaussian_p_Lk_Dk_A_Dk", "Gaussian_p_L_Ck", "Gaussian_p_Lk_Ck")
    # all general
    general=c(general.free,general.equal)
    
    # all models
    all.free=c(spherical.free,diagonal.free,general.free)
    all.equal=c(spherical.equal,diagonal.equal,general.equal)
    all=c(spherical,diagonal,general)
    
    if ( !missing(listModels) ){
      # save the list of models
      .Object@listModels <- listModels
      # set family 
      if ( missing(family) ){
        if ( sum(listModels %in% spherical) & sum(listModels %in% diagonal) & sum(listModels %in% general) ){
          .Object@family<-"all"
        }
        else{
          family<-character(0)
          if ( sum(listModels %in% spherical) ){ family<-c(family,"spherical") }
          if ( sum(listModels %in% diagonal) ){ family<-c(family,"diagonal") }
          if ( sum(listModels %in% general) ){ family<-c(family,"general") }
          .Object@family<-family
        }
      }
      else { .Object@family<-family }
      
      # set free.proportions
      if ( missing(free.proportions) ){
        if ( sum(listModels %in% all.free) ){ .Object@free.proportions<-TRUE }
        else{ .Object@free.proportions<-FALSE }
      }
      else{ .Object@free.proportions<-free.proportions }
      # set equal.proportions
      if ( missing(equal.proportions) ){
        if ( sum(listModels %in% all.equal) ){ .Object@equal.proportions<-TRUE }
        else{ .Object@equal.proportions<-FALSE }
      }
      else{ .Object@equal.proportions<-equal.proportions }
    }
    else{
      # check free.proportions option
      if ( missing(free.proportions) ){ .Object@free.proportions<-TRUE }
      else{ .Object@free.proportions<-free.proportions }
      
      # check equal.proportions option
      if ( missing(equal.proportions) ){ .Object@equal.proportions<-TRUE }
      else{ .Object@equal.proportions<-equal.proportions }
      
      # define an empty list of models
      list<-character(0)
      
      # set family as "all" if missing
      if ( missing(family) ){
        .Object@family<-"all"
        if ( .Object@free.proportions ){ list<-c(list,all.free) }
        if ( .Object@equal.proportions ){ list<-c(list,all.equal) }
      }
      else{
        # all gaussian models
        if (sum(family=="all")){
          # set family label in case of multiple entries
          .Object@family = "all"
          if ( .Object@free.proportions ){ list<-c(list,all.free) }
          if ( .Object@equal.proportions ){ list<-c(list,all.equal) }
        }
        else{
          # all spherical models
          if ( sum(family=="spherical") ){
            if ( .Object@free.proportions ){ list<-c(list,spherical.free) }
            if ( .Object@equal.proportions ){ list<-c(list,spherical.equal) }
          }
          # all diagonal models
          if ( sum(family=="diagonal") ){
            .Object@family = "diagonal"
            if ( .Object@free.proportions ){ list<-c(list,diagonal.free) }
            if ( .Object@equal.proportions ){ list<-c(list,diagonal.equal) }
          }
          # all general models
          if ( sum(family=="general") ){
            .Object@family = "general"
            if ( .Object@free.proportions ){ list<-c(list,general.free) }
            if ( .Object@equal.proportions ){ list<-c(list,general.equal) }
          }
          # set family label in case of multiple entries
          .Object@family = family
        }
      }
      # create the list of models depending on the proportions option
      .Object@listModels<-list
    }
    
    validObject(.Object)        
    return(.Object)
  }
)
###################################################################################


###################################################################################
##' Create an instance of the [\code{\linkS4class{GaussianModel}}] class
##'
##' Define a list of Gaussian model to test in MIXMOD.
##' 
##' In the Gaussian mixture model, following Banfield and Raftery (1993) and Celeux and Govaert (1995), we consider a parameterization of the variance matrices of the mixture components consisting of expressing the variance matrix \eqn{\Sigma_{k}} in terms of its eigenvalue decomposition \deqn{ \Sigma_{k}= \lambda_{k} D_{k} A_{k}D'_{k}} where \eqn{\lambda_{k}=|\Sigma_{k}|^{1/d}, D_{k}} is the matrix of eigenvectors of \eqn{\Sigma_{k}} and \eqn{A_{k}} is a diagonal matrix, such that \eqn{| A_{k} |=1}, with the normalized eigenvalues of \eqn{\Sigma_{k}} on the diagonal in a decreasing order. The parameter \eqn{\lambda_{k}} determines the \emph{volume} of the \eqn{k}th cluster, \eqn{D_{k}} its \emph{orientation} and \eqn{A_{k}} its \emph{shape}.  By allowing some but not all of these quantities to vary between clusters, we obtain parsimonious and easily interpreted models which are appropriate to describe various clustering situations.
##' 
##' In general family, we can allow the volumes, the shapes and the orientations of clusters to vary or to be equal between clusters. Variations on assumptions on the parameters \eqn{\lambda_{k}, D_{k}} and \eqn{A_{k}} \eqn{(1 \leq k \leq K)} lead to 8 general models of interest. For instance, we can assume different volumes and keep the shapes and orientations equal by requiring that \eqn{A_{k}=A} (\eqn{A} unknown) and \eqn{D_{k}=D} (\eqn{D} unknown) for \eqn{k=1,\ldots,K}. We denote this model \eqn{[\lambda_{k}DAD']}. With this convention, writing \eqn{[\lambda D_{k}AD'_{k}]} means that we consider the mixture model with equal volumes, equal shapes and different orientations.
##' In diagonal family, we assume that the variance matrices \eqn{\Sigma_{k}} are diagonal. In the parameterization, it means that the orientation matrices \eqn{D_{k}} are permutation matrices. We write \eqn{\Sigma_{k}=\lambda_{k}B_{k}} where \eqn{B_{k}} is a diagonal matrix with \eqn{| B_{k}|=1}.  This particular parameterization gives rise to 4 models: \eqn{[\lambda B]}, \eqn{[\lambda_{k}B]}, \eqn{[\lambda B_{k}]} and \eqn{[\lambda_{k}B_{k}]}.
##' 
##' In spherical family, we assume spherical shapes, namely \eqn{A_{k}=I}, \eqn{I} denoting the identity matrix. In such a case, two parsimonious models are in competition: \eqn{[\lambda I]} and \eqn{[\lambda_{k}I]}.
##'
##' @param family character defining a family of models. "general" for the general family, "diagonal" for the diagonal family, "spherical" for the spherical family and "all" for all families. Default is "general".
##' @param listModels a list of characters containing a list of models. It is optional.
##' @param free.proportions logical to include models with free proportions. Default is TRUE.
##' @param equal.proportions logical to include models with equal proportions. Default is TRUE.
##'
##' @return an object of [\code{\linkS4class{GaussianModel}}] which contains some of the 28 Gaussian Models:
##' \tabular{rlllll}{
##'     Model  \tab Family \tab Prop. \tab Volume \tab Shape \tab Orient. \cr
##'     Gaussian_p_L_C         \tab General \tab Equal \tab Equal \tab Equal  \tab Equal \cr
##'     Gaussian_p_Lk_C        \tab \tab \tab Free \tab Equal \tab Equal \cr
##'     Gaussian_p_L_D_Ak_D    \tab  \tab \tab Equal \tab Free \tab Equal \cr
##'     Gaussian_p_Lk_D_Ak_D   \tab  \tab \tab Free \tab Free \tab Equal \cr
##'     Gaussian_p_L_Dk_A_Dk   \tab  \tab \tab Equal \tab Equal \tab Free \cr
##'     Gaussian_p_Lk_Dk_A_Dk  \tab  \tab \tab Free \tab Equal \tab Free \cr
##'     Gaussian_p_L_Ck        \tab  \tab \tab Equal \tab Free \tab Free \cr
##'     Gaussian_p_Lk_Ck       \tab  \tab \tab Free \tab Free \tab Free \cr
##'     Gaussian_p_L_B         \tab Diagonal  \tab Equal \tab Equal \tab Equal \tab Axes \cr
##'     Gaussian_p_Lk_B        \tab  \tab \tab Free \tab Equal \tab Axes \cr
##'     Gaussian_p_L_Bk        \tab  \tab \tab Equal \tab Free \tab Axes \cr
##'     Gaussian_p_Lk_Bk       \tab  \tab \tab Free \tab Free \tab Axes \cr
##'     Gaussian_p_L_I         \tab Spherical \tab Equal \tab Equal \tab Equal \tab NA \cr
##'     Gaussian_p_Lk_I        \tab \tab \tab Free \tab Equal \tab NA \cr
##'     Gaussian_pk_L_C        \tab General  \tab Free \tab Equal \tab Equal \tab Equal \cr
##'     Gaussian_pk_Lk_C       \tab \tab \tab Free \tab Equal \tab Equal \cr
##'     Gaussian_pk_L_D_Ak_D   \tab \tab \tab Equal \tab Free \tab Equal \cr
##'     Gaussian_pk_Lk_D_Ak_D  \tab \tab \tab Free \tab Free \tab Equal \cr
##'     Gaussian_pk_L_Dk_A_Dk  \tab \tab \tab Equal \tab Equal \tab Free \cr
##'     Gaussian_pk_Lk_Dk_A_Dk \tab \tab \tab Free \tab Equal \tab Free \cr
##'     Gaussian_pk_L_Ck       \tab \tab \tab Equal \tab Free \tab Free \cr
##'     Gaussian_pk_Lk_Ck      \tab \tab \tab Free \tab Free \tab Free \cr
##'     Gaussian_pk_L_B        \tab  Diagonal  \tab Free \tab Equal \tab Equal \tab Axes \cr
##'     Gaussian_pk_Lk_B       \tab \tab \tab Free \tab Equal \tab Axes \cr
##'     Gaussian_pk_L_Bk       \tab \tab \tab Equal \tab Free \tab Axes \cr
##'     Gaussian_pk_Lk_Bk      \tab \tab \tab Free \tab Free \tab Axes \cr
##'     Gaussian_pk_L_I        \tab Spherical  \tab Free \tab Equal \tab Equal \tab NA \cr
##'     Gaussian_pk_Lk_I       \tab \tab \tab Free \tab Equal \tab NA \cr 
##' }
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##' @examples
##'   mixmodGaussianModel()
##'   # all Gaussian models with equal proportions
##'   mixmodGaussianModel(family="all",free.proportions=FALSE)
##'   # Diagonal and Spherical Gaussian models
##'   mixmodGaussianModel(family=c("diagonal","spherical"))
##'   # Gaussian models with a pre-defined list
##'   mixmodGaussianModel(listModels=c("Gaussian_p_L_C","Gaussian_p_L_Ck","Gaussian_pk_L_I"))
##'
##' @author Remi Lebret and Serge Iovleff and Florent Langrognet, with contributions from C. Biernacki and G. Celeux and G. Govaert \email{contact@@mixmod.org}
##' @export
##'
mixmodGaussianModel<- function( family="all", listModels=NULL, free.proportions=TRUE, equal.proportions=TRUE ){
    if ( is.null(listModels) ){
      new("GaussianModel", family=family, free.proportions=free.proportions, equal.proportions=equal.proportions)
    }else{
      new("GaussianModel", listModels=listModels)
    }
}
###################################################################################



###################################################################################
##' @rdname extract-methods
##' @aliases [,GaussianModel-method
##'
setMethod(
  f="[", 
  signature(x = "GaussianModel"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "listModels"={return(x@listModels)},
        "free.proportions"={return(x@free.proportions)},
        "equal.proportions"={return(x@equal.proportions)},
        "family"={return(x@family)},
        stop("This attribute doesn't exist !")
      )
    }else{
      switch(EXPR=i,
        "listModels"={return(x@listModels[j])},
        stop("This attribute doesn't exist !")
      )
    }
  }
)
###################################################################################



###################################################################################
##' @name [
##' @rdname extract-methods
##' @aliases [<-,GaussianModel-method
##'
setReplaceMethod(
  f="[", 
  signature(x = "GaussianModel"), 
  definition=function(x,i,j,value){
    if ( missing(j) ){
      switch(EXPR=i,
        "listModels"={x@listModels<-value},
        "free.proportions"={x@free.proportions<-value},
        "equal.proportions"={x@equal.proportions<-value},
        "family"={return(x@family)},
        stop("This attribute doesn't exist !")
      )
    }else{
      switch(EXPR=i,
        "listModels"={x@listModels[j]<-value},
        stop("This attribute doesn't exist !")
      )
    }
    validObject(x)
    return(x)
  }
)
###################################################################################
