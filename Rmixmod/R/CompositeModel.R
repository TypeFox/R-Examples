###################################################################################
##                          CompositeModel.R                                      ##
###################################################################################

###################################################################################
##' @include global.R
##' @include Model.R
##' @include GaussianModel.R
##' @include MultinomialModel.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{CompositeModel}}] class
##'
##' This class defines a Composite Model. Inherits the [\code{\linkS4class{Model}}] class.
##' 
##' \describe{
##'   \item{variable.independency}{logical}
##'   \item{component.independency}{logical}
##' }
##'
##' @examples
##'   new("CompositeModel")
##'   new("CompositeModel", listModels=c("Heterogeneous_pk_E_L_B","Heterogeneous_pk_Ekj_L_B") )
##'   new("CompositeModel", free.proportions=FALSE, variable.independency=TRUE )
##'
##'   getSlots("CompositeModel")
##' 
##' @name CompositeModel-class
##' @rdname CompositeModel-class
##' @exportClass CompositeModel
##'
setClass(
  Class="CompositeModel",
  representation=representation(
      variable.independency = "logical",
      component.independency = "logical"
  ),
  contains=c("Model"),
  prototype=prototype(
      variable.independency = logical(0),
      component.independency = logical(0)
  ),    
  validity=function(object){
        
    # define list of models
    vcf<-c("Heterogeneous_pk_E_L_B",
           "Heterogeneous_pk_E_Lk_B",
           "Heterogeneous_pk_E_L_Bk",
           "Heterogeneous_pk_E_Lk_Bk")
    vce<-c("Heterogeneous_p_E_L_B",
           "Heterogeneous_p_E_Lk_B",
           "Heterogeneous_p_E_L_Bk",
           "Heterogeneous_p_E_Lk_Bk")
    f<-c("Heterogeneous_pk_Ekj_L_B",
         "Heterogeneous_pk_Ekj_Lk_B",
         "Heterogeneous_pk_Ekj_L_Bk",
         "Heterogeneous_pk_Ekj_Lk_Bk",
         "Heterogeneous_pk_Ekjh_L_B",
         "Heterogeneous_pk_Ekjh_Lk_B",
         "Heterogeneous_pk_Ekjh_L_Bk",
         "Heterogeneous_pk_Ekjh_Lk_Bk")
    e<-c("Heterogeneous_p_Ekj_L_B",
         "Heterogeneous_p_Ekj_Lk_B",
         "Heterogeneous_p_Ekj_L_Bk",
         "Heterogeneous_p_Ekj_Lk_Bk",
         "Heterogeneous_p_Ekjh_L_B",
         "Heterogeneous_p_Ekjh_Lk_B",
         "Heterogeneous_p_Ekjh_L_Bk",
         "Heterogeneous_p_Ekjh_Lk_Bk")
    cf<-c("Heterogeneous_pk_Ej_L_B",
          "Heterogeneous_pk_Ej_Lk_B",
          "Heterogeneous_pk_Ej_L_Bk",
          "Heterogeneous_pk_Ej_Lk_Bk")
    ce<-c("Heterogeneous_p_Ej_L_B",
          "Heterogeneous_p_Ej_Lk_B",
          "Heterogeneous_p_Ej_L_Bk",
          "Heterogeneous_p_Ej_Lk_Bk")
    vf<-c("Heterogeneous_pk_Ek_L_B",
          "Heterogeneous_pk_Ek_Lk_B",
          "Heterogeneous_pk_Ek_L_Bk",
          "Heterogeneous_pk_Ek_Lk_Bk")
    ve<-c("Heterogeneous_p_Ek_L_B",
          "Heterogeneous_p_Ek_Lk_B",
          "Heterogeneous_p_Ek_L_Bk",
          "Heterogeneous_p_Ek_Lk_Bk")
    all.free<-c(vcf,f,cf,vf)
    all.equal<-c(vce,e,ce,ve)
    variable.free<-c(vcf,vf)
    variable.equal<-c(vce,ve)
    variable<-c(variable.free,variable.equal)
    component.free<-c(vcf,cf)
    component.equal<-c(vce,ce)
    component<-c(component.free,component.equal)
    
    # all models
    all=c(all.free,all.equal)
              
    # check listModels validity
    if ( sum(object@listModels %in% all) != length(object@listModels) )
      stop("At least one model is not a valid model. See ?mixmodCompositeModel for the list of all composite models.")

    # check proportions parameters validity
    if ( !object@equal.proportions & !object@free.proportions )
      stop("equal.proportions and free.porportions cannot be both as FALSE !")
    
    if ( !object@free.proportions & (sum(object@listModels %in% all.free)>0) )
      stop("At least one model has a free proportions but free.proportions is set as FALSE. See ?mixmodCompositeModel for the list of models with equal proportions.")
    
    if ( !object@equal.proportions & (sum(object@listModels %in% all.equal)>0) )
      stop("At least one model has an equal proportions but equal.proportions is set as FALSE. See ?mixmodCompositeModel for the list of models with free proportions.")

    # check independencies parameters
    # for variable
    if ( length(object@variable.independency) ){
      if ( object@variable.independency & sum(object@listModels %in% variable) != length(object@listModels) )
        stop("At least one model is not independent of the variable j. See ?mixmodCompositeModel for the list of all composite models.")
    }
    # for component
    if ( length(object@component.independency) ){
      if ( object@component.independency & sum(object@listModels %in% component) != length(object@listModels) )
        stop("At least one model is not independent of the variable j. See ?mixmodCompositeModel for the list of all composite models.")
    }
  }
)
###################################################################################

###################################################################################
##' Create an instance of the [\code{\linkS4class{CompositeModel}}] class using new/initialize.
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
  signature=c("CompositeModel"),
  definition=function(.Object, listModels, free.proportions, equal.proportions, variable.independency, component.independency){

      # define list of models
    vcf<-c("Heterogeneous_pk_E_L_B",
           "Heterogeneous_pk_E_Lk_B",
           "Heterogeneous_pk_E_L_Bk",
           "Heterogeneous_pk_E_Lk_Bk")
    vce<-c("Heterogeneous_p_E_L_B",
           "Heterogeneous_p_E_Lk_B",
           "Heterogeneous_p_E_L_Bk",
           "Heterogeneous_p_E_Lk_Bk")
    f<-c("Heterogeneous_pk_Ekj_L_B",
         "Heterogeneous_pk_Ekj_Lk_B",
         "Heterogeneous_pk_Ekj_L_Bk",
         "Heterogeneous_pk_Ekj_Lk_Bk",
         "Heterogeneous_pk_Ekjh_L_B",
         "Heterogeneous_pk_Ekjh_Lk_B",
         "Heterogeneous_pk_Ekjh_L_Bk",
         "Heterogeneous_pk_Ekjh_Lk_Bk")
    e<-c("Heterogeneous_p_Ekj_L_B",
         "Heterogeneous_p_Ekj_Lk_B",
         "Heterogeneous_p_Ekj_L_Bk",
         "Heterogeneous_p_Ekj_Lk_Bk",
         "Heterogeneous_p_Ekjh_L_B",
         "Heterogeneous_p_Ekjh_Lk_B",
         "Heterogeneous_p_Ekjh_L_Bk",
         "Heterogeneous_p_Ekjh_Lk_Bk")
    cf<-c("Heterogeneous_pk_Ej_L_B",
          "Heterogeneous_pk_Ej_Lk_B",
          "Heterogeneous_pk_Ej_L_Bk",
          "Heterogeneous_pk_Ej_Lk_Bk")
    ce<-c("Heterogeneous_p_Ej_L_B",
          "Heterogeneous_p_Ej_Lk_B",
          "Heterogeneous_p_Ej_L_Bk",
          "Heterogeneous_p_Ej_Lk_Bk")
    vf<-c("Heterogeneous_pk_Ek_L_B",
          "Heterogeneous_pk_Ek_Lk_B",
          "Heterogeneous_pk_Ek_L_Bk",
          "Heterogeneous_pk_Ek_Lk_Bk")
    ve<-c("Heterogeneous_p_Ek_L_B",
          "Heterogeneous_p_Ek_Lk_B",
          "Heterogeneous_p_Ek_L_Bk",
          "Heterogeneous_p_Ek_Lk_Bk")
    
    all.free<-c(vcf,f,cf,vf)
    all.equal<-c(vce,e,ce,ve)
    
    variable.free<-c(vcf,vf)
    variable.equal<-c(vce,ve)
    variable<-c(variable.free,variable.equal)
    
    component.free<-c(vcf,cf)
    component.equal<-c(vce,ce)
    component<-c(component.free,component.equal)
    
    if ( !missing(listModels) ){
      # save the list of models
      .Object@listModels <- listModels
      
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
      
      # set variable.independency
      if ( missing(variable.independency) ){
        if ( sum(listModels %in% variable) == length(listModels) ){ .Object@variable.independency<-TRUE }
      }
      else{ .Object@variable.independency<-variable.independency }
      # set component.independency
      if ( missing(component.independency) ){
        if ( sum(listModels %in% component) == length(listModels)  ){ .Object@component.independency<-TRUE }
      }
      else{ .Object@component.independency<-component.independency }
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
    
      if ( !missing(variable.independency) & !missing(component.independency)){
        if ( variable.independency & component.independency){
          if ( .Object@free.proportions ){ list<-c(list,vcf) }
          if ( .Object@equal.proportions ){ list<-c(list,vce) }
        }else if ( !variable.independency & !component.independency){
          if ( .Object@free.proportions ){ list<-c(list,f) }
          if ( .Object@equal.proportions ){ list<-c(list,e) }
        }else if ( !variable.independency & component.independency){
          if ( .Object@free.proportions ){ list<-c(list,cf) }
          if ( .Object@equal.proportions ){ list<-c(list,ce) }
        }else if ( variable.independency & !component.independency){
          if ( .Object@free.proportions ){ list<-c(list,vf) }
          if ( .Object@equal.proportions ){ list<-c(list,ve) }
        }
        .Object@component.independency <- component.independency
        .Object@variable.independency <- variable.independency
      }
      else if ( !missing(component.independency) ){
        if ( component.independency ){
          if ( .Object@free.proportions ){ list<-c(list,component.free) }
          if ( .Object@equal.proportions ){ list<-c(list,component.equal) }
        }else{
          if ( .Object@free.proportions ){ list<-c(list,f,vf) }
          if ( .Object@equal.proportions ){ list<-c(list,e,ve) }
        }
        .Object@component.independency<-component.independency
        .Object@variable.independency <- logical(0)
      }
      else if ( !missing(variable.independency) ){
        if ( variable.independency ){
          if ( .Object@free.proportions ){ list<-c(list,variable.free) }
          if ( .Object@equal.proportions ){ list<-c(list,variable.equal) }
        }else{
          if ( .Object@free.proportions ){ list<-c(list,f,cf) }
          if ( .Object@equal.proportions ){ list<-c(list,c,ce) }
        }
        .Object@component.independency <- logical(0)
        .Object@variable.independency <- variable.independency
      }
      else{
        # all multinomial models with free proportions
        if ( .Object@free.proportions ){ list<-c(list,all.free) }
        # all multinomial models with equal proportions
        if ( .Object@equal.proportions ){ list<-c(list,all.equal) }
        .Object@component.independency <- logical(0)
        .Object@variable.independency <- logical(0)
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
##' Create an instance of the [\code{\linkS4class{CompositeModel}}] class
##'
##' Define a list of heterogeneous model to test in MIXMOD.
##'
##' In heterogeneous case, Gaussian model can only belong to the diagonal family. We assume that the variance matrices \eqn{\Sigma_{k}} are diagonal. In the parameterization, it means that the orientation matrices \eqn{D_{k}} are permutation matrices. We write \eqn{\Sigma_{k}=\lambda_{k}B_{k}} where \eqn{B_{k}} is a diagonal matrix with \eqn{| B_{k}|=1}.  This particular parameterization gives rise to 4 models: \eqn{[\lambda B]}, \eqn{[\lambda_{k}B]}, \eqn{[\lambda B_{k}]} and \eqn{[\lambda_{k}B_{k}]}.
##' The multinomial distribution is associated to the \eqn{j}th variable of the \eqn{k}th component is reparameterized by a center \eqn{a_k^j} and the dispersion \eqn{\varepsilon_k^j} around this center. Thus, it allows us to give an interpretation similar to the center and the variance matrix used for continuous data in the Gaussian mixture context. In the following, this model will be denoted by \eqn{[\varepsilon_k^j]}. In this context, three other models can be easily deduced. We note \eqn{[\varepsilon_k]} the model where \eqn{\varepsilon_k^j} is independent of the variable \eqn{j}, \eqn{[\varepsilon^j]} the model where \eqn{\varepsilon_k^j} is independent of the component \eqn{k} and, finally, \eqn{[\varepsilon]} the model where \eqn{\varepsilon_k^j} is independent of both the variable $j$ and the component \eqn{k}.  In order to maintain some unity in the notation, we will denote also \eqn{[\varepsilon_k^{jh}]} the most general model introduced at the previous section.
##'
##' @param listModels a list of characters containing a list of models. It is optional.
##' @param free.proportions logical to include models with free proportions. Default is TRUE.
##' @param equal.proportions logical to include models with equal proportions. Default is TRUE.
##' @param variable.independency logical to include models where \eqn{[\varepsilon_k^j]} is independent of the variable \eqn{j}. Optionnal.
##' @param component.independency logical to include models where \eqn{[\varepsilon_k^j]} is independent of the component \eqn{k}. Optionnal.
##'
##' @return an object of [\code{\linkS4class{CompositeModel}}] which contains some of the 40 heterogeneous Models:
##' \tabular{rlllll}{
##'     Model  \tab Prop. \tab Var. \tab Comp. \tab Volume \tab Shape \cr
##'     Heterogeneous_p_E_L_B      \tab Equal  \tab TRUE \tab TRUE \tab Equal  \tab Equal \cr
##'     Heterogeneous_p_E_Lk_B     \tab \tab TRUE \tab TRUE \tab Free  \tab Equal \cr
##'     Heterogeneous_p_E_L_Bk     \tab \tab TRUE \tab TRUE \tab Equal  \tab Free \cr
##'     Heterogeneous_p_E_Lk_Bk    \tab  \tab TRUE \tab TRUE \tab Free  \tab Free \cr
##'     Heterogeneous_p_Ek_L_B     \tab  \tab TRUE \tab FALSE \tab Equal  \tab Equal \cr
##'     Heterogeneous_p_Ek_Lk_B    \tab  \tab TRUE \tab FALSE \tab Free  \tab Equal \cr
##'     Heterogeneous_p_Ek_L_Bk    \tab  \tab TRUE \tab FALSE \tab Equal  \tab Free \cr
##'     Heterogeneous_p_Ek_Lk_Bk   \tab   \tab TRUE \tab FALSE \tab Free  \tab Free \cr
##'     Heterogeneous_p_Ej_L_B     \tab  \tab FALSE \tab TRUE \tab Equal  \tab Equal \cr
##'     Heterogeneous_p_Ej_Lk_B    \tab  \tab FALSE \tab TRUE \tab Free  \tab Equal \cr
##'     Heterogeneous_p_Ej_L_Bk    \tab   \tab FALSE \tab TRUE \tab Equal  \tab Free \cr
##'     Heterogeneous_p_Ej_Lk_Bk   \tab   \tab FALSE \tab TRUE \tab Free  \tab Free \cr
##'     Heterogeneous_p_Ekj_L_B    \tab  \tab FALSE \tab FALSE \tab Equal  \tab Equal \cr
##'     Heterogeneous_p_Ekj_Lk_B   \tab   \tab FALSE \tab FALSE \tab Free  \tab Equal \cr
##'     Heterogeneous_p_Ekj_L_Bk   \tab   \tab FALSE \tab FALSE \tab Equal  \tab Free \cr
##'     Heterogeneous_p_Ekj_Lk_Bk  \tab   \tab FALSE \tab FALSE \tab Free  \tab Free \cr
##'     Heterogeneous_p_Ekjh_L_B   \tab  \tab FALSE \tab FALSE \tab Equal  \tab Equal \cr
##'     Heterogeneous_p_Ekjh_Lk_B  \tab   \tab FALSE \tab FALSE \tab Free  \tab Equal \cr
##'     Heterogeneous_p_Ekjh_L_Bk  \tab  \tab FALSE \tab FALSE \tab Equal  \tab Free \cr
##'     Heterogeneous_p_Ekjh_Lk_Bk \tab  \tab FALSE \tab FALSE \tab Free  \tab Free \cr
##'     Heterogeneous_pk_E_L_B      \tab Free  \tab TRUE \tab TRUE \tab Equal  \tab Equal \cr
##'     Heterogeneous_pk_E_Lk_B     \tab \tab TRUE \tab TRUE \tab Free  \tab Equal \cr
##'     Heterogeneous_pk_E_L_Bk     \tab \tab TRUE \tab TRUE \tab Equal  \tab Free \cr
##'     Heterogeneous_pk_E_Lk_Bk    \tab  \tab TRUE \tab TRUE \tab Free  \tab Free \cr
##'     Heterogeneous_pk_Ek_L_B     \tab  \tab TRUE \tab FALSE \tab Equal  \tab Equal \cr
##'     Heterogeneous_pk_Ek_Lk_B    \tab  \tab TRUE \tab FALSE \tab Free  \tab Equal \cr
##'     Heterogeneous_pk_Ek_L_Bk    \tab  \tab TRUE \tab FALSE \tab Equal  \tab Free \cr
##'     Heterogeneous_pk_Ek_Lk_Bk   \tab   \tab TRUE \tab FALSE \tab Free  \tab Free \cr
##'     Heterogeneous_pk_Ej_L_B     \tab  \tab FALSE \tab TRUE \tab Equal  \tab Equal \cr
##'     Heterogeneous_pk_Ej_Lk_B    \tab  \tab FALSE \tab TRUE \tab Free  \tab Equal \cr
##'     Heterogeneous_pk_Ej_L_Bk    \tab   \tab FALSE \tab TRUE \tab Equal  \tab Free \cr
##'     Heterogeneous_pk_Ej_Lk_Bk   \tab   \tab FALSE \tab TRUE \tab Free  \tab Free \cr
##'     Heterogeneous_pk_Ekj_L_B    \tab  \tab FALSE \tab FALSE \tab Equal  \tab Equal \cr
##'     Heterogeneous_pk_Ekj_Lk_B   \tab   \tab FALSE \tab FALSE \tab Free  \tab Equal \cr
##'     Heterogeneous_pk_Ekj_L_Bk   \tab   \tab FALSE \tab FALSE \tab Equal  \tab Free \cr
##'     Heterogeneous_pk_Ekj_Lk_Bk  \tab   \tab FALSE \tab FALSE \tab Free  \tab Free \cr
##'     Heterogeneous_pk_Ekjh_L_B   \tab  \tab FALSE \tab FALSE \tab Equal  \tab Equal \cr
##'     Heterogeneous_pk_Ekjh_Lk_B  \tab   \tab FALSE \tab FALSE \tab Free  \tab Equal \cr
##'     Heterogeneous_pk_Ekjh_L_Bk  \tab  \tab FALSE \tab FALSE \tab Equal  \tab Free \cr
##'     Heterogeneous_pk_Ekjh_Lk_Bk \tab  \tab FALSE \tab FALSE \tab Free  \tab Free \cr
##' }
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##' @examples
##'   mixmodCompositeModel()
##'   # composite models with equal proportions
##'   mixmodCompositeModel(free.proportions=FALSE)
##'   # composite models with equal proportions and independent of the variable
##'   mixmodCompositeModel(free.proportions=FALSE, variable.independency=TRUE)
##'   # composite models with a pre-defined list
##'   mixmodCompositeModel( listModels=c("Heterogeneous_pk_Ekjh_L_Bk","Heterogeneous_pk_Ekjh_Lk_B") )
##'
##' @author Remi Lebret and Serge Iovleff and Florent Langrognet, with contributions from C. Biernacki and G. Celeux and G. Govaert \email{contact@@mixmod.org}
##' @export
##'
mixmodCompositeModel<- function(listModels=NULL, free.proportions=TRUE, equal.proportions=TRUE, variable.independency=NULL, component.independency=NULL ){
  if ( !is.null(listModels) ){
    new("CompositeModel", listModels=listModels)
  }
  else{
    if ( !is.null(variable.independency) & !is.null(component.independency) ){
      new("CompositeModel", free.proportions=free.proportions, equal.proportions=equal.proportions, variable.independency=variable.independency, component.independency=component.independency)
    }
    else if ( !is.null(variable.independency) & is.null(component.independency) ){
      new("CompositeModel", free.proportions=free.proportions, equal.proportions=equal.proportions, variable.independency=variable.independency)
    }
    else if ( is.null(variable.independency) & !is.null(component.independency) ){
      new("CompositeModel", free.proportions=free.proportions, equal.proportions=equal.proportions, component.independency=component.independency)
    }
    else{
      new("CompositeModel", free.proportions=free.proportions, equal.proportions=equal.proportions)
    }
  }
}

###################################################################################


###################################################################################
##' Get the heterogeneous model name using Gaussian and Multnomial model name
##'
##' @param g_modelname Name of Gaussian model
##' @param m_modelname Name of Multinomial model
##' 
##' @return name of heterogeneous model
##' 
##' @export
##'
composeModelName <- function (g_modelname,m_modelname) {
  gaussianmodel = mixmodGaussianModel(listModels = g_modelname)
  multinomialmodel = mixmodMultinomialModel(listModels = m_modelname)   
  if(gaussianmodel["free.proportions"]!=multinomialmodel["free.proportions"] & gaussianmodel["equal.proportions"]!=multinomialmodel["equal.proportions"])
     stop("Proportions should either be free or equal for both the models.")
  
  if(gaussianmodel["family"]!="diagonal")
     stop("In heterogeneous case, Gaussian model can only belong diagonal family.")
  
  if(gaussianmodel["free.proportions"])
  return(paste("Heterogeneous_",substr(m_modelname,8,nchar(m_modelname)),substr(g_modelname,12,nchar(g_modelname)),sep=""))
  else
    return(paste("Heterogeneous_",substr(m_modelname,8,nchar(m_modelname)),substr(g_modelname,11,nchar(g_modelname)),sep=""))
  
    
}
###################################################################################



###################################################################################
##' @rdname extract-methods
##' @aliases [,CompositeModel-method
##'
setMethod(
  f="[", 
  signature(x = "CompositeModel"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
             "listModels"={return(x@listModels)},
             "free.proportions"={return(x@free.proportions)},
             "equal.proportions"={return(x@equal.proportions)},
             "variable.independency"={return(x@variable.independency)},
             "component.independency"={return(x@component.independency)},
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
##' @aliases [<-,CompositeModel-method
##'
setReplaceMethod(
  f="[", 
  signature(x = "CompositeModel"), 
  definition=function(x,i,j,value){
    if ( missing(j) ){
      switch(EXPR=i,
             "listModels"={x@listModels<-value},
             "free.proportions"={x@free.proportions<-value},
             "equal.proportions"={x@equal.proportions<-value},
             "variable.independency"={x@variable.independency<-value},
             "component.independency"={x@component.independency<-value},
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
