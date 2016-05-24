###################################################################################
##                             MultinomialModel.R                                ##
###################################################################################

###################################################################################
##' @include global.R
##' @include Model.R
NULL
###################################################################################

###################################################################################
##' Constructor of [\code{\linkS4class{MultinomialModel}}] class
##'
##' This class defines a multinomial Model. Inherits the [\code{\linkS4class{Model}}] class.
##' 
##' \describe{
##'   \item{variable.independency}{logical}
##'   \item{component.independency}{logical}
##' }
##'
##' @examples
##'   new("MultinomialModel")
##'   new("MultinomialModel", listModels=c("Binary_pk_E","Binary_p_E") )
##'   new("MultinomialModel", free.proportions=FALSE, variable.independency=TRUE )
##'
##'   getSlots("MultinomialModel")
##' 
##' @name MultinomialModel-class
##' @rdname MultinomialModel-class
##' @exportClass MultinomialModel
##'
setClass(
    Class="MultinomialModel",
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
      vcf<-"Binary_pk_E"
      vce<-"Binary_p_E"
      f<-c("Binary_pk_Ekj","Binary_pk_Ekjh")
      e<-c("Binary_p_Ekj","Binary_p_Ekjh")
      cf<-"Binary_pk_Ej"
      ce<-"Binary_p_Ej"
      vf<-"Binary_pk_Ek"
      ve<-"Binary_p_Ek"
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
        stop("At least one model is not a valid model. See ?mixmodMultinomialModel for the list of all multinomial models.")

      # check proportions parameters validity
      if ( !object@equal.proportions & !object@free.proportions )
        stop("equal.proportions and free.porportions cannot be both as FALSE !")
      
      if ( !object@free.proportions & (sum(object@listModels %in% all.free)>0) )
        stop("At least one model has a free proportions but free.proportions is set as FALSE. See ?mixmodMultinomialModel for the list of models with equal proportions.")
      
      if ( !object@equal.proportions & (sum(object@listModels %in% all.equal)>0) )
        stop("At least one model has an equal proportions but equal.proportions is set as FALSE. See ?mixmodMultinomialModel for the list of models with free proportions.")

      # check independencies parameters
      # for variable
      if ( length(object@variable.independency) ){
        if ( object@variable.independency & sum(object@listModels %in% variable) != length(object@listModels) )
          stop("At least one model is not independent of the variable j. See ?mixmodMultinomialModel for the list of all multinomial models.")
      }
      # for component
      if ( length(object@component.independency) ){
        if ( object@component.independency & sum(object@listModels %in% component) != length(object@listModels) )
          stop("At least one model is not independent of the variable j. See ?mixmodMultinomialModel for the list of all multinomial models.")
      }
      
    }
)
###################################################################################


###################################################################################
##' Create an instance of the [\code{\linkS4class{MultinomialModel}}] class using new/initialize.
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
  signature=c("MultinomialModel"),
  definition=function(.Object, listModels, free.proportions, equal.proportions, variable.independency, component.independency){
    
    # define list of models
    vcf<-"Binary_pk_E"
    vce<-"Binary_p_E"
    f<-c("Binary_pk_Ekj","Binary_pk_Ekjh")
    e<-c("Binary_p_Ekj","Binary_p_Ekjh")
    cf<-"Binary_pk_Ej"
    ce<-"Binary_p_Ej"
    vf<-"Binary_pk_Ek"
    ve<-"Binary_p_Ek"
    
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
##' Create an instance of the [\code{\linkS4class{MultinomialModel}}] class
##'
##' Define a list of multinomial model to test in MIXMOD.
##'
##' In the multinomial mixture model, the multinomial distribution is associated to the \eqn{j}th variable of the \eqn{k}th component is reparameterized by a center \eqn{a_k^j} and the dispersion \eqn{\varepsilon_k^j} around this center. Thus, it allows us to give an interpretation similar to the center and the variance matrix used for continuous data in the Gaussian mixture context. In the following, this model will be denoted by \eqn{[\varepsilon_k^j]}. In this context, three other models can be easily deduced. We note \eqn{[\varepsilon_k]} the model where \eqn{\varepsilon_k^j} is independent of the variable \eqn{j}, \eqn{[\varepsilon^j]} the model where \eqn{\varepsilon_k^j} is independent of the component \eqn{k} and, finally, \eqn{[\varepsilon]} the model where \eqn{\varepsilon_k^j} is independent of both the variable $j$ and the component \eqn{k}.  In order to maintain some unity in the notation, we will denote also \eqn{[\varepsilon_k^{jh}]} the most general model introduced at the previous section.
##' 
##' @param listModels a list of characters containing a list of models. It is optional.
##' @param free.proportions logical to include models with free proportions. Default is TRUE.
##' @param equal.proportions logical to include models with equal proportions. Default is FALSE.
##' @param variable.independency logical to include models where \eqn{[\varepsilon_k^j]} is independent of the variable \eqn{j}. Optionnal.
##' @param component.independency logical to include models where \eqn{[\varepsilon_k^j]} is independent of the component \eqn{k}. Optionnal.
##'
##' @return an object of [\code{\linkS4class{MultinomialModel}}] containing some of the 10 Binary Models:
##' \tabular{rlll}{
##'     Model \tab Prop. \tab Var. \tab Comp. \cr
##'     Binary_p_E     \tab Equal \tab TRUE \tab TRUE \cr
##'     Binary_p_Ej    \tab \tab FALSE \tab TRUE \cr
##'     Binary_p_Ek    \tab \tab TRUE \tab FALSE \cr
##'     Binary_p_Ekj   \tab \tab FALSE \tab FALSE \cr
##'     Binary_p_Ekjh  \tab \tab FALSE \tab FALSE \cr
##'     Binary_pk_E    \tab  Free \tab TRUE \tab TRUE \cr
##'     Binary_pk_Ej   \tab \tab FALSE \tab TRUE \cr
##'     Binary_pk_Ek   \tab \tab TRUE \tab FALSE \cr
##'     Binary_pk_Ekj  \tab \tab FALSE \tab FALSE  \cr
##'     Binary_pk_Ekjh \tab \tab FALSE \tab FALSE \cr
##' }
##'
##' @references 
##'   R. Lebret, S. Iovleff, F. Langrognet, C. Biernacki, G. Celeux, G. Govaert (2015), "Rmixmod: The R Package of the Model-Based Unsupervised, Supervised, and Semi-Supervised Classification Mixmod Library", Journal of Statistical Software, 67(6), 1-29, doi:10.18637/jss.v067.i06
##' @examples
##'   mixmodMultinomialModel()
##'   # multinomial models with equal proportions
##'   mixmodMultinomialModel(equal.proportions=TRUE,free.proportions=FALSE)
##'   # multinomial models with a pre-defined list
##'   mixmodMultinomialModel( listModels=c("Binary_pk_E","Binary_p_E") )
##'   # multinomial models with equal proportions and independent of the variable
##'   mixmodMultinomialModel(free.proportions=FALSE, variable.independency=TRUE)
##'
##' @author Remi Lebret and Serge Iovleff and Florent Langrognet, with contributions from C. Biernacki and G. Celeux and G. Govaert \email{contact@@mixmod.org}
##' @export
##'
mixmodMultinomialModel<- function( listModels=NULL, free.proportions=TRUE, equal.proportions=TRUE, variable.independency=NULL, component.independency=NULL ){

    if ( !is.null(listModels) ){
      new("MultinomialModel", listModels=listModels)
    }else{
      if ( !is.null(variable.independency) & !is.null(component.independency) ){
        new("MultinomialModel", free.proportions=free.proportions, equal.proportions=equal.proportions, variable.independency=variable.independency, component.independency=component.independency)
      }
      else if ( !is.null(variable.independency) & is.null(component.independency) ){
        new("MultinomialModel", free.proportions=free.proportions, equal.proportions=equal.proportions, variable.independency=variable.independency)
      }
      else if ( is.null(variable.independency) & !is.null(component.independency) ){
        new("MultinomialModel", free.proportions=free.proportions, equal.proportions=equal.proportions, component.independency=component.independency)
      }
      else{
        new("MultinomialModel", free.proportions=free.proportions, equal.proportions=equal.proportions)
      }
    }
}
###################################################################################


###################################################################################
##' @rdname extract-methods
##' @aliases [,MultinomialModel-method
##'
setMethod(
  f="[", 
  signature(x = "MultinomialModel"),
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
##' @aliases [<-,MultinomialModel-method
##'
setReplaceMethod(
  f="[", 
  signature(x = "MultinomialModel"), 
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

