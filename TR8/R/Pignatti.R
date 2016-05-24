## this must be checked when the final Pignatti dataset is provided
setClass("ellenberg_pignatti",representation = list(species_list="vector",results="ANY",not_valid="vector",traits="vector"))

setMethod("initialize",
          "ellenberg_pignatti",
          definition = function(.Object,species_list,results,not_valid,traits,TRAITS){
              .Object@species_list<-species_list

              if(length(TRAITS)>0){
                  .Object@traits<-TRAITS
              }else{
                  .Object@traits<-c("L","T","C","U","R","N","S","life_form_P","corotipo") 
          }
              return(.Object)
          }
          )


setGeneric(name= "get_traits", def=function(.Object){standardGeneric("get_traits")} )

setMethod(f="get_traits",
          signature="ellenberg_pignatti",
          definition = function(.Object){
              
              env<-new.env(parent = parent.frame())
              data(pignatti,envir=env)
              pignatti<-get("pignatti",envir=env)
              DF<-pignatti[pignatti$Name.tnrs%in%.Object@species_list,]
              
              ## the following ddply is needed since pignatti dataframe has
              ## double (or triple) entries for some species [eg. Xanthium strumarium (according
              ## to TNRS name) corresponds to 3 different species in Pignatti dataframe)
              ## thus among the 3 retrieved species I select the one which, according to Pignatti,
              ## is Xanthium strumarium]
              DF<-ddply(DF,.(DF$Name.tnrs),function(x){
                  if(nrow(x)>1){
                      
                      try(return(x[as.character(x$Specie.Pignatti)==as.character(x$Name.tnrs),]),silent=FALSE)
##                      return(x[as.character(x$Specie.Pignatti)==as.character(x$Name.tnrs),])
                  }else{return(x)}
              })

              
##              results<-df[,c("Name.tnrs","forma_biologica","corotipo","L","T","C","U","R","N","S")]
              row.names(DF)<-DF$Name.tnrs
              DF<-DF[,c("L","T","C","U","R","N","S","forma_biologica","corotipo")]
              names(DF)<-c("L","T","C","U","R","N","S","life_form_P","corotipo")
              selected<-.Object@traits[.Object@traits%in%names(DF)]
              results<-DF[,selected]
              
              results<-as.data.frame(results)
              names(results)<-selected#.Object@traits
              row.names(results)<-row.names(DF)
              .Object@results<-results
              .Object@not_valid<-.Object@species_list[!.Object@species_list%in%pignatti$Name.tnrs]
              return(.Object)
          }
          )






#' Extracts ellenberg values for the Italia Flora as provided by Pignatti et al. (2005)
#' 
#' This function is not ment to be used by the final user; it's used by the
#' wrapper function \code{tr8()}
#' 
#' @param species a vector containing plant species names
#' @param TRAITS a vector containing the traits to be downloaded (used as a check for
#' tr8_gui() created variables)
#' @return a data frame
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @references Pignatti S., Menegoni P., Pietrosanti S., 2005, Biondicazione
#' attraverso le piante vascolari. Valori di indicazione secondo Ellenberg
#' (Zeigerwerte) per le specie della Flora d'Italia. Braun-Blanquetia 39,
#' Camerino, pp.  97.
pignatti_f<-function(species,TRAITS){
    res<-new("results")
    ## the second condition is needed when the user
    ## use the tr8_config() and chooses only the beginning
    ## and end of the flowering period: this traits are not
    ## included in the Pignatti databases, but are retrieved
    ## from the web, thus the slot "results" should be
    ## set to null
    if(is.null(TRAITS)){
        res@results<-NULL
    }else{
        if(length(TRAITS)>0&&(!TRAITS%in%c("L","T","C","U","R","N","S","life_form_P","corotipo"))){
            res@results<-NULL
                  }else{

                      obj<-new("ellenberg_pignatti",species_list=species,TRAITS=TRAITS)
                      obj<-get_traits(obj)
#        row.names
                      res@results<-obj@results}
    }
    res@bibliography<-"Pignatti, S., Menegoni, P., Pietrosanti. S., 2005. Biondicazione attraverso le piante vascolari.\nValori di indicazione secondo Ellenberg (Zeigerwerte) per le specie della Flora di Italia. \nBraun-Blanquetia 39, Camerino, pp.  97.\n"
    return(res)
        
}
