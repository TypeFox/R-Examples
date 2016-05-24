### * Definition of Ecoflora Class

## Class containing traits from Ecological Flora of the british Isles
##
## The class is a box containing url and traits data for species
## which are present in the Ecoflora website
setClass("Ecoflora",representation=list(species_list="vector",reference="data.frame",df="data.frame",not_valid="vector",results="data.frame",traits="list",double_names="vector",rest="numeric",issues="ANY"))



## * initialization of Ecoflora Class

## During the initialization of an \code{ecoflora} object, the
## method will take care of chacking and correcting url of
## plant species

setMethod('initialize',
          signature="Ecoflora",
          definition = function(.Object,species_list,reference,traits,rest){          
              .Object@species_list<-species_list
              .Object@rest<-rest
              ## list which will contain ecoflora
              ## web link for each species
              base_url<-"http://www.ecoflora.co.uk/search_ecochars.php?"
              lookup=list()
              for(i in species_list){
                  ## is the species of interest either found in ecoflora matchedname, species name or acceptedname?
                  ##results<-with(reference,reference[matchedname==i|species==i|acceptedname==i,])
                  results<-with(reference,reference[species==i,])
                  ## if there's just a single species, everything is fine, just extract the corresponding url
                  if(nrow(results)==1){
                      lookup[[i]]<-paste(base_url,results$web_link,sep="")
                  }
                  else
                      {
                          ## if multiple species are found, just take the url for that species which has acceptedname==ecoflora name
                          if(nrow(results)>1){

                              results2<-with(reference,reference[species==i,])

                              ## in some cases Ecoflora has two species names (not accepted) which
                              ## correspond to a single accepted name eg. but are both different from that
                              ##
                              ##           species           web_link      acceptedname     score     matchedname
                              ## Bromopsis benekenii  plant_no=1930440110  Bromus ramosus     1    Bromopsis benekenii
                              ## Bromopsis ramosa     plant_no=1930440100  Bromus ramosus     1     Bromopsis ramosa
                              ##
                              ## in those cases results2 will have 0 rows => use back'results' and choose the first occurence
                              ## plus put a  warning message in .Object@issues
                              
                              if(nrow(results2)==0){
                                  stringa<-"\n\tBEWARE: Ecoflora contains the following species:\n\t - "
                                  stringa1<-paste(results$species,collapse="\n\t - ")
                                  stringa2<-"for which multiple entries are present. Only one will be included; please check results on Ecolflora web site\n"
                                  ## stringa2<-paste("\n\twhich correspond to the same accepted name: ",unique(results$acceptedname))
                                  ## stringa3<-paste("\n\tIn this case only data corresponding to Ecoflora species",results$species[1],"will be used.\n",sep=" ")
                                  ## stringa4<-paste("\n\tYou may want to re-run tr8 using the other Ecoflora species names to get traits for them\n")
                                  ## tot_alert<-paste(stringa,stringa1,stringa2,stringa3,stringa4)
                                  tot_alert<-paste(stringa,stringa1,stringa2)
                                  
                                  .Object@issues<-c(.Object@issues,tot_alert)
                                  results2<-results[1,]
                              }


                              lookup[[i]]<-paste(base_url,results2$web_link,sep="")
                          }
                          else
                              {
                                  ## if the above doesn't work, put the species in the "non valid" class slot
                                  .Object@not_valid<-c(.Object@not_valid,i)
                                  lookup[[i]]<-"not found"
                              }
                      }
              }
              ## convert the list to a dataframe and store the result in the df slot
              lookup<-ldply(lookup)
              if(nrow(lookup)>0){
                  names(lookup)<-c("species","web_link")
                  .Object@df<-lookup
              }
              .Object@traits<-traits

              return(.Object)
              }

          )
              
### * 'retrieve' method for "Ecoflora"

## The 'retrieve' method for "Ecoflora" objects will download data from
## the Ecoflora website for the passed species
setGeneric( name= "retrieve", def=function(.Object){standardGeneric("retrieve")} )


setMethod(f='retrieve',
          signature='Ecoflora',
          definition = function(.Object){
                  ## eco will contain as many slots as the species passed
                  ## and each slot-species will contain a list of the
                  ## downloaded ecological traits
                  eco<-list()

                  
                  for(species in .Object@df$species){
                      ## url of the web page for the species of interest
                      species_url<-.Object@df$web_link[.Object@df$species==species]
                      if(species_url=="not found"){
                          for(trait in names(.Object@traits)){
                                  eco[[species]][trait]<-NA
                          }
                      }
                      else{
                          ## extract tabe data from the scraped web page
                          Sys.sleep(.Object@rest)
                          eco_data<-readHTMLTable(species_url)[[2]]
                          ##eco_data<-readHTMLTable(species_url)
                          ## for some traits there are several entries (with the same Code), thus
                          ## the retrieved table must be "aggregated" in order to have 1 entr/trait
                          eco_data<-aggregate(eco_data$Value,by=list(eco_data$Number),paste,collapse=';')                  
                          names(eco_data)<-c("Code","Value")
                          
                          ## fill in the list "eco"; NA values are used for those traits
                          ## which do not have values in the Ecoflora database
                          for(trait in names(.Object@traits)){
                              if(.Object@traits[trait]%in%eco_data$Code){
                                  eco[[species]][trait]<-eco_data$Value[eco_data$Code==.Object@traits[trait]]
                              }else{
                                  eco[[species]][trait]<-NA}
                          }
                      }
                  }
                  ## a dataframe is returned
                  ## NB: species names, being
                  
                  eco<-t(as.data.frame(eco))
                  row.names(eco)<-.Object@df$species
                  
                  ##              return(eco)
                  eco<-as.data.frame(eco)
                  .Object@results<-eco
              return(.Object)
              }
          )




### * wrapper function to extract data from Ecoflora

#' Retrieves traits data from Ecoflora website
#' 
#' The function accepts a list of plant species names, tries to download the
#' corresponding functional traits from the Ecoflora website
#' (\samp{http://www.ecoflora.co.uk/}) and return a data.frame with species
#' names as rows and functional traits as columns.
#' 
#' @param species_list a vector containing list of plant species names.
#' @param reference the reference lookup data.frame (this is not ment to be set
#' by users; it is left here for further development)
#' @param TRAITS a vector containing the traits to be downloaded (used as a check for
#' tr8_gui() created variables)
#' @return Return a data.frame with species as rows and traits as columns.
#' Only those species present in the Ecoflora database will be included in this
#' data.frame, other species will be left out.
#' @author Bocci Gionata
#' @references Fitter, A . H. and Peat , H. J., 1994, The Ecological Flora Database, J. Ecol., 82, 415-425.
#' @seealso \code{\link{traits_eco}}
#' @examples \dontrun{
#' #' #My_data<-ecoflora(species_list=c("Abies alba"))
#' }
ecoflora<-function(species_list,TRAITS,rest)
    {
        env<-new.env(parent = parent.frame()) ##        env<-new.env()
        res<-new("results")
        ##data("ECOFLORA_df",envir=env)
        ##data("traits_eco",envir=env)
        ECOFLORA_df<-get("ECOFLORA_df",envir=env)
        traits_eco<-get("traits_eco",envir=env)
        ## test if Ecoflora is providing data (if not the web page
        ##  http://www.ecoflora.co.uk/search_species.php will contain
        ## "No Species currently available"
        #eco_check<-readLines("http://www.ecoflora.co.uk/search_species.php",warn=FALSE)
        #res_check<-length(grep("No Species currently available",eco_check))
        ## if Ecoflora is not working, res_check will be equal to 1
        

        ## if traits is NULL it means that the user did not selected
        ## a subset of traits (by means of the tr8_config function, thus
        ## all the traits should be downloaded
        
        
        if(is.null(TRAITS)){
            res@results<-NULL
        }else{
            ## check that internet connection is working
            ## otherwise it will stop and provide an error 
            ## if(tryCatch(nsl("www.cran.r-project.org"), error =function(e){return(FALSE)},warning=function(w){return(FALSE)})==FALSE){
            ##     stop("You need a working internet connection to download traits from Ecolflora")
            ## }

            ## check that ecoflora is up and working
            
            eco_check<-readLines("http://www.ecoflora.co.uk/search_species.php",warn=FALSE)
            res_check<-length(grep("No Species currently available",eco_check))
            
            
            if(length(TRAITS)==0){
                traits<-traits_eco}else{
                    
                    ## if a subset has been passed, only the
                    ## corresponding codes should be used
                    traits<-traits_eco[names(traits_eco)%in%TRAITS]
                }
            if(res_check==1){
                message("Ecoflora is not working at the moment, please retry later.")
                res@results<-NULL
            }else{
            
                obj<-new("Ecoflora",species_list=species_list,reference=ECOFLORA_df,traits=traits,rest=rest)
            ##        ret<-as.data.frame(ret@results)
            ##remove(list=c("ECOFLORA_df","traits_eco"),pos =".GlobalEnv")
                ret<-retrieve(obj)
                res@results<-ret@results
                res@issues<-ret@issues
            }
        }
        ##remove(list=c("ECOFLORA_df","traits_eco"), envir = env)
        res@bibliography<-"Fitter, A . H. and Peat , H. J., 1994. The Ecological Flora Database,\nJ. Ecol., 82, 415-425.  http://www.ecoflora.co.uk"

        return(res)

    }
