
## list_of_traits_Biolflor=c(
##     "Life form"                        
##     ,"Life span"                              
##     ,"Rosettes"                               
##     ,"Vegetative propagation"                 
##     ,"Storage organs"                         
##     ,"Type of reproduction"                   
##     ,"Type of seed production"                
##     ,"Strategy type"                          
##     ,"Pollen vector"
##     )

traits_special_Biolflor=c(
    "Flower class after MUELLER",
    "Begin of flowering (months)",
    "End of flowering (months)",
    "Duration of flowering (months)",
    "Number of flowering phases"
)

traits_pollen_Biolflor=c("Pollen vector")


## biolflor_traits a class for containing traits retrieved from Biolflor
##
## This class is used to retrieve traits data from species contained
## in the Biolflor web database
##
## @slot url url of the corresponding web page for the species of interest
## @slot extracted a list of traits
## @slot list_traits a vector containing the list of traits which can be downloaded by the package
## @slot list_special_traits some of the traits require special Xpath rules to be extracted
## @exportClass biolflor_traits
setClass("biolflor_traits",representation=list(url="character",extracted="list",list_traits="vector",list_special_traits="vector",list_pollen_traits="vector",results="data.frame"))


## extract data fro biolflor_traits classes
## @param .Object an object of biolflor_traits signature
setGeneric(name="extract",def=function(.Object){standardGeneric("extract")})

setMethod(f="extract",
          signature="biolflor_traits",
          definition=function(.Object){


              if(.Object@url=="not present"){
                  for(trait in .Object@list_traits){
                      .Object@extracted[[trait]]=NA
                  }
                  for(trait in .Object@list_special_traits){
                      .Object@extracted[[trait]]=NA
                  }
                  for(trait in .Object@list_pollen_traits){
                      .Object@extracted[[trait]]=NA
                      }
              }

              else{
                  ##set language to english
                  form="http://www2.ufz.de/biolflor/index.jsp"
                  ##base_url<-"http://www2.ufz.de"
                  param<-list("language"='en')
                  vai<-getForm(form,.params=list("language"="en"),style="POST")
                  
                  ##get web page
                  base_url<-"http://www2.ufz.de"
                  temp_pag<-htmlParse(getURL(paste(base_url,.Object@url,sep="")))
                  
                  ## parse html and extract data

                  for(trait in .Object@list_traits){
                      query=paste("//*[text()='",trait,"']/following-sibling::td/a",sep="")
                      value=xpathApply(temp_pag,query,xmlValue)

                      if(length(value) > 0) {
                          ##                      print(value)
                          .Object@extracted[[trait]]=paste(unlist(value),collapse = " - ")
                      }else{.Object@extracted[[trait]]=NA}
                  }

                  ## extract special traits (pollen)
                  for(trait in .Object@list_special_traits){
                      query=paste("//*[text()='",trait,"']/following-sibling::td",sep="")
                      value=xpathApply(temp_pag,query,xmlValue)
                      if(length(value) > 0) {
                          ##                      print(value)
                          .Object@extracted[[trait]]=paste(unlist(value),collapse = " - ")
                      }
                      else{.Object@extracted[[trait]]=NA}
                      
                  }



                  for(trait in .Object@list_pollen_traits){
                      query1=paste("//*[text()='",trait,"']/following-sibling::td/a",sep="")
                      value1=xpathApply(temp_pag,query1,xmlValue)


                      query2=paste("//*[text()='",trait,"']/ancestor::*[1]/following-sibling::tr[1]/*[text()='Abundance']/following-sibling::td/text()",sep="")
                      value2=xpathApply(temp_pag,query2,xmlValue)
                      value2<-lapply(value2,function(x){gsub(" \\[.*\\]","",x)})
                      value2<-lapply(value2,function(x){gsub("(.*)","[\\1]",x)})


                      if(length(value2) > 0) {
                          ##                      print(value)
                          .Object@extracted[[trait]]=paste(paste(value1,value2),collapse=" - ")
                      }
                      else{.Object@extracted[[trait]]=NA}
                      
                  }



              }
              return(.Object)

          }
          )



#' Retrieve traits data from the BiolFlor website.
#' 
#' This function allows the user to download some pre-defined traits from the
#' BiolFlor website: the function returs a dataframe with species name in row and traits data in
#' column.
#' 
#' @param list_species vector containing names of those plant species for
#' which traits data need to be downloaded.
#' @param TRAITS a vector containing the traits to be downloaded (used as a check for tr8_gui() created variables)
#' 
#' @return dataframe with species name in row and traits data in
#' column.
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @references BIOLFLOR - Eine Datenbank zu biologisch-ökologischen Merkmalen zur Flora von Deutschland.
#' Schriftenreihe für Vegetationskunde 38: 1-333.  (Bundesamt für. Bonn, Bundesamt für Naturschutz)
#' @keywords traits
#' @usage biolflor(list_species,TRAITS)
#' @examples \dontrun{
#' biolflor(c("Abies alba"))
#' }
biolflor<-function(list_species,TRAITS,rest=NULL){
    res<-new("results")
    env<-new.env(parent = parent.frame())
    data(biolflor_lookup,envir=env)
    biolflor_lookup<-get("biolflor_lookup",envir = env)    
    if(is.null(TRAITS)){
        res@results<-NULL
    }else{
        ## if(tryCatch(nsl("www.cran.r-project.org"), error =function(e){return(FALSE)},warning=function(w){return(FALSE)})==FALSE){
        ##     stop("You do not have a working internet connection.\n  Please re-run tr8() function when your internet connection is working.")
        ## }
        if(length(TRAITS)>0){
            ## otherwise check which of the selected traits are available
            list_of_traits_Biolflor<-list_of_traits_Biolflor[list_of_traits_Biolflor%in%TRAITS]
            ## check also the "special biolflor traits"
            traits_special_Biolflor<-traits_special_Biolflor[traits_special_Biolflor%in%TRAITS]
            traits_pollen_Biolflor<-traits_pollen_Biolflor[traits_pollen_Biolflor%in%TRAITS]
            ## use the user's selected traits to instantiate the class
            ## and retrieve data
            
        }

        tmp_list=list()
        for(cur in list_species){
            if(cur%in%biolflor_lookup$submittedname){
                species_url<-with(biolflor_lookup,biolflor_lookup[submittedname==cur,"V2"])
                ## species_url<-with(biolflor_lookup,biolflor_lookup[acceptedname==cur|submittedname==cur,"V2"])
                ## ## check if 2 species have the same accepted name
                ## if(length(species_url)>1){
                ##     species_url<-with(biolflor_lookup,biolflor_lookup[acceptedname==cur&submittedname==cur,"V2"])
                ## }
            }else{species_url<-"not present"}
            Sys.sleep(rest)
            prova<-new("biolflor_traits",url=species_url,list_traits=list_of_traits_Biolflor,list_special_traits=traits_special_Biolflor,list_pollen_traits=traits_pollen_Biolflor)
            bio_res<-extract(prova)
            go<-bio_res@extracted
            ##            go<-data.frame(go,check_names=F)           
            for(i in names(go)){
                tmp_list[[cur]][[i]]<-go[i][[1]]

            }


        }
        tp<-ldply(tmp_list)
        row.names(tp)<-tp$.id
        names_species<-names(tp)[names(tp)!=".id"]
        ## drop=FALSE is necessary to avoid that single column dataframe (ie only one
        ## trait is selected from biolflor) is converted to a vector (without names
        ## and row.names
        tp<-data.frame(tp[,!names(tp)==".id",drop=FALSE],check.names = F)
        
        res@results<-tp

        }
    remove(list=c("biolflor_lookup"), envir = env)    
    stringa<-"Klotz, S., K\303\274hn, I., Durka, W. (eds), 2002. BIOLFLOR - Eine Datenbank zu \nbiologisch-\303\266kologischen Merkmalen zur Flora von Deutschland. Schriftenreihe\nf\303\274r Vegetationskunde 38: 1-333. (Bundesamt f\303\274r. Bonn, Bundesamt f\303\274r Naturschutz).\n"
    Encoding(stringa)<-"unicode"
    res@bibliography<-stringa
    
    return(res)
}

