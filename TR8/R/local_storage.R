##' local_storage will download traits data from LEDA and Akhmetzhanova databases
##' and store them in a local folder.
##'
##' Downloading data from the web is time consuming, thus a local storage of some
##' traits data will speed up future data requests; this is possible for LEDA and 
##' Akhmetzhanova databases. The function must be run only once (ideally before running
##' the \code{tr8} function for the first time): thanks to the
##' \code{rappdirs} package, the downloaded data will be stored in the directory
##' commonly used for \code{user data} (which depends on the Operatim System where
##' \code{R} is running). Users can change the destination folder through the \code{directory}
##' parameters, passing the full path of the directory to be used by the function.
##'
##' 
##' @title A utility to storage a local copy of traits data 
##' @param directory the directory where the local Rda files will be stored;
##' default is NULL;
##' @return nothing
##' @author Gionata Bocci <boccigionata@@gmail.com>>
local_storage<-function(db=c("LEDA","Akhmetzhanova","MycoFlor","Catminat","BROT","PLANTS","Ecoflora"),directory){

    ##dir.create(directory,showWarnings=FALSE)

    
    ## download AMF data
    if("Akhmetzhanova"%in%db){
        myco_url <- "http://esapubs.org/Archive/ecol/E093/059/myco_db.csv"
        myco<- tryCatch(read.csv(myco_url,sep=",",header=TRUE),
                        error=function(res){
                            message("URL does not seem to exist:")
                            return(NA)},
                        warning=function(res){
                            message("URL does not seem to exist:")
                            return(NA)
                        })
        save(file=file.path(directory,"myco.Rda"),myco,precheck = F) 
    }

    if("MycoFlor"%in%db){
        mycoflor_url<-"http://www.esapubs.org/archive/ecol/E094/123/MycoFlor.txt"
        MycoFlor<- tryCatch(read.delim(mycoflor_url,header=TRUE),
                            error=function(res){
                                message("URL does not seem to exist:")
                                return(NA)},
                            warning=function(res){
                                message("URL does not seem to exist:")
                                return(NA)
                            } )
        save(file=file.path(directory,"MycoFlor.Rda"),MycoFlor,precheck = F)
    }

     
    if("LEDA"%in%db){
        ## download LEDA data
        leda_download_to_local_directory(directory)
    }

    ## download AMF data
    if("Catminat"%in%db){
        catminat_download_to_local_directory(directory)
    }
    if("BROT"%in%db){
        brot_download_to_local_directory(directory)
    }
    if("PLANTS"%in%db){
        PLANTS_download_to_local_directory(directory)
    }
    if("Ecoflora"%in%db){
        ecoflora_download_to_local_directory(directory)
    }

    

}
