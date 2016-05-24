##' The function is used to create the available_traits data frame which is meant to help users
##' in showing which traits (and from which databases) can be downloaded
##'
##' The function shows what data are available for download
##' and decide which one should be passed to the \code{tr8()} function (in the
##'  \code{download_list} argument); the codes to be used as the \code{download_list} argument
##' are those contained \code{short_code} column.
##' 
##' @title available_traits shows which traits are available for download
##' @return a data frame
##' @seealso tr8
##' @author Gionata Bocci <boccigionata@@gmail.com>
##' @examples
##' available_traits()
##' ## If the traits \code{Maximum area}  and \code{Leaf area} from
##' ## Ecoflora are needed for the species Salix alba and Populus nigra, type
##' 
##' \dontrun{tr8(species_list=c("Salix alba","Populus nigra"),download_list=c("h_max","le_area"))}
available_traits<-function(){

    env<-new.env(parent = parent.frame())
    data(column_list,envir = env)
    column_list<-get("column_list",envir=env)
    ## load lookup table
    ## convert it to a data frame
    temp_dframe<-ldply(column_list)
    names(temp_dframe)<-c("long_code","short_code","description","db")
    remove(list=c("column_list"), envir = env)    
    return(temp_dframe[,c("short_code","description","db")])

}

