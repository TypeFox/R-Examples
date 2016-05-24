##' Shortens up columns' names in the intermidiate steps of
##' tr8() retrieved data
##'
##' This function is internally used by \code{tr8()}.
##' @title column_conversion
##' @return a dataframe with shortened col names
##' @param DF a intermediate dataframe built by \code{tr8} 
##' @author Gionata Bocci <boccigionata@@gmail.com>
column_conversion<-function(DF){

    env<-new.env(parent = parent.frame())
    data(column_list,envir=env)
    column_list<-get("column_list",envir=env)
    for(i in names(column_list)){
        if(i%in%names(DF)){
            names(DF)[which(names(DF)==i)]<-column_list[i][[1]][1]
        }
    }
    remove(list=c("column_list"), envir=env)
    return(DF)
}
    
