##' control: a function to check if the user wants to download some traits
##' from a certain database
##'
##' This function check whether the user has run the \code{tr8_config()}
##' function and, in case he did, which traits were selected (i.e. need
##' to be downloaded by the tr8 function) for each database. These variables
##' have the form "res_NAMEDB" (eg. \code{res_Biolflor}) and they contain
##' the "output" of a "gWidget::notebook" window. The values of these
##' variables can be accessed through the \code{svalue}
##' @title control
##' @param name_variable name of the variable set up by tr8_config()
##' @param dframe a dataframe containing traits definition (created by the tr8() function).
##' @param DB name of the database to be used (eg. "Ecoflora")
##' @seealso tr8()
##' @return a vector of seleceted traits (if the variable
##' was set through the \code{tr8_config()} function OR \code{NULL}
##' if \code{tr8_config()} was run, but no traits were chosen for that
##' database OR an empty vector if  \code{tr8_config()} was not run.
##' @author Gionata Bocci <boccigionata@gmail.com>
control<-function(name_variable,dframe,DB){
    
    if(exists(name_variable)){
        ## get the values of the subset
        trait_list<-eval(parse(text=paste("svalue(",name_variable,")",sep="")))
        ## test if the subset contains traits or not
        if(length(trait_list)>0){
            ## extract the short_code for the traits
            df_list<-dframe$long_code[dframe$description%in%trait_list&dframe$db%in%DB]
        }else
            {
                ## if tr8 was not configured, set the variable, which will be 
                ## passed as an argument, to NULL
                df_list<-NULL    
            }
    }else{
        df_list<-character()
                    }
    return(df_list)
}
