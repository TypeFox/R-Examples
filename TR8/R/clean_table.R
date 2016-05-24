#' Removes redundant text in traits collected from Biolflor
#' 
#' BiolFlor tables contains brief explanations of traits: that is ok for the
#' website but tends to produce clumsy tables in dataframes returned by tr8(),
#' thus this extra-text is removed (to improve readibility of such tables).
#'  This function is internally used by \code{tr8()}, users do not need to run
#' it.
#' 
#' @param input a intermediate dataframe retrieved by \code{tr8()} 
#' @return a dataframe with shortened names for the traits levels' values
#' @author Gionata Bocci <boccigionata@@gmail.com>
#' @seealso \code{\link{biolflor}}
#' @keywords traits
#' @examples \dontrun{
#' bolflor_clean(biolflor("Avena sativa"))
#' }
biolflor_clean<-function(input){
    DF<-input


    for(param in c("Life.span","Life.form","Rosettes","Vegetative.propagation","Storage.organs","Type.of.reproduction","Strategy.type","Pollen.vector","Flower.class.after.MUELLER")){
        ## only clean these columns if they are really
        ## present in the retrieved data
        if(param %in% names(DF)){
            ##extra text is included between brackets -> thus here text enclosed in brackets is removed
            eval(substitute((DF$AA<-gsub("\\(.*\\)","",DF$AA)),list(AA=param)))
            ## remove backtick when present
            eval(substitute((DF$AA<-gsub("`","",DF$AA)),list(AA=param)))
            ## for pollinators there's also text which begins with "Typical pollinators..." -> that is removed as well
            eval(substitute((DF$AA<-gsub("Typical pollinators.*","",DF$AA)),list(AA=param)))
        }
       }
    return(DF)
}


