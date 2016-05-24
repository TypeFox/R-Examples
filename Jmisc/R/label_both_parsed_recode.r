#' Combine label_both and label_parsed in \pkg{ggplot2}. Also added a rename function to it
#' see label_both and label_parsed in \pkg{ggplot2} for details.
#' @name label_both_parsed_recode
#' @aliases label_both_parsed_recode
#' @title Combine label_both and label_parsed in \pkg{ggplot2}.
#' @param display_name A vector contains the display name. Names of the vector are the original name.
#' @return A function similar to label_both and label_parsed in \pkg{ggplot2} for details.
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @references \url{http://cran.r-project.org/web/packages/ggplot2/index.html}
#' @export
label_both_parsed_recode <-
function(display_name){
    function(variable, value){
        variable=recode(variable,from=names(display_name),to=display_name)
        v<-paste(variable, value, sep = "==")
        lapply(as.character(v), function(x) parse(text = x))
    }
}
