writeDatasetSDML <- function(x,
                             file,
                             textdata,
                             sep,
                             na.string,
                             null.string,
                             posinf.string,
                             neginf.string,
                             nan.string,
                             true,
                             false
                             )
{
    ## 2nd level function: responsible for the correct 
    ## choice of StatDataML objects, currently only lists and arrays

    if (is.null(x)) return(NULL)
    
    catSDML("<dataset>\n", file = file)
    
    writeListArraySDML(x, file = file, textdata = textdata,
                       sep = sep, na.string = na.string,
                       null.string = null.string, posinf.string = posinf.string,
                       neginf.string = neginf.string, nan.string = nan.string,
                       true = true, false = false)
    
    catSDML("</dataset>\n", file = file)
}




