#' R Interface of PullWord service
#'
#' This function deals with communication with the server.
#' The result will be parsed in a vector or a matrix, depends on show probability or not
#' 
#' @param input The input text.
#' @param file The input file.
#' @param threshold The minimum probability for the words appearing in the result. 
#'   Should be a real value between 0 and 1.
#' @param showProb logical. The return value would be a \code{data.frame} if \code{TRUE}, or a \code{vector} otherwise.
#' @examples
#' require(pullword)
#' pullword('Replace this field with a Chinese sentence.',threshold=0,showProb=TRUE)
#' 
#' @export
#' 
pullword = function(input=NULL,file=NULL,threshold=0,showProb=FALSE)
{
    if (is.null(input) && is.null(file))
        stop('No Input.')
    if (!is.null(file))
        input = readLines(file)
    if (!isUTF8(input[1]))
        input = toUTF8(input)
    if (length(input)==0)
        stop('Empty Input.')
    n = length(input)
    if (n>1)
        input = paste(input,sep='',collapse='')
    
    threshold = max(min(threshold,1),0)
    
    # Upload 
    cat('Uploading\n')
    result = getForm(uri = "http://api.pullword.com/get.php?",
                     style = "get",
                     'source' = input,
                     'param1' = threshold,
                     'param2' = as.numeric(showProb))
    # Parse
    if (showProb)
    {
        result = strsplit(result, '\r\n')[[1]]
        result = result[result!='']
        result = strsplit(result, ':')
        result = do.call(rbind,result)
        ans = data.frame(word = result[,1], prob = as.numeric(result[,2]), 
                         stringsAsFactors = F)
    } 
    else
    {
        result = strsplit(result, '\r\n')[[1]]
        result = result[result!='']
        ans = result
    }
    
    ans
}
