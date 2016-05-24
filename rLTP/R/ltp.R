#' R Interface of LTP-cloud service
#'
#' This function deals with communication with the server.
#' XML result will be parsed if the mission is word-splitting.
#' Else the raw XML texts will be returned for further analysis.
#' 
#' @param input The input text.
#' @param file The input file.
#' @param mission Expected result for the cloud server, may be unfinished. Optional choices are 
#' 'ws' for word-splitting,
#' 'pos' for part-of-speech,
#' 'ner' for named entity recognition,
#' 'dp' for dependency parser,
#' 'srl' for semantic role labeling,
#' 'all' for all missions.
#' @param api_key Your API_Key for the cloud server. Visit http://www.ltp-cloud.com/dashboard/ to get it.
#' @param maxUpload Due to the limitation of the server, we cut the input in pieces.
#' @examples
#' require(rLTP)
#' # This api_key is publicly accessible.
#' # So it is strongly recommended to register for your own key.
#' options(ltp_api_key='l2T9N724koSqEcDJvQHtRGVV2erajgPOgB0FAcLj')
#' ltp('Replace this field with a Chinese sentence.')
#' 
#' @export
#' 
ltp = function(input=NULL,file=NULL,mission='ws',
               api_key = getOption('ltp_api_key'),
               maxUpload=100000)
{
    if (is.null(input) && is.null(file))
        stop('No Input.')
    if (!is.null(file))
        input = readLines(file)
    if (!isUTF8(input[1]))
        input = toUTF8(input)
    input = gsub('\\s+','',input)
    input = input[input!='']
    if (length(input)==0)
        stop('Empty Input.')
    n = length(input)
    if (n==1)
    {
        if (nchar(input)>maxUpload)
            stop('This paragraph is too long to upload, please split it into smaller pieces.')
        inputs = input
    }
    else
    {
        nc = nchar(input[1])
        inputs = input[1]
        L = 1
        for (i in 2:n)
        {
            if (nchar(inputs[L])+nchar(input[i])<maxUpload)
                inputs[L] = paste(inputs[L],input[i],sep='\n')
            else
            {
                L = L+1
                inputs[L] = input[i]
            }
        }
    }
    
    if (length(inputs)==1)
    {
        cat('Uploading\n')
        para = commLTP(inputs,mission=mission,api_key)
        if (mission!='ws')
        {
            cat('Parser for non-Word-Splitting result is not available yet.',
                '\n','Returning Raw json texts for each part.\n')
            return(para)
        }
        ans = para
        cat('Analysis Finished!\n')
    }
    else
    {
        cat('Uploading\n')
        paras = list()
        n = length(inputs)
        p = proc.time()
        for (i in 1:n)
        {
            paras[[i]] = commLTP(inputs[i],mission,api_key)
            message(proc.time()-p)
            cat(paste(i,n,sep='/'),'\n')
            p = proc.time()
        }
        if (mission!='ws')
        {
            cat('Parser for non-Word-Splitting result is not available yet.',
                '\n','Returning Raw json texts for each part.\n')
            return(paras)
        }
        else
        {
            ans = unlist(paras)
            cat('Analysis Finished!\n')
        }
        
    }
    ans
}
