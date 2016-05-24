#' Generate a design matrix from raw RNAi data.
#' 
#' Function takes the raw RNAi data as input and generates a design matrix for regression.
#'
#' @param data Matrix of RNAi measurements; includes columns batch, query_gene and template_gene
#' @param covariates Vector of strings; each string is the name of a covariate
#' @return a design matrix (without the intercept).
#' @export
#'

constructDesignMatrix <- function(data, covariates){
    # initialize design matrix
    nMeasurements <- nrow (data)
    nCovariates <- length (covariates)
    designMatrix <- matrix(0, nrow=nMeasurements, ncol=nCovariates)
    colnames(designMatrix) <- covariates
    
    # assign 1 to corresponding batch
    designMatrix[,"batch1"] <- data$batch=="batch1"
    designMatrix[,"batch2"] <- data$batch=="batch2"
    designMatrix[,"batch4"] <- data$batch=="batch4"
    
    # go through each row of data
    # and assign 1 to corresponding covariate(s)
    for (i in 1:nMeasurements){
        template <- as.character(data$template_gene[i])
        query <- as.character(data$query_gene[i])
        
        # main effect of template
        # if template gene not empty
        # then assign 1 to corresponding covariate
        if (!((template=="empty")|(template=="NT"))){
            designMatrix[i,template] <- 1
        }
        
        # main effect of query
        # if query gene not empty
        # then assign 1 to corresponding covariate
        if (!(query=="empty")|(query=="NT")){
            designMatrix[i,query] <- 1
        }
        
        # interaction between template and query
        # if neither template or query is empty
        # then assign 1 to corresponding pair
        if ((template!="empty") & (template!="NT") & (query!="empty") & (query!="NT")){
            # if template same as query
            # then treat as single knockdown
            if (template == query) {
                designMatrix[i,template] <- 1
            } else {
                # otherwise
                # add 1 to corresponding gene pair
                idc <- ifelse (template<query, paste (template, query, sep=":"), paste (query, template, sep=":"))
                #print (idc)
                designMatrix[i,idc] <- 1
            }
            
        }
    }
    return(designMatrix)
}
