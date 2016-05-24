#' Search for supersets and subsets
#' 
#' A function for searching the supersets and subsets of the binary drug-target interaction data.
#' 
#' @param profile_data the binary drug-target interaction matrix with row indexes as drugs and column
#' indexes as targets.
#' @return A list contains the following components:
#' \item{superset}{all the possible supersets of the input drug-target interaction data}
#' \item{subset}{all the possible subsets of the input drug-target interaction data}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' data(tyner_interaction_binary)
#' sets<-binarySet(tyner_interaction_binary[1, 1:3])
#' 
binarySet <- function(profile_data) {
    # return all possible supersets and subsets of input drug target profile data as binary set 
    # e.g. profile_data=[1 0 0 1], subset=[0 0 0 1], superset=[1 1 1 1] return multiple paramters
    
    
    # [1 1 1 1]: 2^4-1 [1 0 0 1]: 2^4-1-2^2-2^1
    length_profile_data <- length(profile_data)
    dec_zeros <- 2^(length_profile_data - which(profile_data == 0))
    dec_ones <- 2^(length_profile_data - which(profile_data == 1))
    dec_profile_data <- 2^length_profile_data - 1 - sum(dec_zeros)
    
    # gray code for zeros check dec_zeros is empty or not
    len_zeros <- length(dec_zeros)
    if (len_zeros == 0) {
        superset = dec_profile_data
        # dec_zeros_bs<-superset
    } else {
        zeros_graycode <- dec2bin(grays(len_zeros), len_zeros)
        len_zeros_graycode <- nrow(zeros_graycode)
        dec_zeros_bs <- t(matrix(rep(dec_zeros, len_zeros_graycode), ncol = len_zeros_graycode))
        superset <- rowSums(dec_zeros_bs * zeros_graycode) + dec_profile_data
    }
    
    
    
    # gray code for ones check dec_zeros is empty or not
    len_ones <- length(dec_ones)
    if (len_ones == 0) {
        # set subset NULL
        subset <- vector("numeric")
    } else {
        ones_graycode <- dec2bin(grays(len_ones), len_ones)
        # remove first line
        ones_graycode <- ones_graycode[-1, ]
        len_ones_graycode <- nrow(ones_graycode)
        if (is.null(len_ones_graycode)) 
            len_ones_graycode <- 1
        dec_ones_bs <- t(matrix(rep(dec_ones, len_ones_graycode), ncol = len_ones_graycode))
        subset <- dec_profile_data - rowSums(dec_ones_bs * ones_graycode)
    }
    
    return(list(superset = superset, subset = subset))
    #return(new("binaryset", superset = superset, subset = subset))
}

 
