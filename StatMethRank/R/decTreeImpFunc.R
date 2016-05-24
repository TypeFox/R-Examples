#' Decision Tree for Ranking Data-Impurity Function Approach
#'
#' See Section 10.1 of the book for more details
#'
#' @param ITEM_SIZE number of items
#' @param RANK_SIZE top-q size
#' @param CART_SIZE number of cross validation in tree pruning
#' @param DATA_FILE please specify the input file as fullpath/input_file.txt
#' @param INFO_FILE please specify the input file as fullpath/input_file.txt, note that the data should be 
#' seperated by tab.
#' @param OUTPUT_FILE name of the output file
#' @param mode Spliting Criterion: E = Entropy G = Gini S = Statistical Test C = Chi-square test
#' @param TMODE When mode = "C" P = Pearson Chi-square test; L =  Likelihood ratio test
#' @param usePairwise use pairwise comparison model or not
#' @param isPW_DATA use pairwise data or top k-ranked data
#' @param prediction M : mean, F : frequency, C : center
#' @param NODE_SIZE usually one tenth of the number of observations 
#' @param MIN_NODE_SIZE min NODE_DIZE
#' @param RANK_ITEM top-q measure (1-3)
#' @param ST_ALP level of significance of Statistical test
#' @param CHI_ALP level of significance of Chi-square test
#' @param TRAIN_PROP proportion of training data; effective only when CV_TEST = false
#' @param useCV_TEST use 10-fold CV testing or not 
#' @param TEST_STAGE the stage of the 10-fold CV testing; effective only when CV_TEST = true; value starts from 0 to V-1
#' 
#' @return a list contains the information of the tree
#' @examples 
#' #see example 10.1.4.R on the websited for more details.
#' @export
#' @author Li Qinglong <liqinglong0830@@163.com>, William Lai 
#' @references Decision tree modeling for ranking data, Yu, P.L.H.,Wan, W. M.,& Lee, P.H.(2010) 
decTreeImpFunc <- function( ITEM_SIZE, RANK_SIZE, CART_SIZE, 
                            DATA_FILE, INFO_FILE, OUTPUT_FILE,
                            mode = c("E", "G", "S", "C"), TMODE = c("L", "P"), 
                            usePairwise = FALSE, isPW_DATA = FALSE,
                            prediction = c("M", "F", "C"), NODE_SIZE = 400, 
                            MIN_NODE_SIZE = 400, RANK_ITEM = 2, ST_ALP = 0.95, 
                            CHI_ALP = 0.2, TRAIN_PROP = 1.0, useCV_TEST = FALSE, 
                            TEST_STAGE = 0)
{
    if(!(is.logical(usePairwise) && is.logical(isPW_DATA) && is.logical(useCV_TEST)))
    {
        stop("usePairwise, isPW_DATA and useCV_TEST should be in logical type!")    
    }
    mode = match.arg(mode)
    TMODE = match.arg(TMODE)
    prediction = match.arg(prediction)
    str_var_info = readLines(INFO_FILE)
    mytree = decTreeImpFuncCpp(ITEM_SIZE, RANK_SIZE, CART_SIZE, 
                        DATA_FILE, INFO_FILE, OUTPUT_FILE,
                        mode, TMODE, usePairwise, isPW_DATA,
                        prediction, NODE_SIZE, MIN_NODE_SIZE,
                        RANK_ITEM, ST_ALP, CHI_ALP, TRAIN_PROP, 
                        useCV_TEST, TEST_STAGE)
    
    return(c(mytree, Var_Info = list(str_var_info)))
}