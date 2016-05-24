baseLineCalc <-
function(covData){

    #######
    # ----- Parameter details -----------
    #
    #   covData    : vector with dataset (all genome)
    #######
    
    return(sum(as.numeric(covData))/length(covData))

# end of function baseLineCalc
}
