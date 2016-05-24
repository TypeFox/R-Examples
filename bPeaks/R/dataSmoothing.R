dataSmoothing <-
function(vecData, widthValue = 20){

#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######


    #######
    # ----- Parameter details -----------
    #
    #   widthValue : number of bases before and after the analyzed position used 
    #                to calculate the mean value
    #######

    smoothedVec = filter(vecData, rep(1/(2*widthValue+1), (2*widthValue+1)))
    # -> note that this procedure create missing values at the beginning and at the end
    # of the obtained vector
    return(as.vector(smoothedVec))
    
# end of function dataSmoothing()
}
