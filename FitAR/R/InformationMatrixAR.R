`InformationMatrixAR` <-
function(phi){
toeplitz(TacvfAR(phi, length(phi)-1))
}

