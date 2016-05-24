CMP.MSE <-
function(CMPbetas,CMPnu,x,y){

CMPresids <- constantCMPfitsandresids(CMPbetas,CMPnu,x,y)$resid

# CMPresids <- CMPfitsandresids(CMPbetas,CMPnu,x,y)$resid

MSE <- mean(CMPresids^2)

return(MSE)
}

