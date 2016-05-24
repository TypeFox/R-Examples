
translateCoord <- function(gridCoord,mvmat,WIN){
    XT <- mvmat[,1] + gridCoord[1]
    YT <- mvmat[,2] + gridCoord[2]
    
    return(sum(inside.owin(XT,YT,WIN))/nrow(mvmat))
}
