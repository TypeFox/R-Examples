pMME <-
function( mat ){
        cs = colSums( rbind(mat))
        phat = cs[1]/sum(cs)
        return(phat)
        }
