ameva <-
function(tb){
        nr <- dim(tb)[1]
        nc <- dim(tb)[2]
        den <- nr*(nc-1)        
        val=chiSq(tb)/den
        return(val)  
}
