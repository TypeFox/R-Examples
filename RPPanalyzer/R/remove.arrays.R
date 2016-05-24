`remove.arrays` <-
function(x,param="target",arrays2rm=c("protein","blank","housekeeping")){

        cols <- which(is.na(match(x[[3]][param,],arrays2rm)))
        
        x[[1]] <- x[[1]][,cols]
        x[[2]] <- x[[2]][,cols]
        x[[3]] <- x[[3]][,cols]
        
        return(x)
}

