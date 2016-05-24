add.dimension <- function(x,t.o){

    result <- array(dim=c(t.o$count.model,
                               t.o$eval.steps,
                               NCOL(x)))
    if(is.null(dim(x))){
        dim(x) <- c(NROW(x),NCOL(x))
    }

    for(counter in 1:t.o$count.model){
        start <- (counter-1)*t.o$eval.steps+1
        end <- (counter)*t.o$eval.steps
        result[counter,,] <- as.matrix(x[start:end,])
        #achtung: Am Anfang von jedem Lauf NA setzen!
        result[counter,1:floor(t.o$window.size/t.o$step.size),] <- NA
    }
   return(result)

}
