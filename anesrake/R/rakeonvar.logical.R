rakeonvar.logical <-
function(weighton, weightto, 
    weightvec) {
    weighton <- 2 - as.numeric(weighton)
    rakeonvar.numeric(weighton, weightto, weightvec)
}

