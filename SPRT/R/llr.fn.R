llr.fn <-
function(distribution, h0, h1) {
    
    output <- c(C.fn(distribution, h1) - C.fn(distribution, h0),
                (D.fn(distribution, h1) - D.fn(distribution, h0)) * -1)
    
    attr(output, "names") <- c("k","n")
    
    output
}
