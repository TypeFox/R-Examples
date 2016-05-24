calcFlows <-
function(seq1, order.flow=c("T","A","C","G")){
    n.seq <- length(seq1)
    i.seq <- 0
    i.flo <- 0
    flow.seq <- NULL
    while(i.seq < n.seq) {
        use <- FALSE      
        flow.base <- order.flow[i.flo%%4 + 1]
        i.flo <- i.flo + 1
        flow.seq <- c(flow.seq, flow.base)
        while(as.character(seq1[i.seq +1 ]) == flow.base) {
            use <- TRUE
            i.seq <- i.seq + 1             
            if(i.seq == n.seq) break 
        }
        if (use) flow.seq[length(flow.seq)] <- tolower(flow.seq[length(flow.seq)])
    }
    flow.seq <- paste(flow.seq, collapse="")
    return(flow.seq)
}
