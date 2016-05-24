extractFlows <-
function(seqs, order.flows=c("T","A","C","G")){
    n <- length(seqs)
    cat("Calculating flow sequences (. = 10 amplicons):\n")
    flows <- matrix(nrow=n ,ncol=6)
    colnames(flows) <- c("TARGET_SEQUENCE","TARGET_LENGTH","N_FLOW_FORWARD","N_FLOW_REVERSE","FLOW_SEQ_FORWARD","FLOW_SEQ_REVERSE")
    for(i in 1:n) {
        s <- DNAString(seqs[i])
        flows[i, 1] <- as.character(s)
        flows[i, 2] <- length(s)
        flows[i, 5] <- calcFlows(reverse(complement(s)), order.flows)
        flows[i, 6] <- calcFlows(s ,order.flows)
        flows[i, 3] <- nchar(flows[i, 5])
        flows[i, 4] <- nchar(flows[i, 6])
        if(i%%10 == 0) cat(".")
        if(i%%500 == 0 || i == n) cat(i, "/", n ,"\n", sep="")
    }
    return(flows)
}
