flowsPanel <-
function(bed.table, genom, order.flows=c("T","A","C","G"), ret.seq=FALSE) {
# Extract amplicon sequences from genome and calculate flows
    bed.seqs <- extractSeqs(bed.table,genom)
    bed.flows <- extractFlows(bed.seqs, order.flows)
    out.table <- cbind(bed.table, bed.flows)
# Print statistics
    namplicon <- as.numeric(as.matrix(out.table[, "TARGET_LENGTH"]))
    nflows <- as.numeric(as.matrix(out.table[, c("N_FLOW_FORWARD", "N_FLOW_REVERSE")]))
    cat("Amplicon length (bp), average (min-max): " , round(mean(namplicon)), " (", min(namplicon), "-", max(namplicon), ")\n", sep="")
    cat("Required number of flows, average (min-max): " , round(mean(nflows)), " (", min(nflows), "-", max(nflows), ")\n", sep="")
# Return Results
if(ret.seq) return(out.table)
}
