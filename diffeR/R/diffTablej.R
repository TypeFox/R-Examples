diffTablej <- function(ctmatrix){
Quantity <- quantityDj(ctmatrix)
Exchange <- exchangeDj(ctmatrix)
Shift <- shiftDj(ctmatrix)
resT <- as.data.frame(cbind(Quantity, Exchange, Shift))
Overall <- apply(resT, 2, sum)/2
resT <- rbind(resT, Overall)
resT <- cbind(c(rownames(ctmatrix), "Overall"), resT)
colnames(resT) <- c("Category", "Quantity", "Exchange", "Shift")
return(resT)
}