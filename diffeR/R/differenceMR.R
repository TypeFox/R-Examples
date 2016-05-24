differenceMR <- function(comp, ref, eval="multiple", percent=TRUE, fact=2, population=NULL){

ctmatrix <- crosstabm(comp, ref, percent=TRUE, population=population)

if(eval=="original"){
resT <- as.data.frame(cbind(res(comp)[1], overallQtyD(ctmatrix), overallExchangeD(ctmatrix), overallShiftD(ctmatrix), overallDiff(ctmatrix)))
colnames(resT) <- c("Resolution", "Quantity", "Exchange", "Shift", "Overall")
}

if(eval=="multiple"){
res1 <- cbind(overallQtyD(ctmatrix), overallExchangeD(ctmatrix), overallShiftD(ctmatrix), overallDiff(ctmatrix))

# calculation of the factors of aggregation in powers of fact
mdim <- max(ncol(comp), nrow(comp))
maxp <- floor(log(mdim, fact))
if((log(mdim, fact) - round(log(mdim, fact), 0)) != 0) maxp <- maxp+1
factvect <-  c(1, fact^(1:maxp))

if(factvect[length(factvect)] > mdim) factvect[length(factvect)] <- mdim

resa <- data.frame(matrix(nrow=length(factvect), ncol=5))
colnames(resa) <- c("Resolution", "Quantity", "Exchange", "Shift", "Overall")
	
resa[,1] <- factvect*res(comp)[1]
resa[1,2:5] <- res1

for(i in 2:length(factvect)){
factor <- factvect[i]
			
#suppressWarnings(compAg <- aggregate(comp, fact=factor, fun=modal))
#suppressWarnings(refAg <- aggregate(ref, fact=factor, fun=modal))
#ctmatrix <- crosstabm(compAg, refAg, percent)

#compAg <- memberships(comp, fact=factor)
#refAg <- memberships(ref, fact=factor)
ctmatrix <- composite(comp, ref, factor)*100

resi <- cbind(overallQtyD(ctmatrix), overallExchangeD(ctmatrix), overallShiftD(ctmatrix), overallDiff(ctmatrix))

resa[i,2:5] <- resi
}
# add the column with the index of multiple resolutions
if(percent==FALSE) resa[,2:5] <- resa[,2:5]/100
resT <- cbind(Multiples=factvect, resa)
}
return(resT)
}