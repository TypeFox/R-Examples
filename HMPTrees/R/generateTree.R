generateTree <-
function(data, nreads=10000, nsamps=50, theta=0, level="genus", split="."){
if(missing(data))
stop("A valid data set is required.")
if(nreads <= 0)
stop("'nreads' must be positive and greater than 0.")
if(nsamps <= 0)
stop("'nsamps' must be positive and greater than 0.")

tempdata <- trimToTaxaLevel(data, level, FALSE, split=split)
tempdata <- transformHMPTreetoHMP(tempdata, TRUE)

if(theta > 0 && theta < 1){
dirfit <- dirmult::dirmult(tempdata) 
dirgamma <- dirfit$pi * ((1 - theta)/theta)
}else{
dirfit <- HMP::DM.MoM(tempdata)
dirgamma <- dirfit$gamma
}

gendata <- HMP::Dirichlet.multinomial(rep(nreads, nsamps), dirgamma)
colnames(gendata) <- colnames(tempdata)
gendata <- transformHMPtoHMPTree(gendata)

gendata <- buildTree(gendata, level, split)

return(gendata)
}
