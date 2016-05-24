RelCoefLine <-
function(LineGenom=matrix(0,nrow=0,ncol=0) , Freq=matrix(0,nrow=0,ncol=0) , LinePop=rep(0,0) , Combination=NULL , NbInit=5 , Prec=10^(-4) , NbCores=NULL){

#Check if a genotype matrix is provided
if (nrow(LineGenom)==0){
stop("Genotype matrix should be provided")
}

#Check if a allelic frequencies matrix is provided
if (nrow(Freq)==0){
stop("Freq matrix should be provided")
}


#Only one among ParentalLineGenom and IndividualGenom should be provided 
if (nrow(LineGenom) != (nrow(Freq))){
stop("Nrow differs between the Genotype matrix and Pop frequencies.")
}

if (length(LinePop)==0){
LinePop <- rep(1,ncol(LineGenom))
}

#Check if NbParent in LineGenom = NbParent in LinePop
if (length(LinePop) != ncol(LineGenom)){
stop("Nb of parents differ between LinePop and LineGenom")
}

#Check if NbPop and NbFreq correspond
if (max(LinePop) > dim(Freq)[2]){
stop("Number of Pop frequencies < number of Pop")
}

NbLine <- ncol(LineGenom)

#Number of cores selected
if (Sys.info()[['sysname']]=="Windows"){
NbCores <- 1
if (length(NbCores)!=0){
print("NbCores > 1 is not supported on Windows, NbCores is set to 1")
}
} else {
if (length(NbCores)==0){
NbCores <- detectCores()-1
}
}

if (length(Combination)==0){
#Every couples of hybrids are studied
comb <- combn(1:NbLine , 2 , simplify = F)



CoupleTwoLines <- mclapply(comb , function(x) .RelatednessLineCouple(LineGenom[,x],Freq,LinePop[x],NbInit,Prec) , mc.cores=NbCores)
#Formatting the result
MatDelta <- sapply(CoupleTwoLines, function(x) x$Delta)
mat <- matrix(0,NbLine,NbLine)
Delta <- .MatTriSup(mat,MatDelta[2,])
diag(Delta) <- rep(1,NbLine)

} else {
NamesCombination <- sapply(1:length(Combination) , function(x) paste0(Combination[[x]][1],"/",Combination[[x]][2]))
#Only selected couples of hybrids are studied
CoupleTwoLines <- mclapply(Combination , function(x) .RelatednessLineCouple(LineGenom[,x],Freq,LinePop[x],NbInit,Prec) , mc.cores=NbCores)
Delta <- lapply(CoupleTwoLines , function(x) x$Delta[2])
names(Delta) <- NamesCombination
}

    return(Delta)
}
