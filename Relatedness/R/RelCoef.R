RelCoef <-
function(IndividualGenom=matrix(0,nrow=0,ncol=0),ParentalLineGenom=matrix(0,nrow=0,ncol=0) ,Freq=matrix(0,nrow=0,ncol=0),Crossing=matrix(0,nrow=0,ncol=0),ParentPop=rep(0,0),Combination=list(),Phased=FALSE,Details=FALSE,NbInit=5,Prec=10^(-4),NbCores=NULL){
  
#### First check the Inputs

#One Genotype matrix is provided
if ((nrow(IndividualGenom)==0) & (nrow(ParentalLineGenom)==0)){
stop("No genotype matrix provided")
}

#Only one among ParentalLineGenom and IndividualGenom should be provided 
if ((nrow(IndividualGenom)!=0) & (nrow(ParentalLineGenom)!=0)){
stop("Either a Parental genotype matrix OR a hybrid genotype matrix should be provided")
}

#Check if a allelic frequencies matrix is provided
if (nrow(Freq)==0){
stop("Freq matrix should be provided")
}

#Check which Genotype Matrix is provided
if (nrow(IndividualGenom)!=0){
Genom <- IndividualGenom
#Check the other input
if ((nrow(ParentalLineGenom)!=0) | (nrow(Crossing)!=0) | (length(ParentPop!=0))){
stop("If IndividualGenom is provided, ParentalLineGenom, Crossing and ParentPop cannot be provided")
}

Crossing <- matrix(seq(1:ncol(Genom)),ncol=2,byrow=T)

ParentPop <- rep(1,ncol(Genom))

} else {
Genom <- ParentalLineGenom
#Check if crossing is provided
if (nrow(Crossing)==0){
stop("Crossing matrix should be provided")
}

#Check if NbParent in Crossing < NbParent in ParentalLineGenom
if (max(Crossing) > ncol(Genom)){
stop("Nb of parents in crossing > Nb of parents in ParentalLineGenom")
}

#Check if ParentPop is provided
if (length(ParentPop)==0){
ParentPop <- rep(1,ncol(Genom))
}

#Check if NbParent in Crossing = NbParent in ParentPop
if (length(ParentPop) < max(Crossing)){
stop("Nb of parents differ between ParentPop and Crossing")
}

#Check if NbParent in ParentalLineGenom = NbParent in ParentPop
if (length(ParentPop) < ncol(Genom)){
stop("Nb of parents differ between ParentPop and ParentalLineGenom")
}

#Check if NbPop and NbFreq correspond
if (max(ParentPop) > dim(Freq)[2]){
stop("Number of Pop frequencies < number of Pop")
}

#Check if Phased is FALSE
if (Phased==FALSE){
print("Argument Phased is change in TRUE")
Phased=TRUE
}
    }

#Check if NbMarkers in Freq = NbMarkers in Genotype matrix
if (dim(Freq)[1] != dim(Genom)[1]){
stop("Nrow differs between the Genotype matrix and Pop frequencies.")
}

NbIndividual <- dim(Crossing)[1]

if(Phased==F){
    
        NbIBD <- 9
        
} else {
    
        NbIBD <- 15
        
}

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

Crossing <- t(Crossing)

if (length(Combination)==0){
#Every couples of hybrids are studied
comb <- combn(1:NbIndividual , 2 , simplify = F)



CoupleTwoHybrids <- mclapply(comb , function(x) .RelatednessCouple(Genom[,c(Crossing[,x])],Freq,Crossing[,x],ParentPop[c(Crossing[,x])],Phased,NbInit,Prec) , mc.cores=NbCores)
    CoupleOneHybridRepeted <- mclapply(1:NbIndividual , function(x) .RelatednessCouple(Genom[,c(Crossing[,c(x,x)])],Freq,Crossing[,c(x,x)],ParentPop[c(Crossing[,c(x,x)])],Phased,NbInit,Prec) , mc.cores=NbCores)
#Formatting the result
MatDelta <- sapply(CoupleTwoHybrids, function(x) x$Delta)
mat <- matrix(0,NbIndividual,NbIndividual)
  DeltaTri <- lapply(1:NbIBD , function(x) .MatTriSup(mat,MatDelta[x,]))
DeltaDiag <- sapply(CoupleOneHybridRepeted, function(x) x$Delta)
Delta <- lapply(1:NbIBD , function(x) .AddDiag(DeltaTri[[x]],DeltaDiag[x,]))
names(Delta) <- paste0("Delta",1:NbIBD)

} else {
NamesCombination <- sapply(1:length(Combination) , function(x) paste0(Combination[[x]][1],"/",Combination[[x]][2]))
#Only selected couples of hybrids are studied
CoupleTwoHybrids <- mclapply(Combination , function(x) .RelatednessCouple(Genom[,c(Crossing[,x])],Freq,Crossing[,x],ParentPop[c(Crossing[,x])],Phased,NbInit,Prec) , mc.cores=NbCores)
Delta <- lapply(CoupleTwoHybrids , function(x) x$Delta)
names(Delta) <- NamesCombination
}
  
#Print pictures of relatedness
if (Details==TRUE){
if (Phased==TRUE){
.TheFifteenDeltaGraph()
}else{
.TheNineDeltaGraph()
}
}
    return(Delta)
}
