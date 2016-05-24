evalContour <- function(AAMat, BBVec = NULL, IPVec = NULL){
#evalContour <- function(AAMat, BBVec, IPVec), output: CST
#given the set of inequalities AAMat%*%ZVec <= BBVec, this function identifies nonredundant
# constraints and computes the vertices of the resulting contour (with an interior point IPVec)
#  defined by the inequalities, more precisely the vertices of the intersection of the contour
#   with the artificial bounding hypercube [-Bound, Bound]^dim(AAMat)[2]
#
#CST     ... the (contour) list containing the output and some summary statistics
#  CST$Status == 0 ... OK
#  CST$Status == 2 ... the contour seems virtually empty
#  CST$Status == 3 ... a lp numerical failure (currently not possible to achieve)
#  CST$Status == 4 ... the number of input parameters is too low
#  CST$Status == 5 ... AAMat is not a numerical matrix with two to six columns
#  CST$Status == 6 ... BBVec is not a numerical column vector of the same length as the first column of AAMat
#  CST$Status == 7 ... IPVec is not a numerical column vector of the same length as the first row of AAMat
#  CST$TVVMat      ... the vertices of the resulting contour (in rows)
#  CST$TKKMat      ... the facets of the resulting contour (in rows); see help(convhulln) for the meaning of TKKMat
#  CST$Vol         ... the volume of the contour region
#  CST$NumF        ... the number of contour facets
#  CST$NumV        ... the number of contour vertices
#
#the contour can then be plotted, see the examples
CST<-list()
#output initialization
CST$Status <-  0
CST$TVVMat <- NULL
CST$TKKMat <- NULL
CST$Vol    <- NaN
CST$NumF   <- NaN
CST$NumV   <- NaN

#checking the input

if (is.null(BBVec)){CST$Status <- 4; return(CST)}
OutList <- checkArray(AAMat,0,0,NULL,c(1, Inf),c(2, 6),1); AAMat <- OutList[[1]]; Status <- OutList[[2]]
if (Status > 0){CST$Status <- 5; return(CST)}
NColAA <- dim(AAMat)[2]
if (is.null(IPVec)){IPVec <- matrix(0,NColAA,1)}
OutList <- checkArray(BBVec,0,0,NULL,c(dim(AAMat)[1],dim(AAMat)[1]),c(1, 1),1); BBVec <- OutList[[1]]; Status <- OutList[[2]]
if (Status > 0){CST$Status <- 6; return(CST)}
OutList <- checkArray(IPVec,0,0,NULL,c(NColAA, NColAA),c(1, 1),1); IPVec <- OutList[[1]]; Status <- OutList[[2]]
if (Status > 0){CST$Status <- 7; return(CST)}














#OptionsCA ... the options for convhulln used for vertex enumeration
if (NColAA <= 3){OptionsCA <- c('C-0', 'Qt', 'Qs', 'Pp', 'FA')} else {OptionsCA <- c('Qt', 'Qx', 'Qs', 'QbB', 'Pp', 'FA')}

#TVDig ... the almost identical vertices are eliminated by round(VV*10^TVDig)/10^TVDig
TVDig <- 7

#'normalizing'
EpsN <- 1e-10
NormAAVec <- sqrt((AAMat^2)%*%matrix(1,NColAA,1))
if (!all(NormAAVec==1)){
    NormAAVec <- NormAAVec + (NormAAVec <= EpsN)
    AAMat <- AAMat/(NormAAVec%*%matrix(1,1,NColAA))
    BBVec <- BBVec/NormAAVec
}

#adding bounding constraints
#Bound ... one half of the edge length of the artificial bounding box for the contour
Bound <- 1e3
AAMat <- rbind(AAMat, diag(1,NColAA), -diag(1,NColAA))
BBVec <- rbind(BBVec, Bound*matrix(1,NColAA,1), Bound*matrix(1,NColAA,1))
NRowAA <- dim(AAMat)[1]

#finding an interior point IP of the region AAMat%*%Z <= BBVec if IPVec is not
# inside the contour, by solving an auxiliary linear programming problem
#EpsIP ... a small positive number used for identification of the interior points
EpsIP <- 1e-7
if (!all((AAMat%*%IPVec) -BBVec <= -EpsIP)){ #if IPVec is not well inside the contour
    #forming the linear programming matrices and vectors
    ExtAAMat <- cbind(AAMat, -AAMat, matrix(1,dim(AAMat)[1],1), diag(dim(AAMat)[1]))
    VCST <- list()

    VCST$l <- 2*NColAA + 1 + NRowAA #the number of non-negative variables
    CCVec <- matrix(0,VCST$l,1); CCVec[2*NColAA+1] <- -1

    #solving the auxiliary linear programming problem
    lp0 <- lp("min", CCVec, ExtAAMat, array('==', dim=c(NRowAA,1)), BBVec)
    ZZVecP <- lp0$solution[1:NColAA]; ZZVecM <- lp0$solution[(NColAA+1):(2*NColAA)];
    IPVec <- matrix(ZZVecP-ZZVecM,NColAA,1)

    #checking the optimization results
    if ((lp0$status > 0)||(!all((AAMat%*%IPVec) -BBVec <= -EpsIP))){
      CST$Status <- 2; return(CST)

    }
  rm(ExtAAMat, CCVec, VCST)
}

#obtaining vertices VV of the contour (vertex enumeration)

NewBBVec <- BBVec - AAMat%*%IPVec
DDMat <- AAMat / (NewBBVec%*%matrix(1,1,NColAA))
OutList <- convhulln(DDMat, OptionsCA); KKMat <- OutList$hull
NRowKK <- dim(KKMat)[1]
GGMat <- matrix(0,NRowKK,NColAA)
for (i in 1:NRowKK){GGMat[i,] <- solve(DDMat[KKMat[i,],], matrix(1,NColAA,1))}
rm(DDMat, KKMat, NewBBVec)
TVVMat <- GGMat + matrix(1,NRowKK,1)%*%t(IPVec)
rm(GGMat)

#rounding the vertices and preparing the output
CST$TVVMat <- unique(round(TVVMat*10^TVDig)/10^TVDig)
CST$TVVMat <- CST$TVVMat[do.call(order, split(CST$TVVMat, col(CST$TVVMat))),]
OutList <- convhulln(CST$TVVMat, OptionsCA); CST$TKKMat <- OutList$hull; CST$Vol <- OutList$vol; CST$Area <- OutList$area
CST$NumV <- dim(CST$TVVMat)[1]
CST$NumF <- dim(CST$TKKMat)[1]
return(CST)
}
