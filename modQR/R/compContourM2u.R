compContourM2u   <- function(Tau = 0.2, YMat = NULL, XMat = NULL, CTechST = NULL){
#compContourM2u   <- function(Tau = 0.2, YMat = NULL, XMat = NULL, CTechST = NULL), output: COutST
#computing (regression) quantile regions by the algorithm described in the CS paper:
#  Paindaveine, D. and \v{S}iman, M. (2012):
#  Computing Multiple-Output Regression Quantile Regions from Projection Quantiles.
#  Computational Statistics 27, 29--49.
#
#Tau      ... the quantile level in [0, 0.5]
#YMat     ... the response matrix with two to six columns
#XMat     ... the design matrix including the (first) intercept column
#CTechST  ... the list with some parameters regarding the computation, possibly modified by the user; see getCTechSTM2u
#
#COutST   ... the list containing some useful information about the computation, with the following fields:
#  CharST        ... the list with some default or user-defined output, provided by CTechST$getCharST (and initialized by getCharSTM2u)
#  (self-explaining warnings and error messages regarding the input and computation)
#  CTechSTMsgS   ... '' or a warning about the (overcome) problems with input CTechST
#  ProbSizeMsgS  ... '' or a warning about the large size of the problem
#  CompErrMsgS   ... '' or the description of the terminal error interrupting the computation
#  TauMsgS       ... '' or a warning about the perturbation of Tau
#  (scalar characteristics of the computation)
#  NDQFiles      ... the counter of (possible) output files, i.e. as if CTechST$OutSaveI = 1
#  NumB          ... the counter of optimal bases considered
#  MaxLWidth     ... 0 or the maximum width of one layer (if !CTechST$D2SpecI)
#  NIniNone      ... the number of trials when the initial solution(s) could not be found at all
#  NIniBad       ... the number of trials when the found initial solution(s) did not have the right number of nonzero coordinates
#  NSkipCone     ... the number of skipped cones (where the interior point could not be found)
#  (the vector indicating the position of individual observations with respect to the contour)
#  PosVec        ... PosVec[i] <- 0/1/2 if the i-th observation is in/on/out of the contour (after the successful computation)
#see the accompanying paper for further details
#
#Prerequisities:
#- R with packages lpSolve and geometry + package rgl for 3D demo examples


#Decoding variable names:
#suffixes:
# Vec (column vector), Row (row vector), Mat (matrix), I (indicator), S (string),
# F (format), CA (cell array), ST (list), SA (list array); Octave structures are replaced with lists in R
#other possibilities: T ... tilde, H ... hat, I ... inverse, D ... denominator

#Initializing the output I
#=========================

COutST <- list()
#error or warning messages (errors interrupt/terminate the computation but warnings do not)
COutST$CTechSTMsgS  <- ''
COutST$ProbSizeMsgS <- ''
COutST$CompErrMsgS  <- ''
COutST$TauMsgS      <- ''

#Checking the input
#==================

if (is.null(CTechST)){CTechST <- getCTechSTM2u(); Status <- 0
}else{OutList <- checkCTechSTu(CTechST,2); CTechST <- OutList[[1]]; Status = OutList[[2]]}

if (Status == 1){COutST$CTechSTMsgS <- 'CTechST had to be replaced with the default list!'
} else if (Status > 1){COutST$CTechSTMsgS <- 'Some fields of CTechST had to be replaced with the default ones!'}


if (!is.matrix(YMat)){COutST$CompErrMsgS <- 'The number of input parameters is too low!'; return(COutST)}

OutList <- checkArray(Tau,0,0,c(0, 0.5),c(1, 1),c(1, 1),1)
if (OutList[[2]] > 0){COutST$CompErrMsgS <- 'Tau is not a real number in the interval [0, 0.5]!'; return(COutST)}
OutList <- checkArray(YMat,0,0,NULL,c(1, Inf),c(2, 6),1)
if (OutList[[2]] > 0){COutST$CompErrMsgS <- 'YMat is not a numeric matrix with at least one row and two to six columns!'; return(COutST)}
if (!is.matrix(XMat)){XMat <- matrix(1, dim(YMat)[1], 1)
}else{
    OutList <- checkArray(XMat,0,0,NULL,c(1, Inf),c(1, Inf),1)
    if (OutList[[2]] > 0){COutST$CompErrMsgS <- 'XMat is not a numeric matrix with at least one row and one column!'; return(COutST)}
}

N <- dim(YMat)[1]
M <- dim(YMat)[2]
P <- dim(XMat)[2]
if (dim(XMat)[1] != N){COutST$CompErrMsgS <- 'The matrices XMat and YMat differ in the number of rows!'; return(COutST)}
if (N <= M+P-1){COutST$CompErrMsgS <- 'N should be strictly higher than M+P-1!'; return(COutST)}
igProblemI <- FALSE
if (M == 2){BigProblemI <- (N>10000)||(P>10)}
if (M == 3){BigProblemI <- (N>500)||(P>5)}
if (M == 4){BigProblemI <- (N>150)||(P>3)}
if (M >  4){BigProblemI <- TRUE}
if (BigProblemI){COutST$ProbSizeMsgS <- 'The problem is too large (and unsupported), the computation may fail easily!'}









#Initializing the output II
#==========================

COutST$CharST      <- CTechST$getCharST(Tau, N, M, P, NULL, NULL, 1) #initializing the user-defined output
COutST$PosVec      <- matrix(0,N,1) #the vector indicating the position of individual observations with respect to the contour:
                                    # COutST$PosVec[i] = 0/1/2 if the i-th observation is in (or unclassified)/on or out of the contour
COutST$NDQFiles    <- 0             #the counter of (possible) output files, i.e. as if CTechST$OutSaveI = 1
COutST$NumB        <- 0             #the counter of optimal bases considered
COutST$MaxLWidth   <- 0             #the maximum width of one layer (if !CTechST$D2SpecI)
COutST$NIniNone    <- 0             #the number of trials when the initial solution(s) could not be found at all
COutST$NIniBad     <- 0             #the number of trials when the found initial solution(s) did not have the right number of nonzero coordinates
COutST$NSkipCone   <- 0             #the number of skipped cones (where the interior point could not be found)

#correcting the initial setting
if (M > 2){CTechST$D2SpecI   <- 0}
if (M > 3){CTechST$ArchAllFI <- 1}

#Initializing technical parameters
#=================================

Dist                   <- 0.9   #the distance of auxiliary interior points of the cones and of the facet 'centers' from the origin
PFNormI                <- 1     #1/0 ... if all the rows of PFModMat (PF condition: PFModMat[1:(N+M),]%*%u <= 0) should be normalized (1) or not (0)

#accuracy parameters
EpsTau                 <- 2e-6  #the perturbation of Tau to avoid (Tau*N) close to an integer; EpsTau should be less than 1/(2*N)
EpsPF                  <- 5e-11 #the bound for distinguishing zeros in the PF condition, in: if !all(Z1Vec > -EpsPF) ...
if (CTechST$D2SpecI){
    EpsTheta           <- 1e-15 #the minimum change of Theta (= the cone angle) to be considered significant
}else{
    if (M > 2){
        EpsFV          <- 5e-12 #the bound for identifying facet vertices (used for computing facet 'centers'), in: ActFCMat = (abs(ActFNMat%*%t(VVMat)) <= EpsFV)%*%VVMat
        EpsInt         <- 1e-11 #the bound distinguishing the points inside the cones, in: IsIn = all(PFModMat[1:(N+M),1:M]%*%IPVec < -EpsInt)
        FIPrec         <- 1e9   #the accuracy/precision of cone facet identifiers (round(x*FIPrec)/FIPrec is used)
        IniStep        <- 1e-5  #the initial step inwards for the adaptive selection of a point inside the cone
        MinStep        <- 1e-11 #the minimum step inwards for the adaptive selection of a point inside the cone (if still large then the adjacent cone is skipped)
    }else{ #if (M == 2)
        EpsFV          <- 3e-13
        EpsInt         <- 1e-11
        FIPrec         <- 3e10
        IniStep        <- 1e-6
        MinStep        <- 1e-11
    }
}

#capacity constraints
NNInDQ                 <- 4e6   #the desired total number of numbers in DQMat (which is the matrix for keeping the results before storing them on the disk,
                                # processing them, or forgetting them). It should be at least 1000 or so.
MaxNB                  <- 5e7   #the maximum number of optimal bases allowed (when reached then the (faulty or too time-consuming) computation is interrupted)

if (!CTechST$D2SpecI){
    if (!CTechST$ArchAllFI){
        NArchLayer       <- 2*(M+P)+4   #the number of the last layers whose cone facet identifiers are stored. It is wise to set it not too low.
        #MaxNFIPerLayer ... the maximum number of facet identifiers resulting from 1 layer (when reached then the (faulty or too memory-consuming) computation is interrupted)
        if (M == 2){MaxNFIPerLayer <- 10}else{MaxNFIPerLayer <- 2.5e4}
    }else{ #if CTechST$ArchAllFI
        MaxNF    <- 1e6 #the maximum number of used facet identifiers stored (when reached then the (faulty or too memory-consuming) computation is interrupted)
    }
}

#options for qhull/convhulln
if (M <= 3){QHullOptionsCA <- c('C-0', 'Qt', 'QbB', 'Pp')}else{QHullOptionsCA <- c('Qt', 'Qx', 'Qs', 'QbB', 'Pp')}




DataCharF <- paste(CTechST$OutFilePrefS, '%06d.dqo', sep='')   #the format of file output name(s)

#saving and changing warning states
WarnST <- options(warn=-1)$warn

#====================================
#The core parametric programming part
#====================================
#It is highly recommended to read it side by side with the CS paper whose algorithm, notation, and text are followed

#perturbing Tau in the location case (to prevent degeneracy/multiple solutions)
if (P == 1){
    BadTau <- round(Tau*N)/N
    GapTau <- Tau - BadTau
    if (abs(GapTau) <= EpsTau){
        if (GapTau > -1e-14){Tau <- BadTau + EpsTau}else{Tau <- BadTau - EpsTau}
        COutST$TauMsgS <- sprintf('Tau has been perturbed to Tau = % 11.10f!\n', Tau)
        if (CTechST$ReportI){print(COutST$TauMsgS)}



    }
}

#auxiliary constants
PPlusM                <- P+M
P2PlusM2              <- 2*PPlusM
Ones1PRow             <- matrix(1,1,P)
if (!CTechST$D2SpecI || (PFNormI == 1)){
    Ones1MRow         <- matrix(1,1,M)
    OnesM1Vec         <- matrix(1,M,1)
}
#BBVec ... the vector for describing the cones intersected by hypercubes, in: PFModMat%*%UVec <= BBVec
if (!CTechST$D2SpecI){ BBVec          <- rbind(matrix(0,N+M,1), matrix(1,M,1))
}else{                 OneToNPlusMVec <- t(1:(N+M)) #for identifying the tightest constraints
}

#if the intercept is equal to the unit vector, then some formulae can be simplified
NoWeightI <- all(XMat[,1] == matrix(1,N,1))

#preparing the matrix DQMat of size MaxNRowQ x NColQ for storing the results temporarily, before processing them at once with CTechST$getCharST
# and/or storing them on the disk. Then the matrix is refreshed and the process of recording the results continues.
if (CTechST$BriefOutputI){NColQ <- 1 + 1 + M + P*M + M         #c(ConeID, Nu, UVec, vec(ACOMat), MuBRow)
}else{NColQ <- 1 + 1 + M + P*M + M + P + P                     #c(ConeID, Nu, UVec, vec(ACOMat), MuBRow, MuR0Row, IZ)
}
MaxNRowQ  <- ceiling(NNInDQ/NColQ)
DQMat     <- matrix(0, MaxNRowQ, NColQ)

#preparing the initial (corner) directional vectors in the rows of U0Mat (for computing the first solution from scratch)
if (CTechST$D2SpecI || !CTechST$CubRegWiseI){
    U0Mat <- -matrix(1, 1, M)/sqrt(M)
}else{
    U0Mat <- matrix(0, 2^M, M)
    for (i in 0:(2^M-1)){
        for (j in 1:M){
             if (bitwAnd(i, bitwShiftL(1,j-1))){U0Mat[i+1,M-j+1] <- 1}
        }
    }
    U0Mat <- (U0Mat - (U0Mat == 0))/sqrt(M)
}

#auxiliary counters
NARDQ        <- 0   #counter of active (nonzero) rows in DQMat
ConeID       <- 0   #identifier of the cone whose details are stored

for (IndU0Vec in 1:dim(U0Mat)[1]){

    #getting the first directional quantile and optimal basis
    #--------------------------------------------------------

    U0Vec <- t(U0Mat[IndU0Vec,])
    OutList <- findOptimalBasisM2FromScratch(M, N, P, Tau, U0Vec, YMat, XMat); TST <- OutList[[1]]; U0Vec <- matrix(OutList[[2]],M,1); IsFound <- OutList[[3]]; TrialST <- OutList[[4]]
    if (!IsFound){
        COutST$CompErrMsgS <- 'No reliable initial optimal solution with the right number of nonzero coordinates has been found!'
        return(COutST)
    }
    COutST$NIniNone <- COutST$NIniNone + TrialST$NNotFound
    COutST$NIniBad  <- COutST$NIniBad  + TrialST$NBad
    if (CTechST$ReportI){
        print(sprintf('   U0Vec = [%s], (NNotFound = % 2d, NBad = % 2d)', paste(U0Vec, collapse=','), TrialST$NNotFound, TrialST$NBad))

    }
    if (CTechST$D2SpecI){U0Angle <- -pi+atan(U0Vec[2]/U0Vec[1])}

    #initializing storage variables (for monitoring visited cones and the cones to be explored)
    #------------------------------------------------------------------------------------------

    if (CTechST$D2SpecI){
        ThetaNew <- U0Angle         #starting angle
        FCVec <- Dist*U0Vec         #standardized facet center/interior point (but not now)
    }else{
        #the data lists for storing some auxiliary results
        PresDatSA  <- list()  #... from the topical layer
        FutDatSA   <- list()  #... from the future layer
        NARFutDat  <- 0       #... the counter of active items in NARFutDat

        PresDatSA[[1]]        <- list() #initializing PresDatSA
        PresDatSA[[1]]$FNVec  <- matrix(0, M, 1) #the facet normal vector (but not now)
        PresDatSA[[1]]$FCVec  <- Dist*U0Vec      #the standardized facet center (interior point) u_f (but not now)
        PresDatSA[[1]]$TST    <- TST             #the corresponding TST list containing I_B, I_C, I_R and other useful characteristics

        #storing used facet identifiers
        if (CTechST$ArchAllFI){
            #a single array for storing identifiers of all visited facets
            FIArchiveMat  <- matrix(0, MaxNF, M)
            NARFIArchive <- 0      #the counter of active rows/items in FIArchiveMat
        }else{
            #a cell array for storing facet identifiers from the last NArchLayer layers
            #each cell corresponds to one layer
            FIArchMatCA  <- list()
            for (IndL in 1:NArchLayer){FIArchMatCA[[IndL]] <- matrix(0, MaxNFIPerLayer, M)}
            NARFIArchVec <- matrix(0, NArchLayer, 1)   #the vector of counters of active rows/items in FIArchMatCA
        }
    }

    #auxiliary counters
    NARPresDat   <- 1   #the number of facets investigated in the current step
    NumL         <- 0   #the number of layers (i.e. 'steps' of the breadth-first search algorithm)

    JMat <- diag(as.vector(sign(U0Vec)))

    while (COutST$NumB <= MaxNB){ #while the number of bases investigated is not too high

        NumL <- NumL + 1
        COutST$MaxLWidth <- max(COutST$MaxLWidth, NARPresDat)

        if (CTechST$ReportI && (!CTechST$D2SpecI)){
            print(sprintf('         Layer No.:  % 8d  | % 6d investigated facet(s)\n', NumL, NARPresDat))

        }

        for (FacetIndex in 1:NARPresDat){ #for each facet from the topical layer

            COutST$NumB <- COutST$NumB + 1

            #decoding the input information
            #------------------------------

            if (!CTechST$D2SpecI){
                FNVec  <- PresDatSA[[FacetIndex]]$FNVec    #topical FNVec
                FCVec  <- PresDatSA[[FacetIndex]]$FCVec    #topical FCVec
                TST    <- PresDatSA[[FacetIndex]]$TST      #topical TST (containing various index vectors)
                if (!CTechST$CubRegWiseI){JMat <- diag(((TST$Ib %% 2) == 1) - ((TST$Ib %% 2) == 0))}
            }
            Pi   <- TST$Pi                         #the number of positive residuals
            Nu   <- TST$Nu                         #the number of negative residuals

            BasNegAVec <- ((TST$Ia %% 2) == 0)     #indicators of basic a_{i-}'s
            BasPosAVec <- !BasNegAVec              #indicators of basic a_{i+}'s


            #finding auxiliary matrices (PFModMat) and vectors
            #-------------------------------------------------

            EJMat <- YMat[TST$IR,]
            if ((P > 1)||(!NoWeightI)){
                FMat <- as.matrix(XMat[TST$IR,])
                FMat[,BasPosAVec] <- -FMat[,BasPosAVec]
                LMat <- solve(FMat[1:P,1:P])
            }else{ #(P == 1 && NoWeightI)
                LMat <- sign(BasNegAVec[1] - 0.5)
                FMat <- LMat * matrix(1,N,1); LMat <- matrix(LMat,1,1)
            }
            KMat <- -LMat%*%EJMat[1:P,]

            #PFModMat: the (PF) condition can be expressed as PFModMat[1:(N+M),]%*%u <= 0
            if (Pi == 0){
                PFModMat <- rbind(-JMat, -KMat, (EJMat[(P+1):N,] + FMat[(P+1):N,]%*%KMat), JMat)
                if ((P > 1)|| !NoWeightI){MuR0Row <- -(1-Tau)*(matrix(1,1,N-P)%*%FMat[(P+1):N,])%*%LMat}
                MuBRow <- -(1-Tau)*(matrix(1,1,N-P)%*%PFModMat[(PPlusM+1):(M+N),])
            }else if (Nu == 0){
                PFModMat <- rbind(-JMat, -KMat, -(EJMat[(P+1):N,] + FMat[(P+1):N,]%*%KMat), JMat)
                if ((P > 1)|| !NoWeightI){MuR0Row <- Tau*(matrix(1,1,N-P)%*%FMat[(P+1):N,])%*%LMat}
                MuBRow <- -Tau*(matrix(1,1,N-P)%*%PFModMat[(PPlusM+1):(M+N),])
            }else{
                PFModMat <- rbind(-JMat, -KMat, -(EJMat[(P+1):(P+Pi),] + FMat[(P+1):(P+Pi),]%*%KMat), EJMat[(P+Pi+1):N,] + FMat[(P+Pi+1):N,]%*%KMat, JMat)
                if ((P > 1)|| !NoWeightI){MuR0Row <- (Tau*(matrix(1,1,Pi)%*%FMat[(P+1):(P+Pi),]) - (1-Tau)*(matrix(1,1,Nu)%*%FMat[(P+Pi+1):N,]))%*%LMat}
                MuBRow <- -Tau*(matrix(1,1,Pi)%*%PFModMat[(PPlusM+1):(PPlusM+Pi),]) - (1-Tau)*(matrix(1,1,Nu)%*%PFModMat[(PPlusM+Pi+1):(N+M),])
            }
            if ((P == 1) && NoWeightI) MuR0Row <- Tau*Pi - (1-Tau)*Nu

            #L2 normalizing of PFModMat rows
            if (PFNormI){
                AuxPFNormVec <- sqrt((PFModMat[(M+1):(N+M),]*PFModMat[(M+1):(N+M),])%*%OnesM1Vec)
                AuxPFNormVec <- AuxPFNormVec + (AuxPFNormVec == 0)
                PFModMat[(M+1):(N+M),] <- PFModMat[(M+1):(N+M),]/(AuxPFNormVec%*%Ones1MRow)
            }

            #testing if (PF) is satisfied (or if more postoptimization steps are required)
            Z1Vec <- -PFModMat[1:(N+M),]%*%FCVec   #the primal feasibility criterion vector z_1
            if (!all(Z1Vec > -EpsPF)){ #if (PF) fails
                MinZ1 <- min(Z1Vec); ArgMinZ1 <- which.min(Z1Vec)
                COutST$CompErrMsgS <- paste('The (PF) condition ceased to be satisfied : Z1Vec[',ArgMinZ1,'] = ',MinZ1,'!')
                return(COutST)
            }

            D23SubVec <- t(cbind(-MuR0Row - Tau*Ones1PRow, MuR0Row - (1-Tau)*Ones1PRow))

            #preparing the regression coefficient matrix ACOMat, AVec = ACOMat%*%u
            #-------------------------------------------------------------------

            ACOMat <- KMat
            ACOMat[BasNegAVec,] <- -ACOMat[BasNegAVec,]

            #finding IH (the index indicating the violated constraint, i in the text)
            #------------------------------------------------------------------------

            if (CTechST$D2SpecI){

                ThetaOld <- ThetaNew
                B11 <- JMat[1,1]; B22 <- JMat[2,2]

                #each relevant constraint is equipped with the 'right' angle Theta between its border line and the negative x semiaxis
                if ((B11 <= 0) && (B22 <= 0)){
                    IsRCVec  <- (PFModMat[1:(N+M),1]>0)&(PFModMat[1:(N+M),2]<=0)
                    ThetaVec <- atan(abs(PFModMat[IsRCVec,1]/PFModMat[IsRCVec,2]))
                    ThetaVec <- -pi + ThetaVec
                }else if ((B11 >= 0) && (B22 <= 0)){
                    IsRCVec  <- (PFModMat[1:(N+M),1]>=0)&(PFModMat[1:(N+M),2]>0)
                    ThetaVec <- atan(abs(PFModMat[IsRCVec,1]/PFModMat[IsRCVec,2]))
                    ThetaVec <- -ThetaVec
                }else if ((B11 >= 0) && (B22 >= 0)){
                    IsRCVec  <- (PFModMat[1:(N+M),1]<0)&(PFModMat[1:(N+M),2]>=0)
                    ThetaVec <- atan(abs(PFModMat[IsRCVec,1]/PFModMat[IsRCVec,2]))
                }else if ((B11 <= 0) && (B22 >= 0)){
                    IsRCVec  <- (PFModMat[1:(N+M),1]<=0)&(PFModMat[1:(N+M),2]<0)
                    ThetaVec <- atan(abs(PFModMat[IsRCVec,1]/PFModMat[IsRCVec,2]))
                    ThetaVec <- pi - ThetaVec
                }

                #identifying the tightest constraint
                RCIndVec <- OneToNPlusMVec[IsRCVec]
                ThetaVec <- ThetaVec + (2*pi)*(ThetaVec < U0Angle - EpsTheta)
                ThetaH <- min(ThetaVec + 100*(ThetaVec <= ThetaOld + EpsTheta)); IHAux <- which.min(ThetaVec + 100*(ThetaVec <= ThetaOld + EpsTheta))
                IHN <- RCIndVec[IHAux] #index of column to remove in the (N) representation
                FCVec <- rbind(PFModMat[IHN,2], -PFModMat[IHN,1])
                FCVec <- Dist*FCVec/norm(FCVec, type="2")   #the standardized center/interior point of the new facet

                if (ThetaH <= 50){
                    ThetaNew <- ThetaH
                    IHP <- TST$IC[IHN] #index of column to remove in the (P) representation
                    NNewFN <- 1        #the number of new facets from this cone (for the next step)
                    NARPresDat <- 1    #if 0, then the program terminates (successfully)
                }else{ #preparation for the end
                    ThetaNew <- U0Angle+2*pi
                    IHN <- Inf
                    IHP <- Inf
                    NNewFN <- 0
                    NARPresDat <- 0    #if 0, then the program terminates (successfully)
                }
                AngleStart <- ThetaOld - 2*pi*(ThetaOld>pi) #the angle of the old facet
                AngleEnd   <- ThetaNew - 2*pi*(ThetaNew>pi) #the angle of the new facet

                if (COutST$NumB == 1){ #at the very beginning
                    U0Angle <- AngleEnd #change U0Angle so that no artifical/false cone is stored
                    VVMat <- NULL #the vertices of the new adjacent cone to store
                    NRowVV <- 0   #the number of vertices to record/store
                }else{
                    #the vertices of the new adjacent cone to store
                    VVMat <- matrix(c(cos(AngleStart), cos(AngleEnd), sin(AngleStart), sin(AngleEnd)), 2, 2)
                    NRowVV <- 2   #the number of vertices to record/store
                }
            }else{



                #finding an interior point IP of the basic cone
                #----------------------------------------------


                for (i in 0:log2(IniStep/MinStep)){ #usually only 1 iteration is required
                    IPVec <- FCVec + (IniStep/2^i)*FNVec
                    IsIn  <- all(PFModMat[1:(N+M),1:M]%*%IPVec < -EpsInt)
                    #(the other constraints are always satisfied for Dist and IniStep small enough)
                    if (IsIn){break}
                }
                if (!IsIn){COutST$NSkipCone <- COutST$NSkipCone+1; next}
                IPVec <- Dist*IPVec/norm(IPVec, type="2")

                #finding non-redundant facets of the basic cone
                #----------------------------------------------

                #the polytope of our concern is given by PFModMat%*%UVec <= BBVec, BBVec = rbind(zeros[N+M,1], matrix(1,M,1))

                BBNewVec <- BBVec - PFModMat%*%IPVec
                DDMat <- PFModMat / (BBNewVec%*%Ones1MRow)
                KKMat <- convhulln(DDMat, QHullOptionsCA)
                AuxNRVec <- sort(as.vector(KKMat)[!duplicated(as.vector(KKMat))]) #indices of nonredundant constraints
                NRowKK <- dim(KKMat)[1]
                HHMat <- matrix(0,NRowKK,M)
                for (i in 1:NRowKK){HHMat[i,] <- solve(DDMat[KKMat[i,],],OnesM1Vec)}
                VVMat <- HHMat + matrix(1,NRowKK,1)%*%t(IPVec)

                VVMat <- VVMat[apply(abs(VVMat),1,max)>1e-3,]   #eliminating redundant zero vertices
                NRowVV <- dim(VVMat)[1]   #the number of vertices to record/store

                #non-redundant facets of our interest
                if (!CTechST$CubRegWiseI){NRFVec <- AuxNRVec[AuxNRVec <= (N+M)]
                }else{                    NRFVec <- AuxNRVec[(AuxNRVec <= (N+M))&(AuxNRVec > M)]
                }

                if (length(NRFVec)>0){

                    #looking for non-redundant facets not considered before
                    #------------------------------------------------------

                    #ActFNMat ... the matrix of potentially active facet normal (standardized) vectors (in rows)
                    if  (PFNormI){ActFNMat <- PFModMat[NRFVec,]
                    }else{ActFNMat <- PFModMat[NRFVec,]/(sqrt(PFModMat[NRFVec,]*PFModMat[NRFVec,])%*%OnesM1Vec%*%Ones1MRow)
                    }

                    #ActFCMat ... the matrix of standardized centers/interior points of possibly active facets (in rows)
                    ActFCMat <- (abs(ActFNMat%*%t(VVMat)) <= EpsFV)%*%VVMat
                    ActFCMat <- Dist*ActFCMat/sqrt(((ActFCMat*ActFCMat)%*%OnesM1Vec)%*%Ones1MRow)

                    #ActFIMat ... the matrix of facet identifiers for searching/storing in the archive
                    ActFIMat <- round(ActFCMat*FIPrec)/FIPrec
                    IsNewVec <- matrix(0,length(NRFVec),1) #the vector indicating which facet has not been visited yet
                    for (IndAF in 1:length(NRFVec)){
                        IDToFindRow <- ActFIMat[IndAF,]
                        if (CTechST$ArchAllFI){
                            OutList <- findRow(IDToFindRow, FIArchiveMat, NARFIArchive, M); IndFIArch <- OutList[[1]]; Status <- OutList[[2]]
                            if (Status != 1){ #if not found in the archive
                                #adding the facet to the archive of used ones
                                if (NARFIArchive >= MaxNF){
                                    COutST$CompErrMsgS <- 'The number of stored facet identifiers set by MaxNF or MaxNFIPerLayer exceeds the limit!'
                                    return(COutST)
                                }
                                FIArchiveMat <- addRow(IDToFindRow, FIArchiveMat, NARFIArchive, IndFIArch)
                                NARFIArchive <- NARFIArchive + 1
                                IsNewVec[IndAF] <- 1
                            }
                        }else{
                            IndFIVec <- matrix(0, NArchLayer,1)
                            for (IndL in 1:NArchLayer){
                                OutList <- findRow(IDToFindRow, FIArchMatCA[[IndL]], NARFIArchVec[IndL], M); IndFIVec[IndL] <- OutList[[1]]; Status <- OutList[[2]]
                                if (Status == 1){break}
                            }
                            if (Status != 1){ #if not found in the archive
                                if (NARFIArchVec[1] >= MaxNFIPerLayer){
                                    COutST$CompErrMsgS <- 'The number of stored facet identifiers set by MaxNF or MaxNFIPerLayer exceeds the limit!'
                                    return(COutST)
                                }
                                FIArchMatCA[[1]] <- addRow(IDToFindRow, FIArchMatCA[[1]], NARFIArchVec[1], IndFIVec[1])
                                NARFIArchVec[1]  <- NARFIArchVec[1] + 1
                                IsNewVec[IndAF]  <- 1
                            }
                        }
                    }
                    NNewFN <- sum(IsNewVec) #the number of new facets found
                    if (NNewFN > 0){IndNewVec <- which(IsNewVec == 1)} #IndNewVec ... indices of new facets
                    if (CTechST$SkipRedI && (NNewFN == 0)){next}
                }else{
                    NNewFN <- 0; NARPresDat <- 0
                } #if NRFVec empty
            } #if CTechST$D2SpecI

            #refreshing (and possibly storing) DQMat
            #---------------------------------------

            if (NRowVV > 0){
                ConeID <- ConeID+1

                if (NARDQ + NRowVV > MaxNRowQ){ #if DQMat is to overflow
                    COutST$NDQFiles <- COutST$NDQFiles + 1
                    COutST$CharST <- CTechST$getCharST(Tau, N, M, P, DQMat[1:NARDQ,1:(2+2*M+P*M)], COutST$CharST, 0)
                    if (CTechST$OutSaveI){
                        OutDQMat <- DQMat[1:NARDQ,]
                        ROFPath <- sprintf(DataCharF, COutST$NDQFiles)
                        write.table(OutDQMat, file=ROFPath, sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE)
                        rm(OutDQMat)
                    }
                    NARDQ <- 0
                }

                if (CTechST$BriefOutputI){
                    DQMat[(NARDQ+1):(NARDQ+NRowVV),] <- cbind(matrix(1,NRowVV,1)%*%c(ConeID, Nu), VVMat, matrix(1,NRowVV,1)%*%c(matrix(ACOMat,1,P*M), MuBRow))
                }else{
                    DQMat[(NARDQ+1):(NARDQ+NRowVV),] <- cbind(matrix(1,NRowVV,1)%*%c(ConeID, Nu), VVMat, matrix(1,NRowVV,1)%*%c(matrix(ACOMat,1,P*M), MuBRow, MuR0Row, t(TST$IZ)))
                }
                NARDQ <- NARDQ + NRowVV
                COutST$PosVec[TST$ITe] <- 2
                COutST$PosVec[TST$IZ] <- pmax(COutST$PosVec[TST$IZ], Ones1PRow)
            }

            #investigating new facets in turn
            if (NNewFN > 0){ #R can loop from 1 to 0
            for (IndNewF in 1:NNewFN){
                #finding the index I to leave the basis (if not known yet)
                #---------------------------------------------------------

                if (CTechST$D2SpecI == 0){
                    IHN <- NRFVec[IndNewVec[IndNewF]]
                    IHP <- TST$IC[IHN]
                }

                #finding the index JH to enter the basis
                #---------------------------------------

                if (IHN <= PPlusM){
                    JHN <- IHN
                    JHP <- TST$IC[N+M+JHN]
                    if ((IHN <= M) && CTechST$D2SpecI){
                        JMat[IHN,IHN] <- -JMat[IHN,IHN]
                    }
                }else{
                    if ((P > 1) || !NoWeightI){
                        FRowLRow <- FMat[IHN-M,]%*%LMat
                        if (IHN <= PPlusM+Pi){
                            T23SubVec <- t(cbind(-FRowLRow, FRowLRow))
                        }else{
                            T23SubVec <- t(cbind(FRowLRow, -FRowLRow))
                        }
                        IndNegTVec <- which(T23SubVec < 0)
                        if (!length(IndNegTVec)){
                            MinCV <- Inf
                        }else{
                            DFCritVec <- D23SubVec[IndNegTVec]/T23SubVec[IndNegTVec]
                            MinCV <- min(DFCritVec); IndMinCVAux <- which.min(DFCritVec)
                            IndMinCV <- IndNegTVec[IndMinCVAux]
                            if (MinCV < 0){
                                COutST$CompErrMsgS <- paste('Min {Dj/Tj, Tj < 0}, is negative:  j = ',PPlusM+IndMinCV,' Min = ',MinCV,'!')
                                return(COutST)
                            }
                        }
                    }else{ #i.e. if ((P == 1) && NoWeightI)
                        IndMinCV <- 2 - (IHN <= PPlusM+Pi)
                        MinCV <- -D23SubVec[IndMinCV]
                    }
                    if (MinCV <= 1){
                        JHP <- TST$IC[N+M+PPlusM+IndMinCV]
                    }else{
                        JHP <- TST$IC[N+M+2*P+IHN]
                    }
                }
















                CopyTST <- TST

                #updating CopyTST
                #----------------

                #removing IHP
                if (IHP > P2PlusM2+N){
                    AuxInd <- IHP-P2PlusM2-N
                    CopyTST$ITe  <- delItem(AuxInd, CopyTST$ITe)
                    CopyTST$IZ   <- addItem(AuxInd, CopyTST$IZ)
                    CopyTST$Nu   <- CopyTST$Nu-1

                }else if (IHP > P2PlusM2){
                    AuxInd  <- IHP-P2PlusM2
                    CopyTST$Ie   <- delItem(AuxInd, CopyTST$Ie)
                    CopyTST$IZ   <- addItem(AuxInd, CopyTST$IZ)
                    CopyTST$Pi   <- CopyTST$Pi-1

                }else if (IHP > 2*M){
                    AuxInd  <- IHP + ((IHP %% 2) == 1) - ((IHP %% 2) == 0)
                    CopyTST$Ia   <- delItem(IHP, CopyTST$Ia)
                    CopyTST$ITa  <- delItem(AuxInd, CopyTST$ITa)


                }else{ #(IHP <= 2*M)
                    AuxInd  <- IHP + ((IHP %% 2) == 1) - ((IHP %% 2) == 0)
                    CopyTST$Ib  <- delItem(IHP, CopyTST$Ib)
                    CopyTST$ITb <- delItem(AuxInd, CopyTST$ITb)


                }

                #adding JHP
                if (JHP > P2PlusM2+N){
                    AuxInd  <- JHP-P2PlusM2-N
                    CopyTST$ITe <- addItem(AuxInd, CopyTST$ITe)
                    CopyTST$IZ  <- delItem(AuxInd, CopyTST$IZ)
                    CopyTST$Nu  <- CopyTST$Nu+1

                }else if (JHP > P2PlusM2){
                    AuxInd  <- JHP-P2PlusM2
                    CopyTST$Ie  <- addItem(AuxInd, CopyTST$Ie)
                    CopyTST$IZ  <- delItem(AuxInd, CopyTST$IZ)
                    CopyTST$Pi  <- CopyTST$Pi+1

                }else if (JHP > 2*M){
                    AuxInd  <- JHP + ((JHP %% 2) == 1) - ((JHP %% 2) == 0)
                    CopyTST$Ia  <- addItem(JHP, CopyTST$Ia)
                    CopyTST$ITa <- addItem(AuxInd, CopyTST$ITa)


                }else{ #(JHP <= 2*M)
                    AuxInd  <- JHP + ((JHP %% 2) == 1) - ((JHP %% 2) == 0)
                    CopyTST$Ib  <- addItem(JHP, CopyTST$Ib)
                    CopyTST$ITb <- addItem(AuxInd, CopyTST$ITb)


                }

                #permuting the rows (if it is necessary)
                IsRowPerm <- (JHP > P2PlusM2) | (IHP > P2PlusM2)
                if (IsRowPerm){CopyTST$IR <- c(CopyTST$IZ, CopyTST$Ie, CopyTST$ITe)}

                #permuting the columns
                CopyTST$IC <- c(CopyTST$Ib, CopyTST$Ia, P2PlusM2+CopyTST$Ie, (P2PlusM2+N)+CopyTST$ITe,
                                CopyTST$ITb, CopyTST$ITa, P2PlusM2+CopyTST$IZ,
                                (P2PlusM2+N)+CopyTST$IZ, (P2PlusM2+N)+CopyTST$Ie, P2PlusM2+CopyTST$ITe)



















































































                if (CTechST$D2SpecI){
                    TST <- CopyTST
                }else{
                    #updating the list of bases to be investigated (corresponding to adjacent cones)
                    NARFutDat <- NARFutDat + 1; FutDatSA[[NARFutDat]] <- list()
                    FutDatSA[[NARFutDat]]$FNVec <- as.matrix(ActFNMat)[IndNewVec[IndNewF],]
                    FutDatSA[[NARFutDat]]$FCVec <- as.matrix(ActFCMat)[IndNewVec[IndNewF],]
                    FutDatSA[[NARFutDat]]$TST   <- CopyTST
                }
            } } #for IndNewF + if (NNewFN > 0)
        } #for FacetIndex

        #upgrade
        if (!CTechST$D2SpecI){
            NARPresDat  <- NARFutDat
            PresDatSA   <- FutDatSA
            NARFutDat   <- 0
            FutDatSA    <- list()
            if (!CTechST$ArchAllFI){
                FIArchMatCA[2:NArchLayer] <- FIArchMatCA[1:(NArchLayer-1)]
                NARFIArchVec[2:NArchLayer] <- NARFIArchVec[1:(NArchLayer-1)]
                NARFIArchVec[1] <- 0
            }
        }
        if (NARPresDat == 0){break}
    } #while (COutST$NumB >= CTechST$MaxB)
    if (COutST$NumB >= MaxNB){
        COutST$CompErrMsgS <- 'The number of bases exceeds the limit set by MaxNB!'
        return(COutST)
    }
} #for IndU0Vec

#saving all the rest
if ((COutST$NDQFiles == 0) || (NARDQ > 0)){
    COutST$NDQFiles <- COutST$NDQFiles + 1
    COutST$CharST <- CTechST$getCharST(Tau, N, M, P, DQMat[1:NARDQ,1:(2+2*M+P*M)], COutST$CharST, 0)
    if (CTechST$OutSaveI){
        OutDQMat <- DQMat[1:NARDQ,]
        ROFPath <- sprintf(DataCharF, COutST$NDQFiles)
        write.table(OutDQMat, file=ROFPath, sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE)
    }
}

#original setting of the warning states renewed
options(warn=WarnST)
return(COutST)
}