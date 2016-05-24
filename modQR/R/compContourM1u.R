compContourM1u <- function(Tau = 0.2, YMat = NULL, XMat = NULL, CTechST = NULL){
#compContourM1u <- function(Tau = 0.2, YMat = NULL, XMat = NULL, CTechST = NULL), output: COutST
#computing (regression) quantile regions by the algorithm described in the CSDA paper:
#  Paindaveine, D. and \v{S}iman, M. (2012):
#  Computing Multiple-Output Regression Quantile Regions.
#  Computational Statistics & Data Analysis 56, 840--853.
#
#Tau      ... the quantile level in [0, 0.5]
#YMat     ... the response matrix with two to six columns
#XMat     ... the design matrix including the (first) intercept column
#CTechST  ... the list with some parameters regarding the computation, possibly modified by the user; see getCTechSTM1u
#
#COutST   ... the list containing some useful information about the computation, with the following fields:
#  CharST        ... the list with some default or user-defined output, provided by CTechST$getCharST (and initialized by getCharSTM1u)
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
#  PosVec        ... PosVec[i] = 0/1/2 if the i-th observation is in/on/out of the contour (after the successful computation)
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

if (is.null(CTechST)){CTechST <- getCTechSTM1u(); Status <- 0
}else{OutList <- checkCTechSTu(CTechST,1); CTechST <- OutList[[1]]; Status = OutList[[2]]}

if (Status == 1){COutST$CTechSTMsgS <- 'CTechST had to be replaced with the default list!'
} else if (Status > 1){COutST$CTechSTMsgS <- 'Some fields of CTechST had to be replaced with the default ones!'
}

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
BigProblemI <- FALSE
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
COutST$NDQFiles    <- 0      #the counter of (possible) output files, i.e. as if CTechST$OutSaveI = 1
COutST$NumB        <- 0      #the counter of optimal bases considered
COutST$MaxLWidth   <- 0      #the maximum width of one layer (if !CTechST$D2SpecI)
COutST$NIniNone    <- 0      #the number of trials when the initial solution(s) could not be found at all
COutST$NIniBad     <- 0      #the number of trials when the found initial solution(s) did not have the right number of nonzero coordinates
COutST$NSkipCone   <- 0      #the number of skipped cones (where the interior point could not be found)

#correcting the initial setting
if (M > 2){CTechST$D2SpecI   <- 0}
if (M > 3){CTechST$ArchAllFI <- 1}

#Initializing technical parameters
#=================================

Dist                   <- 0.9   #the distance of auxiliary interior points of the cones and of the facet 'centers' from the origin
DFNormI                <- 1     #1/0 ... if the D2345 vector in the postoptimization (DF condition: D2345 <= 0) should be normalized (1) or not (0)

#accuracy parameters
EpsTau                 <- 2e-6  #the perturbation of Tau to avoid (Tau*N) close to an integer; EpsTau should be less than 1/(2*N)
EpsDF                  <- 1e-7  #the bound for distinguishing zeros in the DF condition, in: JHNRow = which(DRow > (+ or -)EpsDF)+ AuxZeta + 1
if (CTechST$D2SpecI){
    EpsTheta <- 1e-15           #the minimum change of Theta (= the cone angle) to be considered significant
}else{
    if (M > 2){
        EpsFV          <- 5e-12 #the bound for identifying facet vertices (used for computing facet 'centers'), in: ActFCMat = (abs(ActFNMat%*%t(VVMat)) <= EpsFV)%*%VVMat
        EpsInt         <- 1e-11 #the bound distinguishing the points inside the cones, in: IsIn = all(QUModMat[1:(P2PlusM2-2+NQUApp),]%*%IPVec < -EpsInt)
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
NNInDQ                 <- 3e6   #the desired total number of numbers in DQMat (which is the matrix for keeping the results before storing them on the disk,
                                # processing them, or forgetting them). It should be at least 1000 or so.
MaxNB                  <- 4e7   #the maximum number of optimal bases allowed (when reached then the (faulty or too time-consuming) computation is interrupted)
MaxNPOIter             <- 500   #the maximum number of postoptimization iterations allowed (when reached then the mis-behaving computation is interrupted)
if (!CTechST$D2SpecI){
    if (!CTechST$ArchAllFI){
        NArchLayer       <- 12  #the number of the last layers whose cone facet identifiers are stored. It is wise to set it at least M+P.
        #MaxNFIPerLayer ... the maximum number of facet identifiers resulting from 1 layer (when reached then the (faulty or too memory-consuming) computation is interrupted)
        if (M == 2){MaxNFIPerLayer <- 10}else{MaxNFIPerLayer <- 2.5e4}
    }else{ #if CTechST$ArchAllFI
        MaxNF    <- 1e6         #the maximum number of used facet identifiers stored (when reached then the (faulty or too memory-consuming) computation is interrupted)
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
#It is highly recommended to read it side by side with the CSDA paper whose algorithm, notation, and text are followed

#perturbing Tau in the location case (to prevent degeneracy/multiple solutions)
if (P == 1){
    BadTau <- round(Tau*N)/N
    GapTau <- Tau - BadTau
    if (abs(GapTau) <= EpsTau){
        if (GapTau > - 1e-14){Tau <- BadTau + EpsTau }else{ Tau <- BadTau - EpsTau}
        COutST$TauMsgS <- sprintf('Tau has been perturbed to Tau = % 11.10f!\n', Tau)
        if (CTechST$ReportI){print(COutST$TauMsgS)}



    }
}

#auxiliary constants
PPlusM                <- P+M
P2PlusM2              <- 2*PPlusM
OnesDefZeta1Vec       <- matrix(1,PPlusM-1,1)
if (!CTechST$D2SpecI){
    Ones1MRow         <- matrix(1,1,M)
    OnesM1Vec         <- matrix(1,M,1)
    #BBVec ... the vector for describing the cones intersected by hypercubes, in: QUModMat%*%UVec <= BBVec
    if (CTechST$CubRegWiseI){BBVec <- rbind(matrix(0,P2PlusM2-2,1), matrix(0,M,1), matrix(1,M,1))
    }else{                   BBVec <- rbind(matrix(0,P2PlusM2-2,1), matrix(1,2*M,1))
    }
}

XYMat <- cbind(XMat, YMat)
rm(YMat); rm(XMat)

#preparing the matrix DQMat of size MaxNRowQ x NColQ for storing the results temporarily, before processing them at once with CTechST$getCharST
# and/or storing them on the disk. Then the matrix is refreshed and the process of recording the results continues.
if (CTechST$BriefOutputI){NColQ <- 1 + 1 + M + M + P + 1         #c(ConeID, Nu, UVec, BDVec, ADVec, LambdaD)
}else{NColQ <- 1 + 1 + M + M + P + 1 + (PPlusM-1)*M + (PPlusM-1) #c(ConeID, Nu, UVec, BDVec, ADVec, LambdaD, vec(VUMat), IZ)
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
    OutList <- findOptimalBasisM1FromScratch(M, N, P, Tau, U0Vec, XYMat); TST <- OutList[[1]]; U0Vec <- matrix(OutList[[2]],M,1); IsFound <- OutList[[3]]; TrialST <- OutList[[4]]
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

    #initializing storage variables (for monitoring visited cones and cones to be explored)
    #--------------------------------------------------------------------------------------

    if (CTechST$D2SpecI){
        ThetaNew <- U0Angle         #starting angle
        FCVec <- Dist*U0Vec         #standardized facet center/interior point (but not now)
    }else{
        #the data lists for storing some auxiliary results
        PresDatSA  <- list() #... from the topical layer
        FutDatSA   <- list() #... from the future layer
        NARFutDat  <- 0      #    the counter of active items in NARFutDat

        PresDatSA[[1]]        <- list() #initializing PresDatSA
        PresDatSA[[1]]$FNVec  <- matrix(0, M, 1) #the facet normal vector (but not now)
        PresDatSA[[1]]$FCVec  <- Dist*U0Vec      #the standardized facet center (interior point) u_f (but not now)
        PresDatSA[[1]]$TST    <- TST             #the corresponding TST list containing I_B, I_C, I_R and other useful characteristics

        #storing used facet identifiers
        if (CTechST$ArchAllFI){
            #a single array for storing identifiers of all visited facets
            FIArchiveMat  <- matrix(0, MaxNF, M)
            NARFIArchive <- 0     #the counter of active rows/items in FIArchiveMat
        }else{
            #a cell array for storing facet identifiers from the last NArchLayer layers
            #each cell corresponds to one layer
            FIArchMatCA  <- list()
            for (IndL in 1:NArchLayer){FIArchMatCA[[IndL]] <- matrix(0, MaxNFIPerLayer, M)}
            NARFIArchVec <- matrix(0, NArchLayer, 1)  #the vector of counters of active rows/items in FIArchMatCA
        }
    }

    #auxiliary counters
    NARPresDat  <- 1   #the number of facets investigated in the current step
    NumL        <- 0   #the number of layers (i.e. 'steps' of the breadth-first search algorithm)

    if (!CTechST$D2SpecI){JMat <- diag(as.vector(sign(U0Vec)))}

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

            }
            Pi   <- TST$Pi                         #the number of positive residuals
            Nu   <- TST$Nu                         #the number of negative residuals
            Zeta <- TST$Zeta                       #the number of zero residuals
            BasNegAVec <- ((TST$Ia %% 2) == 0)     #the indicators of basic a_{i-}'s
            BasPosAVec <- !BasNegAVec              #the indicators of basic a_{i+}'s
            BasNegBVec <- ((TST$Ib %% 2) == 0)     #the indicators of basic b_{i-}'s

            #finding auxiliary matrices (QUMat) and vectors
            #----------------------------------------------

            if (CTechST$D2SpecI && (COutST$NumB > 1)){
                A2NFirstColMat <- A2NStartMat
                RevSignIndVec  <- RevStartSignIndVec
            }else{
                A2NFirstColMat <- XYMat[TST$IR,]
                RevSignIndVec <- c(BasPosAVec, BasNegBVec)
                A2NFirstColMat[,RevSignIndVec] <- -A2NFirstColMat[,RevSignIndVec]
            }

            if ((P == 1) && (M == 2)){
                DetD <- A2NFirstColMat[1,2]*A2NFirstColMat[2,3] - A2NFirstColMat[1,3]*A2NFirstColMat[2,2]
                DIMat <- matrix(c(A2NFirstColMat[2,3], -A2NFirstColMat[2,2], -A2NFirstColMat[1,3], A2NFirstColMat[1,2]), 2, 2)/DetD
            }else{
                DIMat <- solve(A2NFirstColMat[1:(PPlusM-1),2:PPlusM])
            }
            SVec  <- -DIMat%*%A2NFirstColMat[1:(PPlusM-1),1]
            G0Col1Vec <- cbind(c(1,SVec))

            if (Pi == 0){
                HVec <- -(1-Tau)*t(matrix(1, 1, Nu)%*%A2NFirstColMat[PPlusM:N,])
                ZHMultiple   <- rbind(G0Col1Vec, -A2NFirstColMat[PPlusM:N,]%*%G0Col1Vec)
            }else if (Nu == 0){
                HVec <- Tau*t(matrix(1, 1, Pi)%*%A2NFirstColMat[PPlusM:N,])
                ZHMultiple   <- rbind(G0Col1Vec,  A2NFirstColMat[PPlusM:N,]%*%G0Col1Vec)
            }else{
                HVec <- Tau*t(matrix(1, 1, Pi)%*%A2NFirstColMat[PPlusM:(PPlusM - 1 + Pi),]) - (1-Tau)*t(matrix(1, 1, Nu)%*%A2NFirstColMat[(PPlusM + Pi):N,])
                ZHMultiple   <- rbind(G0Col1Vec,  A2NFirstColMat[PPlusM:(PPlusM - 1 + Pi),]%*%G0Col1Vec,-A2NFirstColMat[(PPlusM + Pi):N,]%*%G0Col1Vec)
            }

            VXMat <- matrix(0, PPlusM-1, M)
            for (i in 1:M){
                j <- P-1+i
                HModVec <- HVec; HModVec[P+i] <- 0
                VXMat[,i] <- t(rbind(-DIMat[j,], DIMat*SVec[j]-SVec%*%DIMat[j,]))%*%HModVec
            }
            VXModMat <- VXMat + OnesDefZeta1Vec%*%t(Tau*SVec[P:(PPlusM-1)])
            QXMat <- rbind(-VXModMat, VXModMat - (OnesDefZeta1Vec%*%t(SVec[P:(PPlusM-1)])))

            #preparing the output
            #--------------------

            LambdaD <- t(HVec)%*%G0Col1Vec
            ADVec <- G0Col1Vec[1:P];           ADVec[BasNegAVec]   <- - ADVec[BasNegAVec]
            BDVec <- G0Col1Vec[(P+1):PPlusM];  BDVec[BasNegBVec]   <- - BDVec[BasNegBVec]
            QUMat <- QXMat;                    QUMat[,BasNegBVec]  <- - QUMat[,BasNegBVec]
            if (!CTechST$BriefOutputI){
                VUMat <- VXMat;                VUMat[,BasNegBVec]  <- - VUMat[,BasNegBVec]
            }


            #finding JH (the index indicating the violated constraint, j in the text)
            #------------------------------------------------------------------------

            if (CTechST$D2SpecI){

                ThetaOld <- ThetaNew


                #each relevant constraint is equipped with the 'right' angle Theta between its border line and the negative x semiaxis
                ThetaVec     <- matrix(0, P2PlusM2-2, 1)
                DefThetaVec  <- atan(abs(QUMat[,1]/QUMat[,2]))
                IsZeroRowInQU <- 0
                for (i in 1:(P2PlusM2-2)){
                    QUI1 <- QUMat[i,1]; QUI2 <- QUMat[i,2]
                    if       ((QUI2 <= 0) && (QUI1 > 0)){ThetaVec[i] <- -pi  +DefThetaVec[i]
                    }else if ((QUI2 >= 0) && (QUI1 < 0)){ThetaVec[i] <-       DefThetaVec[i]
                    }else if ((QUI1 >= 0) && (QUI2 > 0)){ThetaVec[i] <-      -DefThetaVec[i]
                    }else if ((QUI1 <= 0) && (QUI2 < 0)){ThetaVec[i] <-  pi - DefThetaVec[i]
                    }else{   IsZeroRowInQU <- TRUE; break #i.e. if (QUI1 == 0 && QUI2 == 0)
                    }
                }
                if (IsZeroRowInQU){
                    COutST$CompErrMsgS <- paste('A 0/0 error appears during the computation of Theta, after Theta = ',ThetaOld)
                    return(COutST)
                }

                #identifying the tightest constraint

                ThetaVec <- ThetaVec + (2*pi)*(ThetaVec < U0Angle - EpsTheta)
                ThetaH <- min(ThetaVec + 100*(ThetaVec <= ThetaOld + EpsTheta)); JHAux <- which.min(ThetaVec + 100*(ThetaVec <= ThetaOld + EpsTheta))
                JHN <- PPlusM + JHAux   #the index of column to add in the (N) representation
                FCVec <- rbind(QUMat[JHAux,2], -QUMat[JHAux,1])
                FCVec <- Dist*FCVec/norm(FCVec, type="2")   #the standardized center/interior point of the new facet

                if (ThetaH <= 50){
                    ThetaNew <- ThetaH
                    JHP <- TST$IC[N+1+JHN] #the index of column to add in the (P) representation
                    NNewFN  <- 1           #the number of new facets from this cone (for the next step)
                    NARPresDat <- 1        #if 0, then the program terminates (successfully)
                }else{ #preparing for the end
                    ThetaNew <- U0Angle+2*pi
                    JHN <- Inf
                    JHP <- Inf
                    NNewFN <- 0
                    NARPresDat <- 0    #if 0, then the program terminates (successfully)
                }
                AngleStart <- ThetaOld - 2*pi*(ThetaOld>pi) #the angle of the old facet
                AngleEnd   <- ThetaNew - 2*pi*(ThetaNew>pi) #the angle of the new facet

                if (COutST$NumB == 1){ #at the very beginning
                    U0Angle <- AngleEnd  #change U0Angle so that no artifical/false cone is stored
                    VVMat <- NULL  #the vertices of the new adjacent cone to store
                    NRowVV <- 0  #the number of vertices to record/store
                }else{
                    #the vertices of the new adjacent cone to store
                    VVMat <- matrix(c(cos(AngleStart), cos(AngleEnd), sin(AngleStart), sin(AngleEnd)), 2, 2)
                    NRowVV <- 2 #the number of vertices to record/store
                }
            }else{
                #bounding and normalizing
                QUModMat <- rbind(QUMat/(sqrt((QUMat*QUMat)%*%OnesM1Vec)%*%Ones1MRow), -JMat, JMat)

                #finding an interior point IP of the basic cone
                #----------------------------------------------

                if (CTechST$CubRegWiseI){NQUApp <- M}else{NQUApp <- 0}
                for (i in 0:log2(IniStep/MinStep)){ #usually only 1 iteration is required
                    IPVec <- FCVec + (IniStep/2^i)*FNVec
                    IsIn  <- all(QUModMat[1:(P2PlusM2-2+NQUApp),]%*%IPVec < -EpsInt)
                    #(the other constraints are always satisfied for Dist and IniStep small enough)
                    if (IsIn){break}
                }
                if (!IsIn){COutST$NSkipCone <- COutST$NSkipCone+1; next}
                IPVec <- Dist*IPVec/norm(IPVec, type="2")

                #finding non-redundant facets of the basic cone
                #----------------------------------------------

                #the polytope of our concern is given by QUModMat%*%UVec <= BBVec

                BBNewVec <- BBVec - QUModMat%*%IPVec
                DDMat <- QUModMat / (BBNewVec%*%Ones1MRow)
                KKMat <- convhulln(DDMat, QHullOptionsCA)
                AuxNRVec <- sort(as.vector(KKMat)[!duplicated(as.vector(KKMat))]) #indices of nonredundant constraints
                NRowKK <- dim(KKMat)[1]
                HHMat  <- matrix(0,NRowKK,M)
                for (i in 1:NRowKK){HHMat[i,] <- solve(DDMat[KKMat[i,],],OnesM1Vec)}
                VVMat  <- HHMat + matrix(1,NRowKK,1)%*%t(IPVec)

                VVMat  <- VVMat[apply(abs(VVMat),1,max)>1e-3,]   #eliminating redundant zero vertices
                NRowVV  <- dim(VVMat)[1]   #the number of vertices to record/store

                #non-redundant facets of our interest
                NRFVec <- AuxNRVec[AuxNRVec <= (P2PlusM2-2)]



                if (length(NRFVec)>0){

                    #looking for non-redundant facets not considered before
                    #------------------------------------------------------

                    #ActFNMat ... the matrix of potentially active facet normal (standardized) vectors (in rows)
                    ActFNMat <- QUModMat[NRFVec,]



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
                                FIArchMatCA[[1]]  <- addRow(IDToFindRow, FIArchMatCA[[1]], NARFIArchVec[1], IndFIVec[1])
                                NARFIArchVec[1] <- NARFIArchVec[1] + 1
                                IsNewVec[IndAF] <- 1
                            }
                        }
                    }
                    NNewFN <- sum(IsNewVec) #the number of new facets found
                    if (NNewFN > 0){IndNewVec <- which(IsNewVec == 1)} #IndNewVec ... indices of new facets

                }else{
                    NNewFN <- 0; NARPresDat <- 0
                } #if !isempty(NRFVec)
            } #if CTechST$D2SpecI

            #refreshing (and possibly storing) DQMat
            #---------------------------------------

            if (NRowVV > 0){
                ConeID <- ConeID + 1

                if (NARDQ + NRowVV > MaxNRowQ){ #if DQMat is to overflow
                    COutST$NDQFiles <- COutST$NDQFiles + 1
                    COutST$CharST <- CTechST$getCharST(Tau, N, M, P, DQMat[1:NARDQ,1:(3+2*M+P)], COutST$CharST,0)
                    if (CTechST$OutSaveI){
                        OutDQMat <- DQMat[1:NARDQ,]
                        ROFPath <- sprintf(DataCharF, COutST$NDQFiles)
                        write.table(OutDQMat, file=ROFPath, sep="\t", append=TRUE, col.names=FALSE, row.names=FALSE)
                        rm(OutDQMat)
                    }
                    NARDQ <- 0
                }

                if (CTechST$BriefOutputI){
                    DQMat[(NARDQ+1):(NARDQ+NRowVV),] <- cbind(matrix(1,NRowVV,1)%*%c(ConeID, Nu), VVMat, matrix(1,NRowVV,1)%*%c(t(BDVec), t(ADVec), LambdaD))
                }else{
                    DQMat[(NARDQ+1):(NARDQ+NRowVV),] <- cbind(matrix(1,NRowVV,1)%*%c(ConeID, Nu), VVMat, matrix(1,NRowVV,1)%*%c(t(BDVec), t(ADVec), LambdaD, matrix(VUMat,1,(PPlusM-1)*M), t(TST$IZ)))
                }
                NARDQ <- NARDQ + NRowVV
                COutST$PosVec[TST$ITe] <- 2
                COutST$PosVec[TST$IZ] <- pmax(COutST$PosVec[TST$IZ], OnesDefZeta1Vec)
            }

            #investigating new facets in turn
            if (NNewFN > 0){ # R can loop from 1 to 0
            for (IndNewF in 1:NNewFN){
                #finding the index J to enter the basis (if not known yet)
                #---------------------------------------------------------

                if (!CTechST$D2SpecI){ 
                    JHN <- PPlusM + NRFVec[IndNewVec[IndNewF]]
                    JHP <- TST$IC[N+1+JHN]
                }

                CopyTST  <- TST
                if (!CTechST$D2SpecI){FCVec <- matrix(ActFCMat[IndNewVec[IndNewF],],M,1)}

                A1NStartRow <- cbind(matrix(0,1,P), t(FCVec))
                A1NStartRow[RevSignIndVec] <- -A1NStartRow[RevSignIndVec]
                A2NStartMat <- A2NFirstColMat
                AuxPi <- Pi; AuxNu <- Nu; AuxZeta <- Zeta

                #postoptimization
                #----------------
                Step   <- 0
                while (Step <= MaxNPOIter){
                    Step <- Step + 1

                    #if Step = 1 then AuxZeta = M+P-1, IsInANBar = 0, ZVec = ZHMultiple

                    #constructing Rho Vector
                    IsInANBar <- (JHN <= P2PlusM2-(AuxZeta+1))
                    if (!IsInANBar){
                        RefInd <- P2PlusM2-(AuxZeta+1)
                        AuxStartVec <- matrix(0, AuxZeta+1,1)
                        if (JHN <= RefInd+AuxZeta){AuxStartVec[JHN - RefInd + 1] <- -1
                        }else{                     AuxStartVec[JHN - RefInd - AuxZeta + 1] <- 1
                        }
                        Rho1Vec <- solve(rbind(A1NStartRow, A2NStartMat[1:AuxZeta,]), AuxStartVec)
                    }else{
                        Rho1Vec <- solve(rbind(A1NStartRow, A2NStartMat[1:AuxZeta,]), ANBarMat[1:(AuxZeta+1),JHN-(AuxZeta+1)])
                    }
                    if (AuxPi == 0){       RhoVec  <- c(Rho1Vec, -A2NStartMat[(AuxZeta+1):N,]%*%Rho1Vec)
                    }else if (AuxNu == 0){ RhoVec  <- c(Rho1Vec,  A2NStartMat[(AuxZeta+1):N,]%*%Rho1Vec)
                    }else{                 RhoVec  <- c(Rho1Vec,  A2NStartMat[(AuxZeta+1):(AuxZeta+AuxPi),]%*%Rho1Vec, -A2NStartMat[(AuxZeta+AuxPi+1):N,]%*%Rho1Vec)
                    }
                    if (IsInANBar){
                        if (AuxZeta+2+AuxPi>N+1){ RhoVec <- RhoVec + c(matrix(0,AuxZeta+1,1), -ANBarMat[(AuxZeta+2):(AuxZeta+1+AuxPi),JHN-(AuxZeta+1)])}else{RhoVec <- RhoVec + c(matrix(0,AuxZeta+1,1), -ANBarMat[(AuxZeta+2):(AuxZeta+1+AuxPi),JHN-(AuxZeta+1)], ANBarMat[(AuxZeta+2+AuxPi):(N+1),JHN-(AuxZeta+1)])}}
                    #constructing Z Vector
                    if (Step == 1){            ZVec <- ZHMultiple
                    }else{
                        Z1Vec <- C1IMat[,1]
                        if (AuxPi == 0){       ZVec    <- c(Z1Vec,   -A2NStartMat[(AuxZeta+1):N,]%*%Z1Vec)
                        }else if (AuxNu == 0){ ZVec    <- c(Z1Vec,    A2NStartMat[(AuxZeta+1):N,]%*%Z1Vec)
                        }else{                 ZVec    <- c(Z1Vec,    A2NStartMat[(AuxZeta+1):(AuxZeta+AuxPi),]%*%Z1Vec, -A2NStartMat[(AuxZeta+AuxPi+1):N,]%*%Z1Vec)
                        }
                    }

                    #evaluating the ZVec./RhoVec criterion
                    PosRhoIndVec <- which(RhoVec > 0)
                    if (length(PosRhoIndVec)==0){
                        MinRho <- max(RhoVec); ArgMinRho <- which.max(RhoVec)
                        COutST$CompErrMsgS <- paste('Max(RhoVec) is negative:  RhoVec[',ArgMinRho,'] = ', MinRho,'!')
                        return(COutST)
                    }else{
                        AuxZOverRhoVec <- ZVec[PosRhoIndVec]/RhoVec[PosRhoIndVec]
                        AuxInd <- which.min(AuxZOverRhoVec)
                        IHN <- PosRhoIndVec[AuxInd]
                        IHP <- CopyTST$IC[IHN]   #the index of original column to remove
                    }

                    #updating CopyTST
                    #----------------

                    #adding JHP
                    if (JHP > P2PlusM2 + N){
                        AuxInd  <- JHP-P2PlusM2-N
                        CopyTST$ITe  <- addItem(AuxInd, CopyTST$ITe)
                        CopyTST$IZ   <- delItem(AuxInd, CopyTST$IZ)
                        CopyTST$Nu   <- CopyTST$Nu+1
                        CopyTST$Zeta <- CopyTST$Zeta-1
                    }else if  (JHP > P2PlusM2){
                        AuxInd  <- JHP-P2PlusM2
                        CopyTST$Ie   <- addItem(AuxInd, CopyTST$Ie)
                        CopyTST$IZ   <- delItem(AuxInd, CopyTST$IZ)
                        CopyTST$Pi   <- CopyTST$Pi+1
                        CopyTST$Zeta <- CopyTST$Zeta-1
                    }else if (JHP > 2*P){
                        AuxInd  <- JHP + ((JHP %% 2) == 1) - ((JHP %% 2) == 0)
                        CopyTST$Ib   <- addItem(JHP, CopyTST$Ib)
                        CopyTST$ITb  <- addItem(AuxInd, CopyTST$ITb)
                        CopyTST$IHb  <- delItem(JHP, CopyTST$IHb)
                        CopyTST$IHb  <- delItem(AuxInd, CopyTST$IHb)
                    }else{ #(JHP <= 2*P)
                        AuxInd  <- JHP + ((JHP %% 2) == 1) - ((JHP %% 2) == 0)
                        CopyTST$Ia   <- addItem(JHP, CopyTST$Ia)
                        CopyTST$ITa  <- addItem(AuxInd, CopyTST$ITa)
                        CopyTST$IHa  <- delItem(JHP, CopyTST$IHa)
                        CopyTST$IHa  <- delItem(AuxInd, CopyTST$IHa)
                    }

                    #removing IHP
                    if (IHP > P2PlusM2+N){
                        AuxInd <- IHP-P2PlusM2-N
                        CopyTST$ITe  <- delItem(AuxInd, CopyTST$ITe)
                        CopyTST$IZ   <- addItem(AuxInd, CopyTST$IZ)
                        CopyTST$Nu   <- CopyTST$Nu-1
                        CopyTST$Zeta <- CopyTST$Zeta+1
                    }else if (IHP > P2PlusM2){
                        AuxInd  <- IHP-P2PlusM2
                        CopyTST$Ie   <- delItem(AuxInd, CopyTST$Ie)
                        CopyTST$IZ   <- addItem(AuxInd, CopyTST$IZ)
                        CopyTST$Pi   <- CopyTST$Pi-1
                        CopyTST$Zeta <- CopyTST$Zeta+1
                    }else if (IHP > 2*P){
                        AuxInd  <- IHP + ((IHP %% 2) == 1) - ((IHP %% 2) == 0)
                        CopyTST$Ib   <- delItem(IHP, CopyTST$Ib)
                        CopyTST$ITb  <- delItem(AuxInd, CopyTST$ITb)
                        CopyTST$IHb  <- addItem(IHP, CopyTST$IHb)
                        CopyTST$IHb  <- addItem(AuxInd, CopyTST$IHb)
                    }else{ #(IHP <= 2*P)
                        AuxInd  <- IHP + ((IHP %% 2) == 1) - ((IHP %% 2) == 0)
                        CopyTST$Ia   <- delItem(IHP, CopyTST$Ia)
                        CopyTST$ITa  <- delItem(AuxInd, CopyTST$ITa)
                        CopyTST$IHa  <- addItem(IHP, CopyTST$IHa)
                        CopyTST$IHa  <- addItem(AuxInd, CopyTST$IHa)
                    }

                    #permuting the rows (if it is necessary)
                    IsRowPerm <- (JHP > P2PlusM2) || (IHP > P2PlusM2)
                    if (IsRowPerm){CopyTST$IR <- cbind(c(CopyTST$IZ, CopyTST$Ie, CopyTST$ITe))}

                    #permuting the columns
                    CopyTST$IC <- cbind(c(CopyTST$Ia, CopyTST$Ib, P2PlusM2+CopyTST$Ie, (P2PlusM2+N)+CopyTST$ITe,
                                        CopyTST$ITa, CopyTST$ITb, CopyTST$IHa, CopyTST$IHb, P2PlusM2+CopyTST$IZ,
                                        (P2PlusM2+N)+CopyTST$IZ, (P2PlusM2+N)+CopyTST$Ie, P2PlusM2+CopyTST$ITe))
                    if (CopyTST$Zeta > PPlusM-1){COutST$CompErrMsgS <- 'Zeta is greater than M+P-1!'; return( COutST)}

                    #computing DVec
                    #--------------

                    AuxPi    <- CopyTST$Pi
                    AuxNu    <- CopyTST$Nu
                    AuxZeta  <- CopyTST$Zeta

                    A1NStartRow <- cbind(matrix(0,1,P), t(FCVec))
                    A2NStartMat <- XYMat[CopyTST$IR,]
                    RevStartSignIndVec  <- c((CopyTST$Ia %% 2) == 1, (CopyTST$Ib %% 2) == 0)
                    if (AuxZeta < PPlusM-1){
                        SelStartIndVec   <- c((CopyTST$Ia + ((CopyTST$Ia %% 2) == 1))/2, (CopyTST$Ib + ((CopyTST$Ib %% 2) == 1))/2)
                        SelBarIndVec     <- c((CopyTST$IHa + ((CopyTST$IHa %% 2) == 1))/2, (CopyTST$IHb + ((CopyTST$IHb %% 2) == 1))/2)
                        RevBarSignIndVec <- c((CopyTST$IHa %% 2) == 1, (CopyTST$IHb %% 2) == 0)
                        ANBarMat         <- rbind(A1NStartRow[SelBarIndVec], XYMat[CopyTST$IR, SelBarIndVec])
                        ANBarMat[,RevBarSignIndVec] <- -ANBarMat[,RevBarSignIndVec]
                        A1NStartRow      <- A1NStartRow[SelStartIndVec]
                        A2NStartMat      <- XYMat[CopyTST$IR, SelStartIndVec]
                    }
                    A1NStartRow[RevStartSignIndVec] <- -A1NStartRow[RevStartSignIndVec]
                    A2NStartMat[,RevStartSignIndVec] <- -A2NStartMat[,RevStartSignIndVec]

                    if (AuxZeta == 0){C1IMat <- 1/A1NStartRow[1]
                    }else{            C1IMat <- solve(rbind(A1NStartRow, A2NStartMat[1:AuxZeta,]))
                    }

                    if (AuxPi == 0){
                        MultStartRow <- -(1-Tau)*(matrix(1,1,AuxNu)%*%A2NStartMat[(AuxZeta+1):N,])%*%C1IMat
                        if (AuxZeta < PPlusM-1) AuxAddD2345Row <- (1-Tau)*(matrix(1,1,AuxNu)%*%ANBarMat[(AuxZeta+2):(N+1),])
                    }else if (AuxNu == 0){
                        MultStartRow <- Tau*(matrix(1,1,AuxPi)%*%A2NStartMat[(AuxZeta+1):N,])%*%C1IMat
                        if (AuxZeta < PPlusM-1) AuxAddD2345Row <- -Tau*(matrix(1,1,AuxPi)%*%ANBarMat[(AuxZeta+2):(N+1),])
                    }else{
                        MultStartRow <- Tau*(matrix(1,1,AuxPi)%*%A2NStartMat[(AuxZeta+1):(AuxZeta+AuxPi),])%*%C1IMat-(1-Tau)*(matrix(1,1,AuxNu)%*%A2NStartMat[(AuxZeta+AuxPi+1):N,])%*%C1IMat

                        if (AuxZeta < PPlusM-1){
                            AuxAddD2345Row <-  -Tau*(matrix(1,1,AuxPi)%*%ANBarMat[(AuxZeta+2):(AuxZeta+1+AuxPi),]) + (1-Tau)*(matrix(1,1,AuxNu)%*%ANBarMat[(AuxZeta+1+AuxPi+1):(N+1),])
                        }
                    }

                    #switch AuxZeta
                    if (AuxZeta == PPlusM-1){D2345Row <- c(-MultStartRow[2:dim(MultStartRow)[2]]-Tau, MultStartRow[2:dim(MultStartRow)[2]]-(1-Tau))
                    }else if (AuxZeta == 0){ D2345Row <- MultStartRow%*%ANBarMat[1:(AuxZeta+1),] + AuxAddD2345Row
                    }else{                   D2345Row <- c(MultStartRow%*%ANBarMat[1:(AuxZeta+1),] + AuxAddD2345Row, -MultStartRow[2:(AuxZeta+1)]-Tau, MultStartRow[2:(AuxZeta+1)]-(1-Tau))
                    }

                    if (DFNormI) D2345Row <- D2345Row/norm(D2345Row, type="2")

                    #evaluating the D2345Row criterion
                    #---------------------------------

                    JHNRow <- which(D2345Row > EpsDF) + AuxZeta + 1
                    LJHN   <- length(JHNRow)
                    if ((LJHN == 0) && (AuxZeta == PPlusM-1)){
                        break #postoptimization terminated successfully
                    }
                    if ((LJHN == 0) && (AuxZeta != PPlusM-1)){
                        JHNRow <- which(D2345Row > -EpsDF) + AuxZeta + 1
                        LJHN   <- length(JHNRow)
                        if (LJHN == 0){
                            COutST$CompErrMsgS <- 'The postoptimization ends with Zeta < M + P - 1!'
                            return(COutST)
                        }
                        JHPVec <- CopyTST$IC[N+1+JHNRow]
                        if ((LJHN == 1) && (JHPVec[1] == IHP)){
                            COutST$CompErrMsgS <- paste('The postoptimization wants to add and remove the same original index ', IHP)
                            return(COutST)
                        }
                        if (JHPVec[1] == IHP){JHN <- JHNRow[2]; JHP <- JHPVec[2]
                        }else{                JHN <- JHNRow[1]; JHP <- JHPVec[1]
                        }
                    }else{ #LJHN > 1
                        JHN <- JHNRow[1]
                        JHP <- CopyTST$IC[N+1+JHN]
                    }
                } #while (Step <= MaxNPOIter)
                if (Step > MaxNPOIter){
                    COutST$CompErrMsgS <- 'The number of postoptimization iterations exceeds the maximum set by MaxNPOIter!'
                   return( COutST)
                }

                if (CTechST$D2SpecI){
                    TST <- CopyTST
                }else{
                    #updating the list of bases to be investigated (corresponding to the adjacent cones)
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
    } #while (COutST$NumB <= MaxNB)
    if (COutST$NumB >= MaxNB){
        COutST$CompErrMsgS <- 'The number of bases exceeds the limit set by MaxNB!'
        return(COutST)
    }
} #for IndU0Vec

#saving all the rest
if ((COutST$NDQFiles == 0) || (NARDQ > 0)){
    COutST$NDQFiles <- COutST$NDQFiles + 1
    COutST$CharST <- CTechST$getCharST(Tau, N, M, P, DQMat[1:NARDQ,1:(3+2*M+P)], COutST$CharST, 0)
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