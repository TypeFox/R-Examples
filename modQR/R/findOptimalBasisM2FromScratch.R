findOptimalBasisM2FromScratch <- function(M, N, P, Tau, U0Vec, YMat, XMat){
#findOptimalBasisM2FromScratch <- function(M, N, P, Tau, U0Vec, YMat, XMat), output: list(TST, NewU0Vec, IsFound, TrialST)
#finding a nondegenerate initial optimal solution with nonzero finite coefficients a and b
#M, N, P, Tau ... m, n, p, \tau from the CS article
#input U0Vec is assumed to be a normalized vertex of the [-1,1]^M (hyper)cube where the computation should start
#YMat ... the response matrix NxM
#XMat ... the design matrix NxP
#IsFound      ... whether the optimal solution has been found (1) or not (0)
#output U0Vec ... the direction corresponding to the final optimal solution
#TrialST      the list providing some information about the search
#  TrialST$NNotFound ... the number of unsuccessful attempts to find any optimal solution
#  TrialST$NBad      ... the number of found optimal solutions not meeting the suitability conditions,
#                         i.e. with the wrong number of zero/infinite coordinates
#                          (within the prespecified numerical error given by NZBound)
#TST          the list containing the information about the optimal basis corresponding to the output U0Vec,
#              i.e. containing the vectors I_Z, I_e, IT_e, I_b, IT_b, I_a, IT_a, I_B, I_R, I_C
#               described in the article (IT ... \tilde{I})

TST <- list(); IsFound <- FALSE; TrialST <- list(); TrialST$NNotFound <- 0; TrialST$NBad <- 0

#setting the necessary technical parameters
NTrial <-  30;    #the maximum number of different choices of u_0 tried in case of troubles,
                  # i.e. the number of attempts to find the suitable optimal solution
NZBound <- 1e-10; #the bound/threshold for determining zero coordinates in the solution vector











#auxiliary constants and variables
P2PlusM2 <- 2*(P+M)
NewU0Vec <- U0Vec

#constructing vectors and matrices appearing in the linear representation of the problem (P):
# min CRVec'*ZHRVec subject to ARMat*ZHRVec = BRVec, ZHRVec >= 0
#(the original representation is modified/reduced by taking u for b)

#DenseMat <- [UMat, -VMat]; UMat and -VMat are submatrices of APMat described in the CS article
DenseMat = matrix(0,N,P2PlusM2)
for (j in 1:M){DenseMat[,2*j-1]     <-  YMat[,j];   DenseMat[,2*j]     <- -YMat[,j]};   rm( YMat)
for (j in 1:P){DenseMat[,2*M+2*j-1] <- -XMat[,j];   DenseMat[,2*M+2*j] <-  XMat[,j]};   rm( XMat)

CRVec <- rbind(matrix(0,2*P,1), Tau*matrix(1,N,1), (1-Tau)*matrix(1,N,1))
ARMat <- cbind(DenseMat[,(2*M+1):P2PlusM2], -diag(1,N), diag(1,N))

for (IndTrial in 1:NTrial){

    #a new choice of U0Vec in the same orthant (if required)
    if (IndTrial > 1){NewU0Vec <- U0Vec + 0.9*runif(1)*(runif(M,-0.5,0.5))/sqrt(M); NewU0Vec <- NewU0Vec/norm(NewU0Vec,type="2")}

    BetaVec <- matrix(0,2*M,1)
    for (j in 1:M){BetaVec[2*j-1] <- max(0,NewU0Vec[j]); BetaVec[2*j] <- max(0,-NewU0Vec[j])}
    BRVec <- -DenseMat[,1:(2*M)]%*%BetaVec

    #computing initial optimal solution ZHVec

    lpO <- lp("min", CRVec, ARMat, array('==', dim=c(N+1,1)), BRVec)
    ZHVec <- rbind(BetaVec, as.matrix(lpO$solution)) #the solution to the original problem (P)

    #checking for optimization errors
    if ((lpO$status > 0) || any(is.infinite(ZHVec))){TrialST$NNotFound <- TrialST$NNotFound + 1; next}
 
    #identifying basic variables (dealing with near-zeros)
    #A coefficients
    ACoefVec <- ZHVec[2*(M+1:P)-1] - ZHVec[2*(M+1:P)]
    for (j in 1:P){
        if (abs(ACoefVec[j]) < NZBound){ACoefVec[j] <- 0}
        ZHVec[2*M+2*j-1] <- max( ACoefVec[j],0)
        ZHVec[2*M+2*j]   <- max(-ACoefVec[j],0)
    } #for j







    #Residuals
    RVec = ZHVec[(P2PlusM2+1):(P2PlusM2+N)] - ZHVec[(P2PlusM2+N+1):(P2PlusM2+2*N)]
    for (j in 1:N){
        if (abs(RVec[j]) < NZBound){RVec[j] <- 0}
        ZHVec[P2PlusM2+j]   <- max( RVec[j],0)
        ZHVec[P2PlusM2+N+j] <- max(-RVec[j],0)
    } #for j

    sortStruct  <- sort(ZHVec, decreasing=TRUE, index.return=TRUE)
    if ((sortStruct$x[N+M+1] != 0)|(sortStruct$x[N+M] == 0)){TrialST$NBad <- TrialST$NBad + 1; next}

    #creating index vectors I_Z, I_e, IT_e, I_b, IT_b, I_a, IT_a, I_B, I_R, I_C
    IB <- sort(sortStruct$ix[1:(N+M)])
    TST$ITe  <- IB[(IB >= P2PlusM2+N+1) & (IB <= P2PlusM2+2*N)] - P2PlusM2 - N
      TST$Nu   <- length(TST$ITe)
    TST$Ie   <- IB[(IB >= P2PlusM2+1) & (IB <= P2PlusM2+N)] - P2PlusM2
      TST$Pi   <- length(TST$Ie)
    TST$IZ   <- setdiff(1:N,union(TST$Ie,TST$ITe))


    TST$Ia  <- IB[(IB >= 2*M+1) & (IB <= P2PlusM2)]
      TST$ITa <- TST$Ia + ((TST$Ia %% 2) == 1) - ((TST$Ia %% 2) == 0)

    TST$Ib  <- IB[(IB <= 2*M)]
      TST$ITb <- TST$Ib + ((TST$Ib %% 2) == 1) - ((TST$Ib %% 2) == 0)

    if ((length(TST$Ib) != M) || (length(TST$Ia) != P)){TrialST$NBad <- TrialST$NBad + 1; next}
    TST$IR  <- c(TST$IZ, TST$Ie, TST$ITe)
    TST$IC  <- c(IB, TST$ITb, TST$ITa, P2PlusM2+TST$IZ, (P2PlusM2+N)+TST$IZ, (P2PlusM2+N)+TST$Ie, P2PlusM2+TST$ITe)
    IsFound <- TRUE
    break
} #for IndTrial
U0Vec <- NewU0Vec
return(list(TST, NewU0Vec, IsFound, TrialST))
}