findOptimalBasisM1FromScratch <-function(M, N, P, Tau, U0Vec, XYMat){
#findOptimalBasisM1FromScratch <-function(M, N, P, Tau, U0Vec, XYMat), output: list(TST, U0Vec, IsFound, TrialST)
#finding a nondegenerate initial optimal solution with nonzero finite coefficients a and b
#M, N, P, Tau ... m, n, p, \tau from the CSDA article
#input U0Vec is assumed to be a normalized vertex of the [-1,1]^M (hyper)cube where the computation should start
#XYMat        = [XMat, YMat], (XMat ... the design matrix NxP, YMat ... the response matrix NxM)
#
#IsFound      ... whether the optimal solution has been found (1) or not (0)
#output U0Vec ... the direction corresponding to the final optimal solution
#TrialST      the list providing some information about the search
#  TrialST$NNotFound ... the number of unsuccessful attempts to find any optimal solution
#  TrialST$NBad      ... the number of found optimal solutions not meeting the suitability conditions,
#                         i.e. with the wrong number of zero/infinite coordinates
#                          (within the prespecified numerical error given by NZBound)
#TST          the list containing the information about the optimal basis corresponding to the output U0Vec,
#              i.e. containing the vectors I_Z, I_e, IT_e, I_b, IT_b, IH_b, I_a, IT_a, IH_a, I_B, I_R, I_C
#               described in the article (IT ... \tilde{I}, IH ... \hat{I})

TST <- list(); IsFound <- FALSE; TrialST <- list(); TrialST$NNotFound <- 0; TrialST$NBad <- 0

#setting the necessary technical parameters
NTrial <-  30;    #the maximum number of different choices of u_0 tried in case of troubles,
                  # i.e. the number of attempts to find the suitable optimal solution
NZBound <- 1e-10; #the bound/threshold for determining zero coordinates in the solution vector











#auxiliary constants and variables
P2PlusM2 <- 2*(P+M)
NewU0Vec <- U0Vec

#constructing vectors and matrices appearing in the linear representation of the problem (P):
# min CPVec'*ZHPVec subject to APMat*ZHPVec = BPVec, ZHPVec >= 0


#DenseMat <- cbind(-VMat, UMat); -VMat and UMat are submatrices of APMat described in the CSDA article
DenseMat <- matrix(0, N,P2PlusM2)
for (j in 1:P){DenseMat[ ,    2*j-1] <- -XYMat[,  j];  DenseMat[ ,    2*j] <-  XYMat[ ,  j]}
for (j in 1:M){DenseMat[ ,2*P+2*j-1] <-  XYMat[,P+j];  DenseMat[ ,2*P+2*j] <- -XYMat[ ,P+j]}
rm(XYMat)
CPVec <- rbind(array(0, dim=c(P2PlusM2,1)), array(Tau, dim=c(N,1)), array(1-Tau, dim=c(N,1)))
BPVec <- array(0, dim=c(N+1,1)); BPVec[1,] <- 1

for (IndTrial in (1:NTrial)){

    #a new choice of U0Vec in the same orthant (if required)
    if (IndTrial > 1){NewU0Vec <- U0Vec + 0.9*runif(1)*(runif(M,-0.5,0.5))/sqrt(M); NewU0Vec <- NewU0Vec/norm(NewU0Vec, type="2")}

    A1PMat <- array(0, dim=c(1,P2PlusM2+2*N))
    for (j in 1:M) {A1PMat[1, 2*P+2*j-1] <- NewU0Vec[j]; A1PMat[1, 2*P+2*j] <- -NewU0Vec[j]}
    APMat <- rbind(A1PMat, cbind(DenseMat, -diag(N), diag(N)))

    #computing initial optimal solution ZHPVec

    lpO <- lp("min", CPVec, APMat, array('==', dim=c(N+1,1)), BPVec)
    ZHPVec <- lpO$solution #the solution to the original problem (P)

    #checking for optimization errors
    if ((lpO$status > 0) || any(is.infinite(ZHPVec))){TrialST$NNotFound <- TrialST$NNotFound + 1; next}

    #identifying basic variables (dealing with near-zeros)
    #A coefficients
    ACoefVec <- ZHPVec[2*(1:P)-1] - ZHPVec[2*(1:P)]
    for (j in 1:P){
        if (abs(ACoefVec[j]) < NZBound){ACoefVec[j] <- 0}
        ZHPVec[2*j-1] <- max( ACoefVec[j],0)
        ZHPVec[2*j]   <- max(-ACoefVec[j],0)
    } #end for j
    #B coefficients
    BCoefVec <- ZHPVec[2*((P+1):(P+M))-1] - ZHPVec[2*((P+1):(P+M))]
    for (j in 1:M){
        if (abs(BCoefVec[j]) < NZBound){BCoefVec[j] <- 0}
        ZHPVec[2*P+2*j-1] <- max( BCoefVec[j],0)
        ZHPVec[2*P+2*j]   <- max(-BCoefVec[j],0)
    } #for
    #Residuals
    RVec <- ZHPVec[P2PlusM2+(1:N)] - ZHPVec[P2PlusM2+N+(1:N)]
    for (j in 1:N){
        if (abs(RVec[j]) < NZBound){RVec[j] <- 0}
        ZHPVec[P2PlusM2+j]   <- max( RVec[j],0)
        ZHPVec[P2PlusM2+N+j] <- max(-RVec[j],0)
    } #for j

    sortStruct  <- sort(ZHPVec, decreasing=TRUE, index.return=TRUE)
    if ((sortStruct$x[N+2] != 0)||(sortStruct$x[N+1]==0)){TrialST$NBad <- TrialST$NBad + 1; next}

    #creating index vectors I_Z, I_e, IT_e, I_b, IT_b, IH_b, I_a, IT_a, IH_a, I_B, I_R, I_C
    IB <- sort(sortStruct$ix[1:(N+1)])
    TST$ITe  <- IB[(IB >= P2PlusM2+N+1) & (IB <= P2PlusM2+2*N)] - P2PlusM2 - N
      TST$Nu   <- length(TST$ITe)
    TST$Ie   <- IB[(IB >= P2PlusM2+1) & (IB <= P2PlusM2+N)] - P2PlusM2
      TST$Pi   <- length(TST$Ie)
    TST$IZ   <- setdiff(1:N,union(TST$Ie,TST$ITe))

      TST$Zeta <- length(TST$IZ)
    TST$Ia  <- IB[IB <= 2*P]
      TST$ITa <- TST$Ia + ((TST$Ia %% 2) == 1) - ((TST$Ia %% 2) == 0)
      TST$IHa <- NULL
    TST$Ib  <- IB[(IB >= 2*P+1) & (IB <= P2PlusM2)]
      TST$ITb <- TST$Ib + ((TST$Ib %% 2) == 1) - ((TST$Ib %% 2) == 0)
      TST$IHb  <- NULL
    if ((length(TST$Ib) != M) || (length(TST$Ia) != P) || (TST$Zeta != P+M-1)){TrialST$NBad <- TrialST$NBad + 1; next}
    TST$IR  <- c(TST$IZ, TST$Ie, TST$ITe)
    TST$IC  <- c(IB, TST$ITa, TST$ITb, TST$IHa, TST$IHb, P2PlusM2+TST$IZ, (P2PlusM2+N)+TST$IZ, (P2PlusM2+N)+TST$Ie, P2PlusM2+TST$ITe)
    IsFound <- TRUE
    break
} #for IndTrial
U0Vec <- NewU0Vec
return(list(TST, NewU0Vec, IsFound, TrialST))
}