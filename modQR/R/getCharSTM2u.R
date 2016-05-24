getCharSTM2u <- function(Tau, N, M, P, BriefDQMat, CharST, IsFirst){
#getCharSTM2u <- function(Tau, N, M, P, BriefDQMat, CharST, IsFirst), output: CharST
#this function computes the quantile characteristics in the output of compContourM2u, namely
# in the list COutST$CharST. It makes possible to obtain some useful information without
#  saving any file on the disk, and it can be easily modified by the user according to his/her wishes.
#
#In the main program compContourM2u, this function is called first with BriefDQMat = NULL, CharST = NULL
# and with IsFirst = 1 to initialize CharST, and then with IsFirst = 0 successively
#  for the content of each potential output file, i.e. even if the output files are not stored on the disk
#   owing to CTechST$OutSaveI = 0. BriefDQMat then successively contains the rows of potential individual
#    output file(s) corresponding to CTechST$BriefOutputI = 1; see compContourM2u or the article for
#     further details.
#
#M ... the dimension of responses
#P ... the dimension of regressors including the intercept
#
#The default CharST has the following fields:
#always:
#  NUESkip ... the number of (skipped) directions artificially induced by the [-1, 1]^M bounding box
#  NAZSkip ... the number of (skipped) hyperplanes with at least one coordinate of a_1, ..., a_p zero
#  NBZSkip ... the number of (skipped) hyperplanes with at least one coordinate of b_1, ..., b_m zero
#
#if (M <= 4):
#  HypMat  ... the matrix with (slightly rounded) quantile hyperplane (b', a') coefficients (in M + P columns)
#                that can be used for the computation of the tau-contour
#if (P = 1):
#  CharMaxMat = [[u', max(Psi)]; [u', max(|MuBVec|)]]
#  CharMinMat = [[u', min(Psi)]; [u', min(|MuBVec|)]]
#
#and if (P > 1):
#CharMaxMat = [[u', max(Psi)]; [u', max(|MuBVec|)]]; [u', max(|(a_2, ..., a_P)'|)]; ...
#              [u', max(|a_2|)]; ...; [u', max(|a_P|)]]
#CharMinMat = [[u', min(Psi)]; [u', min(|MuBVec|)]]; [u', min(|(a_2, ..., a_P)'|)]]
#
#where
#(b', a') = (b_1, ..., b_M, a_1, ..., a_P)  ... the coefficients of quantile hyperplanes
#Psi (= u' * MuBVec)                        ... the optimal value of the objective function
#MuBVec                                     ... the Lagrange multipliers \mu_B corresponding to the constraint b = u

#- the max and min are computed over all quantile hyperplanes passing through M + P -1 observations and corresponding
#  cone vertices. The vector u in CharMaxMat and CharMinMat is (one of the) vectors where the row extreme is achieved.

if (IsFirst == 1){ 
    CharST <- list()
    CharST$CharMaxMat <- matrix(-Inf, 2+(P>1)+(P-1),M+1)
    CharST$CharMinMat <- matrix( Inf, 2+(P>1),M+1)
    CharST$NUESkip <- 0
    CharST$NAZSkip <- 0
    CharST$NBZSkip <- 0
    if (M <= 4){CharST$HypMat <- NULL}
}else{
    UMat <- BriefDQMat[,3:(2+M)]
    UMat <- UMat/matrix(apply(abs(UMat),1,max),dim(UMat)[1],M) #L1-normalization of the directions

    SelNEVec <- (rowSums((abs(UMat)-1) > -5e-9) == 1)          #for eliminating the directions artificially induced by the [-1,1]^M box
    CharST$NUESkip <- CharST$NUESkip + dim(UMat)[1] - sum(SelNEVec)

    UMat <- UMat[SelNEVec,]
    UMat <- UMat/matrix(sqrt((UMat^2)%*%matrix(1,M,1)),dim(UMat)[1],M) #L2-normalization of the directions
    AMat <- matrix(0,dim(UMat)[1],P)
    for (k in 1:M){AMat <- AMat + BriefDQMat[SelNEVec,(2+M+(k-1)*P+1):(2+M+k*P)]*(UMat[,k]%*%matrix(1,1,P))}
    MuBMat <- BriefDQMat[SelNEVec,(2+M+M*P+1):(2+M+M*P+M)]


    rm(BriefDQMat)


    SelBNZVec <- rowSums(abs(UMat) < 5e-9) == 0 #only those hyperplanes with all b_i coefficients nonzero, i = 1, ..., M
    SelANZVec <- rowSums(abs(AMat) < 5e-9) == 0 #only those hyperplanes with all a_j coefficients nonzero, j = 1, ..., P
    CharST$NBZSkip <- CharST$NBZSkip + sum(SelNEVec) - sum(SelBNZVec)
    CharST$NAZSkip <- CharST$NAZSkip + sum(SelNEVec) - sum(SelANZVec)

    SelVec <- SelBNZVec & SelANZVec
    UMat   <- UMat[SelVec,]
    AMat   <- AMat[SelVec,]
    MuBMat <- MuBMat[SelVec,]












    PsiVec <- (MuBMat*UMat)%*%matrix(1,M,1)
    MaxPsi <- max(PsiVec); Ind <- which.max(PsiVec)
    if (MaxPsi   > CharST$CharMaxMat[1,M+1]){         CharST$CharMaxMat[1,] <- c(UMat[Ind,], MaxPsi)}
    MinPsi <- min(PsiVec); Ind <- which.min(PsiVec)
    if (MinPsi   < CharST$CharMinMat[1,M+1]){         CharST$CharMinMat[1,] <- c(UMat[Ind,], MinPsi)}
    rm(PsiVec)

    NormMuBVec <- sqrt((MuBMat*MuBMat)%*%matrix(1,M,1))
    MaxNormMuB <- max(NormMuBVec); Ind <- which.max(NormMuBVec)
    if (MaxNormMuB  > CharST$CharMaxMat[2,M+1]){      CharST$CharMaxMat[2,] <- c(UMat[Ind,], MaxNormMuB)}
    MinNormMuB <- min(NormMuBVec); Ind <- which.min(NormMuBVec)
    if (MinNormMuB  < CharST$CharMinMat[2,M+1]){      CharST$CharMinMat[2,] <- c(UMat[Ind,], MinNormMuB)}
    rm(NormMuBVec)

    if (P > 1){
        A2PNormVec <- sqrt((AMat[,2:P]*AMat[,2:P])%*%matrix(1,P-1,1))
        MaxA2PNorm <- max(A2PNormVec); Ind <- which.max(A2PNormVec)
        if (MaxA2PNorm  > CharST$CharMaxMat[3,M+1]){      CharST$CharMaxMat[3,] <- c(UMat[Ind,], MaxA2PNorm)}
        MinA2PNorm <- min(A2PNormVec); Ind <- which.min(A2PNormVec)
        if (MinA2PNorm  < CharST$CharMinMat[3,M+1]){      CharST$CharMinMat[3,] <- c(UMat[Ind,], MinA2PNorm)}
        rm(A2PNormVec)








        for (i in 2:P){
            AINormVec <- abs(AMat[,i])
            MaxAINorm <- max(AINormVec); Ind <- which.max(AINormVec)
            if (MaxAINorm  > CharST$CharMaxMat[2+i,M+1]){     CharST$CharMaxMat[2+i,] <- c(UMat[Ind,], MaxAINorm)}
            rm(AINormVec)





        } #for i in 2:P
    } #if (P > 1)
    CharST$CharMinMat <- round(1e8*CharST$CharMinMat)/1e8
    CharST$CharMaxMat <- round(1e8*CharST$CharMaxMat)/1e8

    if (M <= 4){
        #rounding, sorting and omitting multiple entries
        BAMat <- round(1e8*cbind(UMat, AMat))/1e8
        CharST$HypMat <- unique(rbind(CharST$HypMat, BAMat))
        CharST$HypMat <- CharST$HypMat[do.call(order, split(CharST$HypMat, col(CharST$HypMat))),]
        colnames(CharST$HypMat) <- NULL
    }
} #if (IsFirst == 1)
return(CharST)
}