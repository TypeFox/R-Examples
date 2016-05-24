getCharSTM1u <- function(Tau, N, M, P, BriefDQMat, CharST, IsFirst){
#getCharSTM1u <- function(Tau, N, M, P, BriefDQMat, CharST, IsFirst), output: CharST
#this function computes the quantile characteristics in the output of compContourM1u, namely
# in the list COutST$CharST. It makes possible to obtain some useful information without
#  saving any file on the disk, and it can be easily modified by the user according to his/her wishes.
#
#In the main program compContourM1u, this function is called first with BriefDQMat = NULL, CharST = NULL
# and with IsFirst = 1 to initialize CharST, and then with IsFirst = 0 successively
#  for the content of each potential output file, i.e. even if the output files are not stored on the disk
#   owing to CTechST$OutSaveI = 0. BriefDQMat then successively contains the rows of potential individual
#    output file(s) corresponding to CTechST$BriefOutputI = 1; see compContourM1u or the article for
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
#  CharMaxMat = [[u', max(|b|)]; [u', max(Lambda)]; [u', max(Lambda/|b|)]]
#  CharMinMat = [[u', min(|b|)]; [u', min(Lambda)]; [u', min(Lambda/|b|)]]
#
#and if (P > 1):
#CharMaxMat = [[u', max(|b|)]; [u', max(Lambda)]; [u', max(Lambda/|b|)]; ...
#              [u', max(|(a_2, ..., a_P)'|)]; [u', max(|(a_2, ... ,a_P)'|/|b|)]; ...
#              [u', max(|a_2|)]; [u', max(|a_2|/|b|)]; ...; [u', max(|a_P|)]; [u', max(|a_P|/|b|)]]
#CharMinMat = [[u', min(|b|)]; [u', min(Lambda)]; [u', min(Lambda/|b|)]; ...
#              [u', min(|(a_2, ..., a_P)'|)]; [u', min(|(a_2, ..., a_P)'|/|b|)]]
#where
#(b', a') = (b_1, ..., b_M, a_1, ..., a_P) ... the coefficients of quantile hyperplanes
#Lambda                                    ... the Lagrange multiplier equal to the optimal value of the objective function
#
#- the max and min are computed over all quantile hyperplanes passing through M + P - 1 observations and corresponding
#- cone vertices. The vector u in CharMaxMat and CharMinMat is (one of the) vectors where the row extreme is achieved.

if (IsFirst == 1){
    CharST <- list()
    CharST$CharMaxMat <-  matrix(-Inf, 3+2*(P>1)+2*(P-1),M+1)
    CharST$CharMinMat <-  matrix( Inf, 3+2*(P>1),M+1)
    CharST$NUESkip <- 0
    CharST$NAZSkip <- 0
    CharST$NBZSkip <- 0
    if (M <= 4){CharST$HypMat <- NULL}
}else{
    UMat  <- BriefDQMat[,3:(2+M)]
    UMat  <- UMat/matrix(apply(abs(UMat),1,max),dim(UMat)[1],M) #L1-normalization of the directions

    SelNEVec <- (rowSums((abs(UMat)-1) > -5e-9) == 1)           #for eliminating the directions artificially induced by the [-1,1]^M box
    CharST$NUESkip <- CharST$NUESkip + dim(UMat)[1] - sum(SelNEVec)

    UMat       <- UMat[SelNEVec,]
    UMat       <- UMat/matrix(sqrt((UMat^2)%*%matrix(1,M,1)),dim(UMat)[1],M) #L2-normalization of the directions
    BDMat      <- BriefDQMat[SelNEVec,(3+M):(2+2*M)]
    ADMat      <- as.matrix(BriefDQMat[SelNEVec,(3+2*M):(2+2*M+P)])
    #ConeIDVec  <- BriefDQMat(SelNEVec,1)
    LambdaDVec <- as.matrix(BriefDQMat[SelNEVec,2*M+P+3])
    BDNormVec  <- sqrt((BDMat^2)%*%matrix(1,M,1))
    rm(BriefDQMat)


    SelBNZVec <- rowSums(abs(BDMat)/(BDNormVec%*%matrix(1,1,M)) < 5e-9) == 0 #only those hyperplanes with all b_i coefficients nonzero, i = 1, ..., M
    SelANZVec <- rowSums(abs(ADMat)/(BDNormVec%*%matrix(1,1,P)) < 5e-9) == 0 #only those hyperplanes with all a_j coefficients nonzero, j = 1, ..., P
    CharST$NBZSkip <- CharST$NBZSkip + sum(SelNEVec) - sum(SelBNZVec)
    CharST$NAZSkip <- CharST$NAZSkip + sum(SelNEVec) - sum(SelANZVec)

    SelVec     <- SelBNZVec & SelANZVec
    UMat       <- UMat[SelVec,]
    BDMat      <- BDMat[SelVec,]
    ADMat      <- ADMat[SelVec,]
    #ConeIDVec  <- ConeIDVec[SelVec,]
    LambdaDVec <- LambdaDVec[SelVec,]
    NominVec   <- (BDMat*UMat)%*%matrix(1,M,1)
    BDNormVec  <- sqrt((BDMat^2)%*%matrix(1,M,1))

    BNormVec         <- BDNormVec/NominVec
    MaxBNorm          <- max(BNormVec); Ind <- which.max(BNormVec)
    if (MaxBNorm    > CharST$CharMaxMat[1,M+1]){CharST$CharMaxMat[1,] <- c(UMat[Ind,], MaxBNorm)}
    MinBNorm          <- min(BNormVec); Ind <- which.min(BNormVec)
    if (MinBNorm    < CharST$CharMinMat[1,M+1]){CharST$CharMinMat[1,] <- c(UMat[Ind,], MinBNorm)}
    rm(BNormVec)

    LambdaVec        <- LambdaDVec/NominVec
    MaxLambda         <- max(LambdaVec); Ind <- which.max(LambdaVec)
    if (MaxLambda   > CharST$CharMaxMat[2,M+1]){CharST$CharMaxMat[2,] <- c(UMat[Ind,], MaxLambda)}
    MinLambda         <- min(LambdaVec); Ind <- which.min(LambdaVec)
    if (MinLambda   < CharST$CharMinMat[2,M+1]){CharST$CharMinMat[2,] <- c(UMat[Ind,], MinLambda)}
    rm(LambdaVec)

    LambdaFVec       <- LambdaDVec/BDNormVec
    MaxLambdaF        <- max(LambdaFVec); Ind <- which.max(LambdaFVec)
    if (MaxLambdaF  > CharST$CharMaxMat[3,M+1]){CharST$CharMaxMat[3,] <- c(UMat[Ind,], MaxLambdaF)}
    MinLambdaF        <- min(LambdaFVec); Ind <- which.min(LambdaFVec)
    if (MinLambdaF  < CharST$CharMinMat[3,M+1]){CharST$CharMinMat[3,] <- c(UMat[Ind,], MinLambdaF)}
    rm(LambdaFVec)

    if (P > 1){
        A2PNormVec       <- sqrt((ADMat[,2:P]^2)%*%matrix(1,P-1,1))/NominVec
        MaxA2PNorm        <- max(A2PNormVec); Ind <- which.max(A2PNormVec)
        if (MaxA2PNorm  > CharST$CharMaxMat[4,M+1]){CharST$CharMaxMat[4,] <- c(UMat[Ind,], MaxA2PNorm)}
        MinA2PNorm        <- min(A2PNormVec); Ind <- which.min(A2PNormVec)
        if (MinA2PNorm  < CharST$CharMinMat[4,M+1]){CharST$CharMinMat[4,] <- c(UMat[Ind,], MinA2PNorm)}
        rm(A2PNormVec)

        A2PNormFVec      <- sqrt((ADMat[,2:P]^2)%*%matrix(1,P-1,1))/BDNormVec
        MaxA2PNormF       <- max(A2PNormFVec); Ind <- which.max(A2PNormFVec)
        if (MaxA2PNormF > CharST$CharMaxMat[5,M+1]){CharST$CharMaxMat[5,] <- c(UMat[Ind,], MaxA2PNormF)}
        MinA2PNormF       <- min(A2PNormFVec); Ind <- which.min(A2PNormFVec)
        if (MinA2PNormF < CharST$CharMinMat[5,M+1]){CharST$CharMinMat[5,] <- c(UMat[Ind,], MinA2PNormF)}
        rm(A2PNormFVec)

        for (i in 2:P){
            AINormVec         <- abs(ADMat[,i])/NominVec
            MaxAINorm         <- max(AINormVec); Ind <- which.max(AINormVec)
            if (MaxAINorm  > CharST$CharMaxMat[3+2*i-1,M+1]){CharST$CharMaxMat[3+2*i-1,] <- c(UMat[Ind,], MaxAINorm)}
            rm(AINormVec)

            AINormFVec       <- abs(ADMat[,i])/BDNormVec
            MaxAINormF       <- max(AINormFVec); Ind <- which.max(AINormFVec)
            if (MaxAINormF > CharST$CharMaxMat[3+2*i,M+1]){CharST$CharMaxMat[3+2*i,] <- c(UMat[Ind,], MaxAINormF)}
            rm(AINormFVec)
        } #for i in 2:P
    } #if (P > 1)
    CharST$CharMinMat <- round(1e8*CharST$CharMinMat)/1e8
    CharST$CharMaxMat <- round(1e8*CharST$CharMaxMat)/1e8

    if (M <= 4){
        #rounding, normalizing, sorting and omitting multiple entries
        BAMat <- round(1e8*cbind(BDMat, ADMat)/(BDNormVec%*%matrix(1,1,M+P)))/1e8
        CharST$HypMat <- unique(rbind(CharST$HypMat, BAMat))
        CharST$HypMat <- CharST$HypMat[do.call(order, split(CharST$HypMat, col(CharST$HypMat))),]
        colnames(CharST$HypMat) <- NULL
    }
} #if (IsFirst == 1)
return(CharST)
}