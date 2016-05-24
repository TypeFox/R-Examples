#Example C: 3D location quantiles
#illustrating the use of compContourM1u/compContourM2u and evalContour
# on 3D location quantiles
try(library('rgl'), silent=TRUE)
if (!exists('open3d')){
    print('OpenGL graphics and the rgl library not found. The output cannot be plotted. Try install.packages("rgl") first.')
    return;
}








N <- 21; M <- 3; Tau <- 0.1


open3d(); rgl.viewpoint(37.5, 30, 0)
title3d(xlab='y_1', ylab='y_2', zlab='y_3'); axes3d()

set.seed(8416); YMat <- matrix(runif(N*M), N, M); XVec <- matrix(1, N, 1)

#Way 1:
COutST <- compContourM1u(Tau, YMat, XVec)
if (length(COutST$CharST$HypMat)){
    CST <- evalContour(-COutST$CharST$HypMat[, 1:M], -COutST$CharST$HypMat[, (M+1)])
    if (CST$Status == 0){
        wire3d(tmesh3d(t(CST$TVVMat), t(CST$TKKMat), homogeneous=FALSE))

    }
}
IsInVec <- (COutST$PosVec == 0)
points3d(YMat[IsInVec, 1], YMat[IsInVec, 2], YMat[IsInVec, 3], size=9, color=c('red'))
IsOnVec <- (COutST$PosVec == 1)
points3d(YMat[IsOnVec, 1], YMat[IsOnVec, 2], YMat[IsOnVec, 3], size=9, color=c('green'))
IsOutVec <- (COutST$PosVec == 2)
points3d(YMat[IsOutVec, 1], YMat[IsOutVec, 2], YMat[IsOutVec, 3], size=9, color=c('blue'))

#Way 2:
ScaleRow <- apply(YMat, 2, max)
CenterRow <- apply(YMat, 2, mean)
YMatAux <- t((t(YMat)-CenterRow)/ScaleRow) #centering and scaling
COutST <- compContourM1u(Tau, YMatAux, XVec)
if (length(COutST$CharST$HypMat)){
    CST <- evalContour(-COutST$CharST$HypMat[, 1:M], -COutST$CharST$HypMat[, (M+1)])
    if (CST$Status == 0){
        TVVMat <- t((t(CST$TVVMat)*ScaleRow) + CenterRow) #decentering and descaling
        wire3d(tmesh3d(t(TVVMat), t(CST$TKKMat), homogeneous=FALSE))

    }
}
IsInVec <- (COutST$PosVec == 0)
points3d(YMat[IsInVec, 1], YMat[IsInVec, 2], YMat[IsInVec, 3], size=5, color=c('red'))
IsOnVec <- (COutST$PosVec == 1)
points3d(YMat[IsOnVec, 1], YMat[IsOnVec, 2], YMat[IsOnVec, 3], size=5, color=c('green'))
IsOutVec <- (COutST$PosVec == 2)
points3d(YMat[IsOutVec, 1], YMat[IsOutVec, 2], YMat[IsOutVec, 3], size=5, color=c('blue'))




