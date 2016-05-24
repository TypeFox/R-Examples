#Example A: interpreting the output
#illustrating the use of compContourM1u/compContourM2u and evalContour
# by means of their output
try(library('rgl'), silent=TRUE)
if (!exists('open3d')){
    print('OpenGL graphics and the rgl library not found. The output cannot be plotted. Try install.packages("rgl") first.')
}









set.seed(110)
N <- 7
M <- 2
P <- 2
Tau <- 0.2
YMat <- matrix(runif(N*M, -1, 1), N, M)
XMat <- cbind(matrix(1, N, 1), matrix(runif(N*(P-1), -1, 1), N, P-1))


print('Method No 1:')
CTechST <- getCTechSTM1u()
CTechST$BriefOutputI <- 0
COutST <- compContourM1u(Tau, YMat, XMat, CTechST); print(COutST)
if (!is.null(COutST$CharST$HypMat)){
    BBVec <- -COutST$CharST$HypMat[, M+1]
    AAMat <- cbind(-COutST$CharST$HypMat[, 1:M], COutST$CharST$HypMat[, (M+2):(M+P)])
    DataMat <- cbind(YMat, XMat[,2:P])
    IPVec <- as.matrix(colMeans(DataMat[(COutST$PosVec <= 1), ], 1))
    CST <- evalContour(AAMat, BBVec, IPVec); print(CST)
    if ((CST$Status == 0) && exists('open3d')){
        open3d(); rgl.viewpoint(-37.5, 30, 0); axes3d()
        title3d(xlab='y_1', ylab='y_2', zlab='x')
        wire3d(tmesh3d(t(CST$TVVMat), t(CST$TKKMat), homogeneous=FALSE))
        points3d(DataMat[, 1], DataMat[, 2], DataMat[, 3], size=16, color=c('blue'), alpha=0.5)
        points3d(CST$TVVMat[, 1], CST$TVVMat[, 2], CST$TVVMat[, 3], size=9, color=c('red'))

        #rgl.postscript('ExampleA1.pdf', fmt="pdf")
    }
}

print('Method No 2:')
CTechST <- getCTechSTM2u()
CTechST.BriefOutputI <- 0
COutST <- compContourM2u(Tau, YMat, XMat, CTechST); print(COutST)
if (!is.null(COutST$CharST$HypMat)){
    BBVec <- -COutST$CharST$HypMat[, M+1]
    AAMat <- cbind(-COutST$CharST$HypMat[, 1:M], COutST$CharST$HypMat[, (M+2):(M+P)])
    DataMat <- cbind(YMat, XMat[, 2:P])
    IPVec <- as.matrix(colMeans(DataMat[(COutST$PosVec <= 1), ], 1))
    CST <- evalContour(AAMat, BBVec, IPVec); print(CST)
    if ((CST$Status == 0) && exists('open3d')){
        open3d(); rgl.viewpoint(-37.5, 30, 0); axes3d()
        title3d(xlab='y_1', ylab='y_2', zlab='x')
        wire3d(tmesh3d(t(CST$TVVMat), t(CST$TKKMat), homogeneous=FALSE))
        points3d(DataMat[, 1], DataMat[, 2], DataMat[, 3], size=16, color=c('blue'), alpha=0.5)
        points3d(CST$TVVMat[, 1], CST$TVVMat[, 2], CST$TVVMat[, 3], size=9, color=c('red'))

        #rgl.postscript('ExampleA2.pdf', fmt="pdf")
    }
}

