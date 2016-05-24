cmat.star <-
function(plist, CorrMat,  no.ord, no.norm){


if (no.norm==0) {
Sigma = IntermediateOO(plist, CorrMat)
}

if (no.ord==0) {
Sigma = CorrMat
}

if (no.norm>0 & no.ord>0) {
if ( validate.target.cormat(plist, CorrMat,  no.ord, no.norm)) { ## target correlation and provabilities are validated here.
k  = length(plist)
OO = IntermediateOO(plist, CorrMat[1:k,1:k])
ON = IntermediateON(plist, CorrMat[(k+1):nrow(CorrMat), 1:k] )
NN = CorrMat[(k+1):ncol(CorrMat), (k+1):ncol(CorrMat) ]
Sigma = cbind(rbind(OO,ON), rbind(t(ON),NN) )

if(!is.positive.definite(Sigma)){
warning( "Intermediate correlation matrix is not positive definite. A nearPD function is applied.")
Sigma=as.matrix(nearPD(Sigma, corr = TRUE, keepDiag = TRUE)$mat)
}
Sigma = ( Sigma+t(Sigma) )/2
}
}
return(Sigma)

}
