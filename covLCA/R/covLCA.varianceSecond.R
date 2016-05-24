covLCA.varianceSecond <-
function(rgivy,x,S1,R,J,K.j,probs,y,S2,N,z)
{
	BB=covLCA.covLBetaBeta(rgivy,x,S1,R)
	GG=covLCA.covLGammaGamma(J,R,K.j,probs,rgivy,y)
	AA=covLCA.covLAlphaAlpha(J,S2,K.j,N,R,probs,rgivy,z)
	AG=covLCA.covLAlphaGamma(J,S2,K.j,R,z,y,probs,N,rgivy)
	BG=covLCA.covLBetaGamma(S1,R,J,K.j,x,y,probs,rgivy)
	BA=covLCA.covLBetaAlpha(S1,R,J,S2,K.j,x,z,probs,rgivy,N)
	
	ret=rbind(cbind(BB,BG,BA),cbind(t(BG),GG,t(AG)),cbind(t(BA),AG,AA))
	
	return(ret)
}
