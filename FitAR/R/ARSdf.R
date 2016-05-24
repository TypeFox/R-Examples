`ARSdf` <-
function(phi, pFFT=8){
pext<-c(1,-phi,rep(0,2^(1+pFFT) -1 -length(phi)))
ft<-fft(pext)
1/(Re(ft*Conj(ft)))[1:2^pFFT]
}

