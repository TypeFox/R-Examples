"SimGLP" <-
function(psi,a){
r<-convolve(a, rev(psi),type="o")
Q<-length(psi)-1
r[-(1:Q)][1:(length(a)-Q)]
}

