Fmatrix <-
function(t1,t2,f0,Sx2,Sy2,Pix,Piy){
F0<-diag(f0)
F1<-t(Pt(Sx2,Pix,t1))%*%F0%*%Pt(Sy2,Piy,t2)
F1
}
