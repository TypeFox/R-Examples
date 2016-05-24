`freq_cukier` <-
function(m, i=1, omega_before=0){
if(i<=1){

 #to reproduce table VI in cukier75
 #cbind(omega_m_cukier75,c(d_m_cukier75,NA),min_runs_cukier75*2)
if(m>NROW(omega_m_cukier75)){
stop("m too high, not implemented")
}
o <- omega_m_cukier75[m]
return(c(o, freq_cukier(m,i+1, o)))
} else {
d_m_cukier75 <- c(4,8,6,10,20,22,32,40,38,26,56,62,46,76,96,60,86,126,134,112,92,128,154,196,34,416,106,208,328,198,382,88,348,186,140,170,284,568,302,438,410,248,448,388,596,217,100,488,166)
o <- omega_before + d_m_cukier75[m+1-i]
if (i==m){
return(o)
} else {
return(c(o, freq_cukier(m,i+1, o)))
}
}
}

