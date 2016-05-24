`freq_mcrae82` <-
function(m, i=1, omega_before=0){
if(i<=1){
if(m>NROW(omega_m_mcrae82)){
print("m too high, not implemented")
return
}
o <- omega_m_mcrae82[m]
return(c(o, freq_mcrae82(m,i+1, o)))
} else {
o <- omega_before + d_m_mcrae82[m+1-i]
if (i==m){
return(o)
} else {
return(c(o, freq_mcrae82(m,i+1, o)))
}
}
}

