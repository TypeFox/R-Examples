z2theta <-
function(z,family){
    if (family==1){
	#cat(paste("z: ",z,"\n"))
	if (z>2.646652) out=0.99
		else if (z<(-2.646652)) out=-0.99
		else 
        out<-(exp(2*z)-1)/(exp(2*z)+1)
	if (out>0.99) out=0.99
	if (out<(-0.99)) out=-0.99
    }
    if (family==3){
        out<-max(exp(z),0.0001)

    }
    if (family==4){
        out<-max(exp(z)+0.99999,1)
    }
    if (family==5){
        out<-z
        if (abs(out)<0.0001) out<-0.0001*sign(out)
    }
    return(out)
}
