CalculPWM <-
function(matriu){
    number<-matrix(,4,ncol(matriu))
	  
       for (ii in c(1:(ncol(matriu)))) {
               number[1,ii]<-sum(as.numeric(matriu[,ii]=="A" | matriu[, ii]=="a"));number[2,ii]<-sum(as.numeric(matriu[,ii]=="C" | matriu[, ii]=="c"))
               number[3,ii]<-sum(as.numeric(matriu[,ii]=="G" | matriu[, ii]=="g"));number[4,ii]<-sum(as.numeric(matriu[,ii]=="T" | matriu[, ii]=="t"))
	       
       }
      
       transnum<-t(number)
      print(transnum)
}

