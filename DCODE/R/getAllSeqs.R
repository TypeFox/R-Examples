getAllSeqs <-
function(pep,frame,code){
  i = 1	
  mot = code[[s2c(pep)[1]]]
  while(i<nchar(pep)){
    i = i+1	
    mot =kronecker(mot,code[[s2c(pep)[i]]],paste, sep="")
  } 
  if(frame==1){
    return(kronecker(s2c("acgt"),mot,paste,sep=""))
    }
  if(frame==0){  
    return(kronecker(mot,s2c("acgt"),paste,sep=""))
    }
}
