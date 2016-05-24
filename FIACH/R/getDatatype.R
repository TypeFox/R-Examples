getDatatype<-function(input,type="RNiftyReg"){
  if(mode(input)=="list"||mode(input)=="character"){
    outype<-max(unlist(lapply(input,function(x){RNiftyReg::dumpNifti(x)$datatype})))
  }else{
    outype<-RNiftyReg::dumpNifti(image = input)$datatype
    
  }
  if(type=="RNiftyReg"){
    outcode<-as.character(datatypes[datatypes[,2]==outype,1])
  }else{
    outcode<-datatypes[datatypes[,2]==outype,2]  
  }
  if(length(outcode)!=1){outcode<-"auto"}
  return(outcode)
}
