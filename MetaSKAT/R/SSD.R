

Print_Error_SSD<-function(code){

	if(code == 0 ){
		return(0);
	} else if(code == 1){
		stop("Error Can't open BIM file")
	} else if(code == 2){
		stop("Error Can't open FAM file")
	} else if(code == 3){
		stop("Error Can't open BED file")
	} else if(code == 4){
		stop("Error Can't open SETID file")
	} else if(code == 5){
		stop("Error Can't write SSD file")
	} else if(code == 6){
		stop("Error Can't read SSD file")
	} else if(code == 7){
		stop("Error Can't write INFO file")
	} else if(code == 8){
		stop("Error Can't read INFO file")
	} else if(code == 9){
		stop("Error Can't write INFO file")
	} else if(code == 13){
		stop("Error Wrong SNP or Individual sizes")
	} else if(code == 14){
		stop("Error SetID not found")
	} else {
		MSG<-sprintf("Error [%d]\n",code)
		stop(MSG)
	}
	
	return(1)
}


Check_File_Exists<-function(FileName){
	
	if(!file.exists(FileName)){
		Msg<-sprintf("File %s does not exist\n",FileName)
		stop(Msg)
	}

}

MetaSKAT_Is_IsLittleEndian<-function(){

	re<-0
	temp<-.C("IsLittleEndian", as.integer(re))
	re1<-temp[[1]]
	
	if(re1 == 0){
		stop("MSSD files can only be used on a little endian machine")
	}
}



