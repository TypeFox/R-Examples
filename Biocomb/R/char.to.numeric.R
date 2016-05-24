 do.numeric<-function(char){
 
 num<-NULL
 
    if(grepl("^\\d+\\,\\d+$", char)){
	    num<-as.numeric(sub(",",".",char))
	 }
	 
	 if(grepl("^\\d+\\/\\d+$", char)){
	    str<-unlist(strsplit(char,"[ /]"))
        num<-as.numeric(str[1])/as.numeric(str[2])
	 }
	 
     if(grepl("^\\d+$", char)|grepl("^\\d+\\.\\d+$",char)){
	   num<-as.numeric(char)
	 }
	 
return(num)	 
}

char.to.numeric<-function(char){

   if(substring(char,1,1)=="-"){
     char<-substring(char,2)
     num<-do.numeric(char)
	 num<-num*(-1)
   }else{
      num<-do.numeric(char)
   }
   
   return(num)
}