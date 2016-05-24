"PPT.Close" <-function(ppt){

if(ppt$method=="rcom"){

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")

}

#comInvoke(ppt$pres,"Close")
ppt$pres$Close()


return(invisible(ppt))
}

