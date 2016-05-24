"PPT.SaveAs" <-function(ppt,file=stop("filename must be specified.")){

if (ppt$method == "rcom") {  

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")

}

file<-PPT.getAbsolutePath(file)
file<-gsub("/","\\\\",as.character(file)) #character should be \\\\ and not / some compatability issues in rcom. 

ppt$pres$SaveAs(file)
#comInvoke(ppt$pres,"SaveAs",gsub("/","\\\\",as.character(file)))

return(invisible(ppt))
}

