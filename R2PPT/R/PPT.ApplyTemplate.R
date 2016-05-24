"PPT.ApplyTemplate" <-function(ppt,file){

if(ppt$method=="rcom"){

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")

}

file<-PPT.getAbsolutePath(file[1]) #New in Version 1.1
if(!file.exists(file)) stop(paste(file, "does not exist"))
file<-gsub("/","\\\\",as.character(file))

#comInvoke(ppt$pres,"ApplyTemplate",file)
ppt$pres$ApplyTemplate(file)

return(invisible(ppt))
}

