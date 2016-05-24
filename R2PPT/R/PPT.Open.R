"PPT.Open" <- function(file=stop("filename must be specified"),method=c("rcom","RDCOMClient")){

file<-PPT.getAbsolutePath(file[1]) #New in Version 1.1
if(!file.exists(file[1])) stop(paste(file[1],"does not exist")) 

ppt<-PPT.Init(visible=TRUE,method=method,addPres=FALSE)

ppt$pres <-ppt$ppt[["Presentations"]]$Open(file)

#ppt$pres <- comInvoke(comGetProperty(ppt$ppt, "Presentations"),"Open",file)

#    if(!require(rcom))  stop("library rcom unavailable")
    
#    file<-PPT.getAbsolutePath(file[1]) #New in Version 1.1
#    if(!file.exists(file[1])) stop(paste(file[1],"does not exist")) 
    

#    file <- gsub("/", "\\\\", as.character(file[1]))

#    ppt <- list()
#    ppt$ppt <- comCreateObject("PowerPoint.Application")
#    if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
#    comSetProperty(ppt$ppt, "Visible", TRUE)

#    ppt$pres <- comInvoke(comGetProperty(ppt$ppt, "Presentations"),"Open",file)
#    if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")


    return(invisible(ppt))
}

