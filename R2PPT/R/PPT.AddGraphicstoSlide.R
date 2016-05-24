"PPT.AddGraphicstoSlide" <-function(ppt,file=NULL,size=c(10,10,700,500),dev.out.type = "jpeg"){

#Check a slide exists

if(ppt$method=="rcom"){

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")
	if(!comIsValidHandle(ppt$Current.Slide))  stop("Invalid handle for presentation slide. Make sure you add a slide before adding graphic.")

}

if(length(size)!=4) stop("Graphic size to export to PowerPoint must be a vector of length 4")


if(!is.null(file[1])){file[1]<-PPT.getAbsolutePath(file[1])} #New in Version 1.1

if(!is.null(file[1]) && !file.exists(file[1])) stop(paste(file[1],"does not exist"))

if(!is.null(file[1]) && file.exists(file[1])){

file<-gsub("/","\\\\",as.character(file[1]))

#myShapes<-comGetProperty(ppt$Current.Slide,'Shapes')
myShapes<-ppt$Current.Slide[["Shapes"]]

#myPicture<-comInvoke(myShapes,'AddPicture',file[1],0,-1,size[1],size[2],size[3],size[4])  # xmargin from left (720), y margin from top.
myShapes$AddPicture(file[1],0,-1,size[1],size[2],size[3],size[4]) 

return(invisible(ppt))

}


if(is.na(match(dev.cur(), dev.list(), NA))){stop("No active graphic device to export to PowerPoint.")}

PPTtemp<-tempfile()
PPTtemp<-paste(tempfile(),dev.out.type,sep=".")
savePlot(PPTtemp,type=dev.out.type)
if(file.exists(paste(PPTtemp,dev.out.type,sep="."))){file.rename(paste(PPTtemp,dev.out.type,sep="."),PPTtemp)}
## Different versions of function savePlot sometimes appends the file extension and sometimes not hence above!

#myShapes<-comGetProperty(ppt$Current.Slide,'Shapes')
myShapes<-ppt$Current.Slide[["Shapes"]]

#myPicture<-comInvoke(myShapes,'AddPicture',PPTtemp,0,-1,size[1],size[2],size[3],size[4]) 
myShapes$AddPicture(PPTtemp,0,-1,size[1],size[2],size[3],size[4]) 
 
unlink(PPTtemp)
return(invisible(ppt))    
}

