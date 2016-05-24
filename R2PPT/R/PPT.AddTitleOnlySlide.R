"PPT.AddTitleOnlySlide" <-function(ppt,title=NULL,title.fontsize=NULL,title.font=NULL){

if(ppt$method=="rcom"){

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")

}

#ppt$Current.Slide <- comInvoke(comGetProperty(ppt$pres,"Slides"),"Add",comGetProperty(comGetProperty(ppt$pres,'Slides'),'Count')+1,11)
ppt$Current.Slide<-ppt$pres[["Slides"]]$add(as.integer(max(1,ppt$pres[["Slides"]][["Count"]]+1)),as.integer(11))

#comInvoke(ppt$Current.Slide,'Select')
ppt$Current.Slide$Select()

#mainseg<-comGetProperty(comGetProperty(comGetProperty(comGetProperty(ppt$Current.Slide,"Shapes"),"Title"),"TextFrame"),"TextRange")
mainseg<-ppt$Current.Slide[["Shapes"]][["Title"]][["TextFrame"]][["TextRange"]]

#if(!is.null(title))          comSetProperty(mainseg,"Text",title) 
if(!is.null(title))           mainseg[["Text"]]<-as.character(title)

#if(!is.null(title.fontsize)) comSetProperty(comGetProperty(mainseg,"Font"),"Size",as.numeric(title.fontsize))
if(!is.null(title.fontsize)){tmp<-mainseg[["Font"]]; tmp[["Size"]]<-as.numeric(title.fontsize)}

#if(!is.null(title.font))     comSetProperty(comGetProperty(mainseg,"Font"),"Name",as.character(title.font))
if(!is.null(title.font))    {tmp<-mainseg[["Font"]]; tmp[["Name"]]<-as.character(title.font)  }   



return(invisible(ppt))

}

