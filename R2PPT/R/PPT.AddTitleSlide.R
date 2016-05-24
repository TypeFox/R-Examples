"PPT.AddTitleSlide" <-function(ppt,title=NULL,subtitle=NULL,title.font=NULL,title.fontsize=NULL,subtitle.font=NULL,subtitle.fontsize=NULL){

if(ppt$method=="rcom"){

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")

}

#ppt$Current.Slide <- comInvoke(comGetProperty(ppt$pres,"Slides"),"Add",comGetProperty(comGetProperty(ppt$pres,'Slides'),'Count')+1,1)
#comInvoke(ppt$Current.Slide,'Select')

ppt$Current.Slide<-ppt$pres[["Slides"]]$add(as.integer(max(1,ppt$pres[["Slides"]][["Count"]]+1)),as.integer(1))
ppt$Current.Slide$Select()

#mainseg<-comGetProperty(comGetProperty(comGetProperty(comGetProperty(ppt$Current.Slide,"Shapes"),"Title"),"TextFrame"),"TextRange")
mainseg<-ppt$Current.Slide[["Shapes"]][["Title"]][["TextFrame"]][["TextRange"]]

#if(!is.null(title))          comSetProperty(mainseg,"Text",title) 
if(!is.null(title))           mainseg[["Text"]]<-as.character(title)

#if(!is.null(title.fontsize)) comSetProperty(comGetProperty(mainseg,"Font"),"Size",as.numeric(title.fontsize))
if(!is.null(title.fontsize)){tmp<-mainseg[["Font"]]; tmp[["Size"]]<-as.numeric(title.fontsize)}

#if(!is.null(title.font))     comSetProperty(comGetProperty(mainseg,"Font"),"Name",as.character(title.font))
if(!is.null(title.font))    {tmp<-mainseg[["Font"]]; tmp[["Name"]]<-as.character(title.font)  }   


#subseg<-comGetProperty(comGetProperty(comInvoke(comGetProperty(ppt$Current.Slide,"Shapes"),"Item",2),"TextFrame"),"TextRange")
subseg<-ppt$Current.Slide[["Shapes"]]$Item(2)[["TextFrame"]][["TextRange"]]

#if(!is.null(subtitle))          comSetProperty(subseg,"Text",subtitle)
if(!is.null(subtitle))           subseg[["Text"]]<-as.character(subtitle)

#if(!is.null(subtitle.fontsize)) comSetProperty(comGetProperty(subseg,"Font"),"Size",as.numeric(subtitle.fontsize))
if(!is.null(subtitle.fontsize))  {tmp<-subseg[["Font"]]; tmp[["Size"]]<-as.numeric(subtitle.fontsize)}

#if(!is.null(subtitle.font))     comSetProperty(comGetProperty(subseg,"Font"),"Name",as.character(subtitle.font))
if(!is.null(subtitle.font))      {tmp<-subseg[["Font"]];tmp[["Name"]]<-as.character(subtitle.font)}


return(invisible(ppt))
}

