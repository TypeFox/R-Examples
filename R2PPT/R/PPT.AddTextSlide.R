"PPT.AddTextSlide" <-function(ppt,title=NULL,title.fontsize=NULL,title.font=NULL,text=NULL,text.fontsize=NULL,text.font=NULL){

if(ppt$method=="rcom"){

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")

}


#ppt$Current.Slide <- comInvoke(comGetProperty(ppt$pres,"Slides"),"Add",comGetProperty(comGetProperty(ppt$pres,'Slides'),'Count')+1,as.integer(2))
ppt$Current.Slide<-ppt$pres[["Slides"]]$add(as.integer(max(1,ppt$pres[["Slides"]][["Count"]]+1)),as.integer(2))

ppt$Current.Slide$Select()
#ppt<<-ppt
#stop()

#mainseg<-comGetProperty(comGetProperty(comGetProperty(comGetProperty(ppt$Current.Slide,"Shapes"),"Title"),"TextFrame"),"TextRange")
mainseg<-ppt$Current.Slide[["Shapes"]][["Title"]][["TextFrame"]][["TextRange"]]

#if(!is.null(title))          comSetProperty(mainseg,"Text",title) 
if(!is.null(title))           mainseg[["Text"]]<-as.character(title)


#if(!is.null(title.fontsize)) comSetProperty(comGetProperty(mainseg,"Font"),"Size",as.numeric(title.fontsize))
if(!is.null(title.fontsize)){tmp<-mainseg[["Font"]]; tmp[["Size"]]<-as.numeric(title.fontsize)}

#if(!is.null(title.font))     comSetProperty(comGetProperty(mainseg,"Font"),"Name",as.character(title.font))
if(!is.null(title.font))    {tmp<-mainseg[["Font"]]; tmp[["Name"]]<-as.character(title.font)  }   

#textseg<-comGetProperty(comGetProperty(comInvoke(comGetProperty(ppt$Current.Slide,"Shapes"),"Item",2),"TextFrame"),"TextRange")
textseg<-ppt$Current.Slide[["Shapes"]]$Item(2)[["TextFrame"]][["TextRange"]]


#if(!is.null(text))          comSetProperty(textseg,"Text",text) 
if(!is.null(text))           textseg[["Text"]]<-as.character(text)

#if(!is.null(text.fontsize)) comSetProperty(comGetProperty(textseg,"Font"),"Size",as.numeric(text.fontsize))
if(!is.null(text.fontsize))  {tmp<-textseg[["Font"]]; tmp[["Size"]]<-as.numeric(text.fontsize)}

#if(!is.null(text.font))     comSetProperty(comGetProperty(textseg,"Font"),"Name",as.character(text.font))
if(!is.null(text.font))      {tmp<-textseg[["Font"]];tmp[["Name"]]<-as.character(text.font)}

return(invisible(ppt))

}

