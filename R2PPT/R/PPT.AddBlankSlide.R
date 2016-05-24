"PPT.AddBlankSlide" <- function(ppt){

if(ppt$method=="rcom"){

	if(!comIsValidHandle(ppt$ppt))   stop("Invalid handle for powerpoint application")
	if(!comIsValidHandle(ppt$pres))  stop("Invalid handle for powerpoint presentation")

}


ppt$Current.Slide<-ppt$pres[["Slides"]]$add(as.integer(max(1,ppt$pres[["Slides"]][["Count"]]+1)),as.integer(12))
ppt$Current.Slide$Select()
#ppt$Current.Slide <- comInvoke(comGetProperty(ppt$pres,"Slides"),"Add",comGetProperty(comGetProperty(ppt$pres,'Slides'),'Count')+1,12)
#comInvoke(ppt$Current.Slide,'Select')

return(invisible(ppt))

}

