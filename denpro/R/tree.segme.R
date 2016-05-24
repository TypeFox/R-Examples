tree.segme<-function(tt,paletti=NULL,pcf=NULL)
{

if (is.null(paletti))
 paletti<-c("red","blue","green",
 "orange","navy","darkgreen",
 "orchid","aquamarine","turquoise",
 "pink","violet","magenta","chocolate","cyan",
 colors()[50:657],colors()[50:657])

colors<-colobary(tt$parent,paletti)
if (is.null(pcf)) segme<-colors
else{
  lenni<-length(pcf$value)
  segme<-matrix(0,lenni,1)
}
for (i in 1:length(colors)) segme[tt$infopointer[i]]<-colors[i]

return(segme)
}






