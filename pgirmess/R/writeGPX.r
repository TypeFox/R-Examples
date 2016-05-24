writeGPX<-function(x,filename="",type="w") {
 if(toupper(substr(filename,nchar(filename)-3,nchar(filename)))!=".GPX" & filename!="") filename<-paste(filename,".gpx",sep="")
   options(digits=7)
   cat('<?xml version="1.0"?>\n<gpx version="1.1" creator="pgirbric">\n',file=filename,sep="") # header
if (type=="w") {
   for(i in 1:nrow(x)){
         cat('<wpt lat="',x[i,3],'" lon="',x[i,2],'">\n',file=filename,sep="",append=TRUE)
         cat('<name>',as.character(x[i,1]),'</name>\n',file=filename,sep="",append=TRUE)
         if(ncol(x)>3) cat('<ele>',x[i,4],'</ele>\n',file=filename,sep="",append=TRUE)
         cat('<sym>Flag</sym>\n',file=filename,sep="",append=TRUE)
         cat('</wpt>\n',file=filename,sep="",append=TRUE)
   }
   }
else {
cat('<trk>\n',file=filename,sep="",append=TRUE)
cat('<name>',filename,'</name>\n',file=filename,sep="",append=TRUE)
cat('<trkseg>\n',file=filename,sep="",append=TRUE)
     for(i in 1:nrow(x)){
         cat('<trkpt lat="',x[i,3],'" lon="',x[i,2],'">\n',file=filename,sep="",append=TRUE)
         if(ncol(x)>2) cat('<ele>',x[i,4],'</ele>\n',file=filename,sep="",append=TRUE)
         cat('</trkpt>\n',file=filename,sep="",append=TRUE)
         }
cat('</trkseg>\n</trk>\n',file=filename,sep="",append=TRUE)
}
   cat('</gpx>\n',file=filename,sep="",append=TRUE) # bottom
   options(digits=3)
}


