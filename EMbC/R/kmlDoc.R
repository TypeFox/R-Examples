# The EMbC Package for R
# 
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#   
#   EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
# 
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.


# Class: kmlDoc
# -------------

setClass("kmlDoc",
	representation(
		bCP="binClstPath",
		fName="character",
		head="character",
		pthH="character",
		pthB="character",
		pthM="character",
		end="character",
		type="character")
	)

# initialize
setMethod("initialize","kmlDoc",function(.Object,bCP,bwise,folder,markerRadius){
	dir.create(folder,showWarnings=FALSE)  
	fName <- format(Sys.time(),"%Y-%m-%d_%H:%M:%S")
	if(.Platform$OS.type=='windows') fName <- gsub(':','-',fName)
	.Object@fName <- paste(getwd(),'/',folder,'/',fName,'.kml',sep='')
	.Object@head <- kmlHead()         					# document head
	.Object@pthH <- kmlPthH()         					# path head
	if(bwise)											# pah body:
		.Object@pthB <- kmlBurstTracks(bCP)				# bursted
	else
		.Object@pthB <- kmlPthB(bCP)     				# pointwise
	.Object@pthM <- kmlPthM(bCP,bwise,markerRadius)		# path markers
	.Object@end  <- kmlEnd()        					# document end
	kmlFile <- file(.Object@fName)
	writeLines(c(.Object@head,.Object@pthH,.Object@pthB,.Object@pthM,.Object@end),sep='',con=kmlFile)
	close(kmlFile)
	.Object
	})


# kmlDoc internal functions
# -------------------------

bgrC <- function(hexrgb){
	return(paste('7f',substr(hexrgb,6,7),substr(hexrgb,4,5),substr(hexrgb,2,3),sep=''))}

kmlHead <- function() {
	txt <- character(0)
	txt <- c(txt, '<?xml version="1.0" encoding="UTF-8"?>\n')
	txt <- c(txt, '<kml xmlns="http://www.opengis.net/kml/2.2">\n')
	txt <- c(txt, 'xmlns:gx="http://www.google.com/kml/ext/2.2"\n')
	txt <- c(txt, 'xmlns:kml="http://www.opengis.net/kml/2.2"\n')
	txt <- c(txt, 'xmlns:atom="http://www.w3.org/2005/Atom">\n')
	txt <- c(txt, '<Document>\n')
	return (txt)}

kmlPthH <- function() {
	txt <- character(0)
	txt <- c(txt, '<name>Paths</name>\n')
	txt <- c(txt, '<description>kmlPath</description>\n')
	txt <- c(txt, '<Style id="kmlStyles">\n')
	txt <- c(txt, '<LineStyle>\n')
	txt <- c(txt, '<color>7f00ffff</color>\n')
	txt <- c(txt, '<width>2</width>\n')
	txt <- c(txt, '</LineStyle>\n')
	txt <- c(txt, '<PolyStyle>\n')
	txt <- c(txt, '<color>7f00ff00</color>\n')
	txt <- c(txt, '</PolyStyle>\n')
	txt <- c(txt, '<IconStyle>\n')
	txt <- c(txt, '<Icon>\n')
	txt <- c(txt, '<href>http://maps.google.com/mapfiles/kml/pal2/icon26.png</href>\n')
	txt <- c(txt, '</Icon>\n')
	txt <- c(txt, '<scale>0.4</scale>\n')
	txt <- c(txt, '</IconStyle>\n')
	txt <- c(txt, '<BalloonStyle>\n')
	txt <- c(txt, '<text>$[description]</text>\n')
	txt <- c(txt, '<bgColor>7fcccccc</bgColor>\n')
	txt <- c(txt, '</BalloonStyle>\n')
	txt <- c(txt, '</Style>\n')
	return (txt)}

kmlBurstTracks <- function(bCP){
	return(unlist(lapply(1:length(bCP@tracks),
		 function(i){ 
			kmlBurstPthB(coordinates(bCP@tracks)[[i]][[1]],bgrC(bCP@C[bCP@tracks@data$A[i]]), i)})))}

kmlBurstPthB <- function(bCPMat,burstColor,drawOrder) {
	txt <- character(0)
	txt <- c(txt, '<Placemark>\n')
	txt <- c(txt, '<LineString>\n')
	txt <- c(txt, '<gx:drawOrder>', drawOrder,'</gx:drawOrder>\n')
	txt <- c(txt, '<extrude>1</extrude>\n')
	txt <- c(txt, '<tessellate>1</tessellate>\n')
	txt <- c(txt, '<altitudeMode>clampToGround</altitudeMode>\n')
	txt <- c(txt, '<coordinates>\n')
	txt <- c(txt, apply(bCPMat, 1, function(x) paste(x["lon"], ",", x["lat"],"\n", sep="")))
	txt <- c(txt, '</coordinates>\n')
	txt <- c(txt, '</LineString>\n')
	txt <- c(txt, '<Style>\n')
	txt <- c(txt, '<LineStyle>\n')
	txt <- c(txt, '<color>', burstColor, '</color>\n')
	txt <- c(txt, '<gx:outerColor>64000000</gx:outerColor>\n')
	txt <- c(txt, '<width>2</width>\n')
	txt <- c(txt, '<gx:outerWidth>.2</gx:outerWidth>\n')
	txt <- c(txt, '</LineStyle>\n')
	txt <- c(txt, '</Style>\n')
	txt <- c(txt, '</Placemark>\n')
	return (txt)}

kmlPthB <- function(bCP) {
	txt <- character(0)
	txt <- c(txt, '<Placemark>\n')
	txt <- c(txt, '<LineString>\n')
	txt <- c(txt, '<extrude>0</extrude>\n')
	txt <- c(txt, '<tessellate>1</tessellate>\n')
	txt <- c(txt, '<altitudeMode>clampToGround</altitudeMode>\n')
	txt <- c(txt, '<coordinates>\n')
	txt <- c(txt,apply(bCP@pth,1,function(x) paste(x["lon"],",",x["lat"],"\n",sep="")))
	txt <- c(txt, '</coordinates>\n')
	txt <- c(txt, '</LineString>\n')
	txt <- c(txt, '</Placemark>\n')
	return (txt)}

kmlEnd <- function(){
	txt <- character(0)
	txt <- c(txt, '</Document>\n')
	txt <- c(txt, '</kml>\n')
	return (txt)}

mrkHead <- function(bCP,bwise,i,loc,lbl,bgr,markerSize=NA){
	txt <- character(0)
	txt <- c(txt, '<Placemark id="point">\n')
	txt <- c(txt, '<visibility>1</visibility>\n')
	txt <- c(txt, '<styleUrl>#kmlStyles</styleUrl>\n')
	txt <- c(txt, '<description>\n')
	txt <- c(txt, '<![CDATA[<TABLE style="font-size:8pt;background-color:#cccccc;border-style:solid;border-color:white;border-collapse:collapse;" border="1" cellspacing="2" cellpadding="2">')
	txt <- c(txt, paste('<TR align="center"><TD colspan=2>',lbl,'</TD></TR>',sep=''))
	if(bwise){
		txt <- c(txt, paste('<TR align="center"><TD>recIdx.	 </TD><TD>',paste(bCP@midPoints@data[i, "first"],bCP@midPoints@data[i, "last"],sep="-"),'</TD></TR>',sep=''))
		txt <- c(txt, paste('<TR align="center"><TD>duration	 </TD><TD>',formatSecs(bCP@midPoints@data[i, "duration"]),'</TD></TR>',sep=''))
		txt <- c(txt, paste('<TR align="center"><TD>distance	 </TD><TD>',formatMeters(bCP@midPoints@data[i, "distance"]),'</TD></TR>',sep=''))
		}		
	else {
		txt <- c(txt, paste('<TR align="center"><TD>recIdx.	 </TD><TD>',i,'</TD></TR>',sep=''))
		txt <- c(txt, paste('<TR align="center"><TD>date	 </TD><TD>',substr(loc$dTm,1,10),'</TD></TR>',sep=''))
		if (nchar(substr(loc$dTm,12,19))>0)
			txt <- c(txt, paste('<TR align="center"><TD>time	 </TD><TD>',substr(loc$dTm,12,19),'</TD></TR>',sep=''))
		else
			txt <- c(txt, paste('<TR align="center"><TD>time	 </TD><TD>',"00:00:00",'</TD></TR>',sep=''))
		for (m in seq(bCP@m))
			txt <- c(txt, paste('<TR align="center"><TD>X',m,'</TD><TD>',round(bCP@X[,m][i],4),'</TD></TR>',sep=''))
		txt <- c(txt, paste('<TR align="center"><TD>dltT.(s)	</TD><TD>',round(bCP@spn[i],1),'</TD></TR>',sep=''))
		txt <- c(txt, paste('<TR align="center"><TD>dltS.(m)	</TD><TD>',round(bCP@dst[i],1),'</TD></TR>',sep=''))
		}
	txt <- c(txt, '</TABLE>]]>\n')
	txt <- c(txt, '</description>\n')
	txt <- c(txt, '<Style>\n')
	txt <- c(txt, '<IconStyle>\n')
	txt <- c(txt, paste('<color>',bgr,'</color>\n',sep=''))
	txt <- c(txt, paste('<scale>',markerSize,'</scale>\n',sep=''))
	txt <- c(txt, '</IconStyle>\n')
	txt <- c(txt, '</Style>\n')
	return (txt)}

mrkBody <- function(loc){
	txt <- character(0)
	txt <- c(txt, '<Point id="point_point">\n')
	txt <- c(txt, '<altitudeMode>relativeToGround</altitudeMode>\n')
	txt <- c(txt, '<tessellate>1</tessellate>\n')
	txt <- c(txt, '<extrude>1</extrude>\n')
	txt <- c(txt, '<coordinates>\n')
	txt <- c(txt, paste(loc$lon,loc$lat,'0\n',sep=','))
	txt <- c(txt, '</coordinates>\n')
	txt <- c(txt, '</Point>\n')
	return (txt)}

mrkTail <- function(){
	return ('</Placemark>\n')}

kmlPthM <- function(bCP,bwise,markerRadius){
	locLbl <- getkLbls(bCP,kNmbrs=TRUE)
	txt <- character(0)	
	if(bwise){
		# treat the result as a pixel diameter and scale it for a 16-pixel radius icon:
		markerSizes <- setMarkerSizes(bCP,nMarkerSizeClasses=4,minMarkerRadius=5,maxMarkerRadius=markerRadius)/16
		for (i in (1:nrow(bCP@midPoints))){
			loc <- bCP@midPoints[i,]
			txt <- c(txt,mrkHead(bCP,bwise,i,loc,locLbl[bCP@midPoints@data$A[i]],bgrC(bCP@C[bCP@midPoints@data$A[i]]),markerSize=markerSizes[i]))
			txt <- c(txt,mrkBody(loc))
			txt <- c(txt,mrkTail())    
			}		
		}
	 else { 
		# markerRadius is in pixels and we need to scale it for a 16-pixel radius icon so:
		for (i in 1:bCP@n){
			loc <- bCP@pth[i,]
			txt <- c(txt,mrkHead(bCP,bwise,i,loc,locLbl[bCP@A[i]],bgrC(bCP@C[bCP@A[i]]),markerSize=markerRadius/16))
			txt <- c(txt,mrkBody(loc))
			txt <- c(txt,mrkTail())    
			}
		}
	return (txt)}

