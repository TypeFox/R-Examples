# The EMbC Package for R
# 
# Copyright 2013, 2014, 2015 Joan Garriga <jgarriga@ceab.csic.es>, Aitana Oltra <aoltra@ceab.csic.es>, John R.B. Palmer <johnrbpalmer@gmail.com>, Frederic Bartumeus <fbartu@ceab.csic.es>
#   
#   EMbC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or  (at your option) any later version.
# 
# EMbC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along with this program.  If not, see http://www.gnu.org/licenses.

# Class: webMap
# -------------

setClass("webMap",
	representation(
		fName="character",
		head="character",
		path="character",
		points="character",
		end="character")
	)


# initialize
setMethod("initialize","webMap",function(.Object,bCP,bwise,folder,apiKey,mapType,markerRadius){
	dir.create(folder,showWarnings=FALSE)  
	fName <- format(Sys.time(),"%Y-%m-%d_%H:%M:%S")
	if(.Platform$OS.type=='windows') fName <- gsub(':','-',fName)
	.Object@fName <- paste(getwd(),'/',folder,'/',fName,'.html',sep='')
	if (bwise) {
		sps <- bCP@pth
		coordinates(sps) <- bCP@pth[,c('lon','lat')]}
	else {
		jttPth <- jitter(bCP@pth)
		sps <- jttPth
		coordinates(sps) <- jttPth[,c('lon','lat')]}
	proj4string(sps) <- CRS("+proj=longlat +datum=WGS84")
	.Object@head <- webMapHead(sps,apiKey,mapType)
	if (bwise) {
		.Object@path <- bwWebMapPath(bCP)
		.Object@points <- bwWebMapPoints(bCP,markerRadius)}
	else {
		.Object@path <- pwWebMapPath(jttPth)
		.Object@points <- pwWebMapPoints(bCP,jttPth,markerRadius)}
	.Object@end <- webMapEnd()
	webMFile <- file(.Object@fName)
	writeLines(paste(.Object@head,.Object@path,.Object@points,.Object@end,sep=""),sep="",con=webMFile)
	close(webMFile)
	.Object
	})


# WebMap internal functions
# -------------------------

jitter <- function(pth){
	ovrlap <- which(diff(pth[,'lon'])==0&diff(pth[,'lat'])==0)
	pth[ovrlap,'lon'] <- pth[ovrlap,'lon'] + runif(length(ovrlap),0.2,0.8)*10e-5
	pth[ovrlap,'lat'] <- pth[ovrlap,'lat'] + runif(length(ovrlap),0.2,0.8)*10e-5
	return(pth)}

webMapHead <- function(sps,apiKey,mapType){
	head <- paste('<!DOCTYPE html>
	<html style="width:100%; height:100%;">
	<head>
	<meta name="viewport" content="initial-scale=1.0, user-scalable=no" />
	<meta http-equiv="content-type" content="text/html; charset=UTF-8"/>
	<title>EMbC Webmap</title>
	<script type="text/javascript" src="http://maps.googleapis.com/maps/api/js?',apiKey,'&sensor=false">
	</script>
	<script type="text/javascript">
		function load() {
			var map=new google.maps.Map(document.getElementById("map"), {
				mapTypeId:google.maps.MapTypeId.',mapType,
				'});
			google.maps.event.addListener(map,"tilesloaded",tilesLoaded);
			function tilesLoaded() {
				google.maps.event.clearListeners(map,"tilesloaded");
				};
			var pLCoordinates;
			var pLine;
			var bounds = new google.maps.LatLngBounds();
			bounds.extend(new google.maps.LatLng(',sps@bbox[2,"min"],',',sps@bbox[1,"min"],'));
			bounds.extend(new google.maps.LatLng(',sps@bbox[2,"max"],',', sps@bbox[1,"max"],'));
			map.fitBounds(bounds);
			',sep="")
	return(head)}

pwWebMapPath <- function(jttPth){
	path <- paste('pLCoordinates=[', 
			paste(apply(jttPth,1,function(loc){
			paste('new google.maps.LatLng(',loc["lat"],',',loc["lon"],')', sep="")}),collapse=","),'];
			pLine = new google.maps.Polyline({
				path: pLCoordinates,
				geodesic: true,
				strokeColor: "white",
				strokeOpacity: 1.0,
				strokeWeight: 2.0});
			pLine.setMap(map);
			',sep="")
	return(path)}

bwWebMapPath <- function(bCP){
	bursts <- unlist(coordinates(bCP@tracks),recursive=FALSE)
	path <- paste(unlist(lapply(1:length(bursts),function(i){
		paste('pLCoordinates=[',
			paste(apply(bursts[[i]],1,function(loc)
			paste('new google.maps.LatLng(',loc["lat"],',',loc["lon"],')',sep="")),collapse=","),'];
			pLine = new google.maps.Polyline({
				path: pLCoordinates,
				geodesic: true,
				strokeColor: "',bCP@C[bCP@midPoints@data[i,"A"]],'",
				strokeOpacity: 1.0,
				strokeWeight: 2.0});
			pLine.setMap(map);
			',sep="")}
		)),collapse="")
	return(path)}

pwWebMapPoints <- function(bCP,jttPth,markerRadius){
	locLbl <- c('LL','LH','HL','HH')
	points <- paste('
		var marker;
		var infoWindow = new google.maps.InfoWindow;
		var html;
		',paste(sapply(1:nrow(bCP@pth), function(i) {
			paste('
			marker = new google.maps.Marker({
				position:new google.maps.LatLng(',jttPth[i,"lat"],',',jttPth[i,"lon"],'),
				icon:{
					path:google.maps.SymbolPath.CIRCLE,
					scale:',markerRadius,',
					fillColor:"',bCP@C[bCP@A[i]],'",
					fillOpacity:',strtoi("0xff")/255,',
					strokeColor:"black",
					strokeWeight:1.0},
				draggable:false,
				map:map});
				html = \'<TABLE style="font-size:8pt;background-color:#cccccc;border-style:solid;border-color:white;border-collapse:collapse;" border="1" cellspacing="2" cellpadding="2">',
				'<TR align="center"><TD colspan=2>',locLbl[bCP@A[i]],'</TD></TR>',
				'<TR align="center"><TD>recIdx.</TD><TD>',i,'</TD></TR>',
				'<TR align="center"><TD>date</TD><TD>',substr(bCP@pth[i,"dTm"],1,10),'</TD></TR>',
				'<TR align="center"><TD>time</TD><TD>',substr(bCP@pth[i,"dTm"],12,19),'</TD></TR>',
				'<TR align="center"><TD>speed(m/s)</TD><TD>',round(bCP@X[i,1],4),'</TD></TR>',
				'<TR align="center"><TD>turn(r)</TD><TD>',round(bCP@X[i,2],4),'</TD></TR>',
				'<TR align="center"><TD>dltT.(s)</TD><TD>',round(bCP@spn[i],1),'</TD></TR>',
				'<TR align="center"><TD>dltS.(m)</TD><TD>',round(bCP@dst[i],1),'</TD></TR>',
				'</TABLE>\';
				bindInfoWindow(marker,map,infoWindow,html);
				',sep="")}),
			collapse=""),
			'function bindInfoWindow(marker,map,infoWindow,html){
				google.maps.event.addListener(marker,"click",function(){
					infoWindow.setContent(html);
					infoWindow.open(map,marker);
					});
				}
			',sep="")
	return(points)}

bwWebMapPoints <- function(bCP,markerRadius){
	locLbl <- c('LL','LH','HL','HH')
	markerSizes <- setMarkerSizes(bCP,nMarkerSizeClasses=3,minMarkerRadius=5,maxMarkerRadius=10)
	points <- paste('
		var marker;
		var infoWindow=new google.maps.InfoWindow;
		var html;
		',paste(sapply(1:length(bCP@midPoints), function(i) {
			paste('
			marker = new google.maps.Marker({
				position:new google.maps.LatLng(',coordinates(bCP@midPoints)[i,"lat"],',',coordinates(bCP@midPoints)[i,"lon"],'),
				icon:{
					path:google.maps.SymbolPath.CIRCLE,
					scale:',markerSizes[i],',
					fillColor:"',bCP@C[bCP@midPoints@data[i,"A"]],'",
					fillOpacity:',strtoi("0xff")/255,',
					strokeColor:"black",
					strokeWeight:1.0},
				draggable:false,
				map:map});
				html = \'<TABLE style="font-size:8pt;background-color:#cccccc;border-style:solid;border-color:white;border-collapse:collapse;" border="1" cellspacing="2" cellpadding="2">',
				'<TR align="center"><TD colspan=2>',locLbl[bCP@midPoints@data[i,"A"]],'</TD></TR>',
				'<TR align="center"><TD>recIdx.</TD><TD>',paste(bCP@midPoints@data[i,"first"], bCP@midPoints@data[i,"last"], sep="-"),'</TD></TR>',
				'<TR align="center"><TD>duration</TD><TD>',formatSecs(bCP@midPoints@data[i,"duration"]),'</TD></TR>',
				'<TR align="center"><TD>distance</TD><TD>',formatMeters(bCP@midPoints@data[i,"distance"]),'</TD></TR>',
				'</TABLE>\';
				bindInfoWindow(marker,map,infoWindow,html);
				',sep="")}),
			collapse=""),
			'function bindInfoWindow(marker,map,infoWindow,html) {
				google.maps.event.addListener(marker,"click",function() {
					infoWindow.setContent(html);
					infoWindow.open(map, marker);
					});
				}
			', sep="")
	return(points)}

webMapEnd <- function(){
	end <- '}
	</script>
	</head>
	<body onload="load()"
		style="width:100%;
		height:100%;
		position:relative;
		margin:0px">
	<div id="map" style="width:100%;height:100%;">
	</div>
	</body>
	</html>
	'
	return(end)}

