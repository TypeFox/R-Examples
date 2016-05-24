Overlay2GE = function(coords, image = "/home/Images/myimage.jpg", goo = "Overlay2GE.kml"){
cat("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n
<kml xmlns=\"http://earth.google.com/kml/2.1\">
<GroundOverlay>
	<name>Image</name>
	<color>b2ffffff</color>
	<Icon>
	  <href>",image,"</href>
	</Icon>
	<LatLonBox>
		<north>",coords[1],"</north>
		<south>",coords[2],"</south>
		<east>",coords[3],"</east>
		<west>",coords[4],"</west>
	</LatLonBox>
</GroundOverlay>
</kml>", file = goo)
}

