###############################################
# NAME: 
# PURPOSE:
# INPUTS:
# OUTPUTS:
# CHANGES: 20/04/2010 - bug en el computo de Lat
# CHANGES: 20/04/2010 - cambio en el computo de dia
# CHANGES: 20/04/2010 - quitar valor por defecto de day
# CHANGES: 23/05/2010 - stop si no es latlong
# CHANGES: 07/02/2013 - adecuar a sp > 0.9-90 (slot coord no existe)
# TODO: si no es latlong proyectar a latlong calcular rad y reproyectar a original
# TODO: Justificar el 898
# TODO: ampliar para calculo simplificado de la FAO
###############################################
solarRad = function (img, day) {	
	
	if (is.projected(img)) stop('Cant calculate extraterrestial radiation over projected images')
	
	DTOR=0.0174533 #cte de grados a radianes
	
	nRow=img@grid@cells.dim[2]
	nCol=img@grid@cells.dim[1]
	rowSize=img@grid@cellsize[2]
	# CHANGE lowerY=a@coords[3] no funciona a partir de version de sp 0.9-91
	# uso coordinates
	lowerY=min(coordinates(img)[,2]) #Y del punto central
	lat=as.matrix(img)
	
	lat[1:nCol,]=rep((0:(nRow-1)*rowSize)+lowerY,each=nCol) #imagen donde cada pixel tiene el valor central de su coord Y 
	
	lat=lat*DTOR
	
	dia=2*pi/365*(day-1)
	DST=1.00011+0.034221*cos(dia)+0.00128*sin(dia)+0.000719*cos(2*dia)+0.000777*sin(2*dia)
	DEC=0.006918-0.399912*cos(dia)+0.070257*sin(dia)-0.006758*cos(2*dia)+0.000907*sin(2*dia)-0.002697*cos(3*dia)+0.00148*sin(3*dia)
	AGH=acos(-tan(lat)*tan(DEC))
	RAD=898*DST*(sin(lat)*sin(DEC)*AGH+cos(lat)*cos(DEC)*sin(AGH))
	
	img$band1=rev(as.vector(RAD))
	
	img	
}

###############################################
# NAME: 
# PURPOSE:
# INPUTS:
# OUTPUTS:
###############################################
solarRad12M = function (img, outFl, ...) {
	
	if (length(outFl)!=12) stop('Fl must be 12 file names')
	
	# valor juliano del dia 15 de cada mes
	#TODO: comprobarlo
	DDA=c(15,45,75,106,136,167,197,228,259,289,320,350)
	
	for (i in 1:12) writeGDAL(solarRad(img,DDA[i]),outFl[i],...)
	outFl
}
