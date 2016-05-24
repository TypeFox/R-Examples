# PACKAGE fisheyeR
# Aux functions


puntosMedios = function(Pcoords, detalle = 5){
  for (i in 1:detalle){
    new_pcoords = matrix(rep(0,4*nrow(Pcoords)), nrow = 2* nrow(Pcoords), byrow = T )
    cont = 0
    for (i in 1:nrow(Pcoords)){
	   if (i == nrow(Pcoords)) {
  		cont = cont + 1
  		new_pcoords[cont,] = Pcoords[i,]
  		cont = cont + 1
  		new_pcoords[cont,] = Pcoords[i,] - ((Pcoords[i,]-Pcoords[1,])/2)
  	}else{
  		cont = cont + 1
  		new_pcoords[cont,] = Pcoords[i,]
  		cont = cont + 1
  		new_pcoords[cont,] = Pcoords[i,] - ((Pcoords[i,]-Pcoords[i+1,])/2)}}
    Pcoords = new_pcoords}
    return(Pcoords)
   }

fishIout = function(x, value){
  d = value
 	if (x > 0){
		signo = 1
	}else{
		signo = -1
	}
	x = abs(x)
	return(signo*(-(x/((d*x)-d-1))))
  }

fishIin = function(x, value){
  d = value
	if (x > 0){
		signo = 1
	}else{
		signo = -1
	}
	x = abs(x)

	return(signo*(((d+1)*x)/(d*x+1)))
  }

toPolar = function(x, y){
	t1 = atan2(y,x)
	rP = sqrt(x^2+y^2)
	return(c(t1 = t1,rP = rP))
  }

toCartesian = function(t1, rP){
	x1 = rP*cos(t1)
	y1 = rP*sin(t1)
	return(c(x = x1,y = y1))
  }

circulo =  function(cx, cy, r, circleCol, PLOT = TRUE){
	t = seq(0,2*pi,length=100)
	circle = t(rbind(cx+sin(t)*r,cy+cos(t)*r))
	if (PLOT == TRUE) plot(circle,type='l',,ylim=c(-1.15,1.15),xlim=c(-1.15,1.15),
		ann=FALSE, axes=F, col = circleCol)
	return(circle)
  }

circulin = function(cx, cy, r = 0.045, objeto, col = 'blue', PLOT = TRUE, label = 0){
	t = seq(0,2*pi,length=100)
	circle = t(rbind(cx+sin(t)*r,cy+cos(t)*r))
	points(circle,type='l', col = col)
	if (label != 0) text(cx,cy,label,cex = .7)
	insiders <- apply(objeto,1,function(co)(cx-co[1])^2+(cy-co[2])^2<r^2)
  assign('insiders', insiders , envir = POI.env)
  }

addNoise = function(m, tamanyo = 0.01){
	noise = function(m, t = tamanyo){
		ruido = rnorm(length(m), 0,t)
		return(m+ruido)
	}
	noised = noise(m)
	unicos = which(duplicated(m) == FALSE)
	m[-unicos,] = noised[-unicos,]
	return(m)
  }

toHiperbolico = function(objeto, M = 1 , cx = 0, cy = 0, r = 1){
	insiders = apply(objeto,1,function(co)(cx-co[1])^2+(cy-co[2])^2<r^2)
	outers = which(insiders < 1)
	objetoP = matrix(toPolar(objeto[,1],objeto[,2]),nc=2)
	if (length(outers)){
			objetoP[outers,2] = 1
	}
	objetoP[,2] = sapply(objetoP[,2],fishIin,M)
	objetoC = matrix(toCartesian(objetoP[,1],objetoP[,2]),nc=2)
  return(list(objetoC = objetoC,
              objetoP = objetoP))
  }

resetear = function (x,y) {
  assign('xClick_old',x , envir = POI.env)
  assign('yClick_old',y , envir = POI.env)
  }

OnClickMotion <- function(x,y)
{
  xClick <- x
  yClick <- y
  img = get('img',envir = POI.env)
  width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
  height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
  if (exists( 'xClick_old', env = POI.env) & exists( 'yClick_old', env = POI.env)) {
     xClick_old <- get('xClick_old', env = POI.env)
     yClick_old <- get('yClick_old', env = POI.env)
  } else {
     xClick_old <- xClick
     yClick_old <- yClick
  }
  POI = get('POI',envir = POI.env)
  incrementoScale = POI@IncVscale
  POI@newcoords <- (as.numeric(xClick) - as.numeric(xClick_old)) * incrementoScale
  POI@newcoords_1  <- - ((as.numeric(yClick) - as.numeric(yClick_old))) * incrementoScale
  POI@selected <- 1
  assign('POI', POI, envir = POI.env)
  assign('xClick_old', xClick , envir = POI.env)
  assign('yClick_old', yClick , envir = POI.env)
  tkrreplot(img)
}

OnDobleClick <- function(x,y)
{
  img = get('img',envir = POI.env)
  POI = get('POI',envir = POI.env)
  xCoords <- POI@objetoC[,1]
  yCoords <- POI@objetoC[,2]
  usrCoords   <- c(-1.242, 1.242, 1.242, -1.242 )
  xClick <- as.numeric(x)
  yClick <- as.numeric(y)
  width  <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
  height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))
  xMin <- 0
  xMax <- width
  yMin <- 0
  yMax <- height
  rangeX <- usrCoords[2] - usrCoords[1]
  rangeY <- usrCoords[4] - usrCoords[3]
  imgXcoords <-  (xCoords-usrCoords[1])*(xMax-xMin)/rangeX + xMin
  imgYcoords <- (yCoords-usrCoords[3])*(yMax-yMin)/rangeY + yMin
  xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
  yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)
  squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
  indexClosest <- which.min(squared.Distance)
  assign('indexClosest',indexClosest , envir = POI.env)
  print(POI@docs[indexClosest])
  }

OnMiddleClick = function (x,y) {
  POI = get('POI',envir = POI.env)
  img = get('img',envir = POI.env)
  POI@selected <- 0
  POI@newcoords <- 0
  POI@newcoords_1 <- 0
  assign('POI',POI, envir = POI.env)
  tkrreplot(img)
  }

MadjustPlus = function(){
  POI = get('POI',envir = POI.env)
  img = get('img',envir = POI.env)
  POI@M <- POI@M + 1
  POI@newcoords <- 0
  POI@newcoords_1 <- 0
  assign('POI',POI, envir = POI.env)
  tkrreplot(img)
}

MadjustMinus = function(){
   POI = get('POI',envir = POI.env)
   img = get('img',envir = POI.env)
   if(POI@M >= 0) {POI@M <- POI@M - 1 } else { POI@M <- 0}
   POI@newcoords <- 0
   POI@newcoords_1 <- 0
   assign('POI',POI, envir = POI.env)
   tkrreplot(img)
}

IncVadjustPlus = function(){
   POI = get('POI',envir = POI.env)
   img = get('img',envir = POI.env)
   POI@IncVscale <- POI@IncVscale * 0.1
   POI@newcoords <- 0
   POI@newcoords_1 <- 0
   assign('POI',POI, envir = POI.env)
   tkrreplot(img)
}

IncVadjustMinus = function(){
   POI = get('POI',envir = POI.env)
   img = get('img',envir = POI.env)
   POI@IncVscale <- POI@IncVscale * 10
   POI@newcoords <- 0
   POI@newcoords_1 <- 0
   assign('POI',POI, envir = POI.env)
   tkrreplot(img)
}

MouseWheel = function(D){
   if(as.numeric(D)/120 == 1) {MadjustPlus()}
   if(as.numeric(D)/120 == -1) {MadjustMinus()}
}

centrarSalida = function() {
   POI = get('POI',envir = POI.env)
   img = get('img',envir = POI.env)
   POI@objeto[,1] <- rnorm(nrow(POI@objeto), 0, 0.005)
   POI@Pcoords[1,] <- c(0,0)
   POI@newcoords <- 0
   POI@newcoords_1 <- 0
   assign('POI',POI, envir = POI.env)
   tkrreplot(img)
   }

HeavyWavoidCluttering = function(object, value = 3) {
   m = round(object, value)
   duplis = duplicated(m)
   duplisb = duplicated(m,  fromLast=T)
   clust = duplis + duplisb
   unic = which( duplis == FALSE)
   return(list(newobjeto = matrix(object[unic,], ncol = 2), uniques = unic, clusters = clust))
}

H1WavoidCluttering = function(object, value = 3) {
   m = signif(object, value)
   duplis = duplicated(m)
   duplisb = duplicated(m,  fromLast=T)
   clust = duplis + duplisb
   unic = which( duplis == FALSE)
   return(list(newobjeto = matrix(object[unic,], ncol = 2), uniques = unic, clusters = clust))
}


plotPOI = function(POI){
  try(tkdestroy(get('tt', envir = POI.env)), silent = T)
  if (!exists('POI.env')){
     POI.env <<- new.env()
  }
  if (require(tkrplot)) {
   require(tkrplot)
   tt <- tktoplevel()
   img <- tkrplot(tt,function() {POIPlot(POI)}, hscale = POI@hscale , vscale = POI@vscale )
   assign('tt', tt, envir = POI.env)
   assign('img', img, envir = POI.env)
   rm(img)
   rm(tt)
   tkgrid(get('img',envir = POI.env)) 
   tkbind(get('img',envir = POI.env), "<B1-Motion>", OnClickMotion )  
   tkbind(get('img',envir = POI.env), "<Motion>", resetear)
   tkbind(get('img',envir = POI.env), "<Double-Button-1>",OnDobleClick)
   tkbind(get('img',envir = POI.env), "<Button-2>",OnMiddleClick)
   tkbind(get('tt',envir = POI.env),'+', MadjustPlus)
   tkbind(get('tt',envir = POI.env),'-', MadjustMinus)
   tkbind(get('tt',envir = POI.env),'0', IncVadjustPlus)
   tkbind(get('tt',envir = POI.env),'.', IncVadjustMinus)
   tkbind(get('tt',envir = POI.env),'<MouseWheel>', MouseWheel )
   } else {
     POIPlot(POI)
  }
}


plotPOIGraph = function(POI){
  try(tkdestroy(get('tt', envir = POI.env)), silent = T)
  if (!exists('POI.env')){
     POI.env <<- new.env()
  }
  if (require(tkrplot)) {
    require(tkrplot)
    tt <- tktoplevel()
    img <- tkrplot(tt,function() {POIPlot(POI)}, hscale = POI@hscale , vscale = POI@vscale)
    assign('tt', tt, envir = POI.env)
    assign('img',img, envir = POI.env)
    rm(img)
    rm(tt)
    tkgrid(get('img',envir = POI.env))
    tkbind(get('img',envir = POI.env), "<B1-Motion>", OnClickMotion )
    tkbind(get('img',envir = POI.env), "<Motion>", resetear)
    tkbind(get('img',envir = POI.env), "<Double-Button-1>",OnDobleClick)
    tkbind(get('img',envir = POI.env), "<Button-2>",OnMiddleClick)
    tkbind(get('tt',envir = POI.env),'+', MadjustPlus)
    tkbind(get('tt',envir = POI.env),'-', MadjustMinus)
    tkbind(get('tt',envir = POI.env),'0', IncVadjustPlus)
    tkbind(get('tt',envir = POI.env),'.', IncVadjustMinus)
    tkbind(get('tt',envir = POI.env),'<MouseWheel>', MouseWheel )    
  } else {
     POIPlot(POI)
  }
}

POICreate = function(type = 'POI',...) { 
   new(type,...)
   }
