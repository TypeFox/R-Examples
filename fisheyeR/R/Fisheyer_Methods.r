# PACKAGE fisheyeR
# Methods
require(methods)

 setGeneric("POIcalculate<-" , function(object, value){standardGeneric ("POIcalculate<-")})

 setReplaceMethod( f ="POIcalculate",
                   signature = 'POI',
                   definition = function(object, value){
                                   object@matrizSim <- value$XhihatQ
                                   object@cos.query.docs <- value$cos.query.docs
                                   object@wordsInQuery <- value$wordsInQuery
                                   object@docs <- value$docs
                                   return(object)
                                }
 )

 setReplaceMethod( f ="POIcalculate",
                   signature = 'multiPOI',
                   definition = function(object, value){
                                   object@matrizSim <- value$XhihatQ
                                   object@cos.query.docs <- value$cos.query.docs
                                   object@wordsInQuery <- value$wordsInQuery
                                   object@docs <- value$docs
                                   object@wordsInQueryFull <- value$wordsInQueryFull
                                   return(object)
                                }
 )

 setGeneric("POIcolors<-" , function(object, value){standardGeneric ("POIcolors<-")})

 setReplaceMethod( f ="POIcolors",
                   signature = 'POI',
                   definition = function(object, value){
                                   object@colores <- value
                                   return(object)
                                }
 )

 setGeneric("query2Cols" ,
             function(object, value){standardGeneric ("query2Cols")})

 setMethod("query2Cols" ,
            signature = "POI",
            function(object, value){

   cosquery = object@cos.query.docs
   paleta = value
   cosquery[cosquery < 0] <- 0
   cosquery[cosquery > 1] <- 1
	 z <- 100 * matrix(cosquery)        
	 x <- 1:nrow(z)
	 y <- 1:ncol(z)
	 zlim <- c(0,100) 
	 zlen <- round(zlim[2] - zlim[1])
	 if (paleta == 'heat'){
	  	colorlut <- heat.colors(zlen)
	 }
	 if (paleta == 'topo'){
  		colorlut <- topo.colors(zlen)
 	 }
	 if (paleta == 'cm'){
	  	colorlut <- cm.colors(zlen)
	 }
	 if (paleta == 'terrain'){
  		colorlut <- terrain.colors(zlen)
	 }
	 z[z < 1] <- 1
	 col <- colorlut[ z-zlim[1] ]
   return(col)
   }
 )

 setGeneric("POIcoords<-" , function(object, value){standardGeneric ("POIcoords<-")})

 setReplaceMethod( f ="POIcoords",
                   signature = 'POI',
                   definition = function(object, value){
                                   object@Pcoords <- value$Pcoords
                                   object@PcoordsFI <- value$PcoordsFI
                                   object@newPcoords <- value$newPcoords
                                   object@objeto <- value$objeto
                                   return(object)
                                }
 )

 setGeneric("POICalc" ,
             function(objeto, NC, cx=0, cy=0, r=1, ...){standardGeneric ("POICalc")})

 setMethod("POICalc" ,
            signature = "POI",
            function(objeto, NC, cx=0, cy=0, r=1, ...){

   MatrizSim = objeto@matrizSim
   secuencia = seq(2/NC,2,2/NC)
   Pcoords = matrix(rep(0,NC*2),nc=2)
   n = 1
   for (i in secuencia){
      Pcoords[n,] = c(r * cos(i*pi), r * sin(i*pi))
      n = n+1
   }
   PcoordsFI = matrix(toPolar(Pcoords[,1],Pcoords[,2]),nc=2)
   PcoordsFI[,2] = PcoordsFI[,2]+.15
   PcoordsFI = matrix(toCartesian(PcoordsFI[,1],PcoordsFI[,2]),nc=2)

   if (nrow(Pcoords) != 1){
      newPcoords = puntosMedios(Pcoords)
   } else {
      newPcoords = Pcoords
   }

   MatrizSim[is.nan(MatrizSim/rowSums(MatrizSim))] <- 0

   W = MatrizSim / rowSums(MatrizSim)
   W[is.nan(W)] <- 0
   nwords = nrow(W)
   objeto = matrix(rep(0,2*nwords),nc=2)
   for (j in 1:nwords){
      for (nPOI in 1:NC){
         objeto[j,1] = objeto[j,1]+(W[j,nPOI]*Pcoords[nPOI,1])
         objeto[j,2] = objeto[j,2]+(W[j,nPOI]*Pcoords[nPOI,2])
      }
   }

   objeto = addNoise(objeto)

   return(list(Pcoords = Pcoords,
               PcoordsFI = PcoordsFI,
               newPcoords = newPcoords,
               objeto = objeto))

   }
 )

 setGeneric("POIPlot" ,
             function(POI){standardGeneric ("POIPlot")})

 setMethod("POIPlot" ,
            signature = "POI",
            function(POI){

   par(bg=POI@plotCol, mar = c(0.1,0.1,0.1,0.1), family = POI@itemsFamily)

   if (exists('POI.env')) {
      if (exists('POI', envir = POI.env)) {
        POI <- get('POI', envir = POI.env)
      }
   }
   
   clustered = POI@clustered 
   selected = POI@selected
   objeto = POI@objeto
   newcoords = POI@newcoords
   newcoords_1 = POI@newcoords_1
   NC = length(POI@wordsInQuery)
   cx=0
   cy=0
   r=1
   etiq2 = POI@docs[,1]
   etiq = POI@wordsInQuery
   fishEYE = TRUE
   M = POI@M
   poisTextCol = POI@poisTextCol
   colores = POI@colores[as.numeric(POI@docs[,2])]
   poisCircleCol = POI@poisCircleCol
   linesCol = POI@linesCol
   itemsCol = POI@itemsCol
   circleCol = POI@circleCol
   LABELS =  POI@LABELS
   Pcoords = POI@Pcoords
   newPcoords = POI@newPcoords
   cgnsphrFont = POI@cgnsphrFont

   newcoords_par = newcoords

   newcoords_Pcoords = matrix(rep( c(newcoords,newcoords_1 ),
                              nrow(Pcoords)),nc=2,byrow=TRUE)

   newcoords_puntosMediosPcoords = matrix(rep( c(newcoords,newcoords_1),
                                          nrow(newPcoords)),nc=2,byrow=TRUE)

   newcoords = matrix(rep( c(newcoords,newcoords_1),
                      nrow(objeto)),nc=2,byrow=TRUE)

   objeto = objeto+newcoords
   objetoH = toHiperbolico(objeto, M)
   objetoC = objetoH$objetoC
   objetoP = objetoH$objetoP

   Pcoords = Pcoords + newcoords_Pcoords
   PcoordsH = toHiperbolico(Pcoords, M)
   PcoordsC = PcoordsH$objetoC
   PcoordsP = PcoordsH$objetoP

   newPcoords = newPcoords + newcoords_puntosMediosPcoords
   newPcoordsH = toHiperbolico(newPcoords, M)
   Pcoords_objetoC = newPcoordsH$objetoC

   if (LABELS) {
      PcoordsFI = matrix(toPolar(PcoordsC[,1],PcoordsC[,2]),nc=2)
      PcoordsFI[,2] = 1 +.15
      PcoordsFI = matrix(toCartesian(PcoordsFI[,1],PcoordsFI[,2]),nc=2)
   }

   if (clustered == TRUE) {
      objetoC.bk <- objetoC
      newobject <- H1WavoidCluttering(objetoC, 2)
      objetoC  <- newobject$newobjeto
      etiq2 <- etiq2[newobject$uniques]

      if (sum(is.na(itemsCol[newobject$uniques])) == 0) {
         itemsCol <- itemsCol[newobject$uniques]
      }
      colores <- colores[newobject$uniques]
      objetoP <- matrix(objetoP[newobject$uniques,], nc = 2)
   
      pchs = newobject$uniques %in% which(newobject$clusters == 0)
      pchs[pchs > 0] <- 19
      pchs[pchs == 0] <- 23
      Pminus = (newobject$uniques %in% which(newobject$clusters != 0))/2
   } else {
      Pminus = 0
      pchs <- 19
   }
   
   plot(circulo(0,0,1, circleCol, PLOT = FALSE),cex=.5,ylim=c(-1.15,1.15),xlim=c(-1.15,1.15),
                  ann=FALSE, axes=F,type='l', col = circleCol)

   points(objetoC, 
          pch = pchs, col = colores, bg = colores, cex = 1.5 - (objetoP[,2] - Pminus) )

   text(objetoC[,1], objetoC[,2], labels = etiq2, cex = cgnsphrFont - objetoP[,2],
        pos = 3, col = itemsCol)

   abline(h = cx, col = 'grey', lty = 'dashed')
   abline(v = cy, col = 'grey', lty = 'dashed')

   points(PcoordsC, cex = 2, col = poisCircleCol)

   lines(Pcoords_objetoC, col = linesCol)

   segments(Pcoords_objetoC[nrow(Pcoords_objetoC),1],Pcoords_objetoC[nrow(Pcoords_objetoC),2],
            Pcoords_objetoC[1,1],Pcoords_objetoC[1,2], col = linesCol)

   if (LABELS) {
      text(PcoordsFI[,1],PcoordsFI[,2],toupper(etiq),cex=.75, col = poisTextCol)
   }

   if (clustered == TRUE) {
      objetoC <- objetoC.bk
   }
   
   if (selected != 1) {
      circulin(0,0, POI@circRadio, objeto = objetoC)
   }

   if (!exists('POI.env')){
      POI.env <<- new.env()
   }

   poiCOPY = POI
   poiCOPY@objeto <- objeto
   poiCOPY@objetoC <- objetoC
   poiCOPY@newPcoords <- newPcoords
   poiCOPY@Pcoords <- Pcoords
   assign('POI',poiCOPY , envir = POI.env)

   }
 )

 setMethod("POIPlot" ,
            signature = "multiPOI",
            function(POI){

   par(bg=POI@plotCol, mar = c(0.1,0.1,0.1,0.1), family = POI@itemsFamily)

   if (exists('POI.env')) {
      if (exists('POI', envir = POI.env)) {
        POI <- get('POI', envir = POI.env)
      }
   }

   clustered = POI@clustered
   selected = POI@selected
   objeto = POI@objeto 
   newcoords = POI@newcoords 
   newcoords_1 = POI@newcoords_1  
   NC = length(POI@wordsInQuery) 
   cx=0
   cy=0
   r=1
   etiq2 = POI@docs[,1]
   etiq = POI@wordsInQuery 
   fishEYE = TRUE
   M = POI@M 
   poisTextCol = POI@poisTextCol 
   if (class(POI) == 'mPOIAnd') {
      colores = POI@colores[as.numeric(POI@docs[,2])] 
   } else {
      colores = POI@colores
   }
   poisCircleCol = POI@poisCircleCol 
   linesCol = POI@linesCol 
   itemsCol = POI@itemsCol 
   circleCol = POI@circleCol
   LABELS =  POI@LABELS 
   Pcoords = POI@Pcoords
   newPcoords = POI@newPcoords
   cgnsphrFont = POI@cgnsphrFont

   newcoords_par = newcoords

   newcoords_Pcoords = matrix(rep( c(newcoords,newcoords_1 ),
                              nrow(Pcoords)),nc=2,byrow=TRUE)

   newcoords_puntosMediosPcoords = matrix(rep( c(newcoords,newcoords_1),
                                          nrow(newPcoords)),nc=2,byrow=TRUE)

   newcoords = matrix(rep( c(newcoords,newcoords_1),
                      nrow(objeto)),nc=2,byrow=TRUE)

   objeto = objeto+newcoords
   objetoH = toHiperbolico(objeto, M)
   objetoC = objetoH$objetoC
   objetoP = objetoH$objetoP

   Pcoords = Pcoords + newcoords_Pcoords
   PcoordsH = toHiperbolico(Pcoords, M)
   PcoordsC = PcoordsH$objetoC
   PcoordsP = PcoordsH$objetoP

   newPcoords = newPcoords + newcoords_puntosMediosPcoords
   newPcoordsH = toHiperbolico(newPcoords, M)
   Pcoords_objetoC = newPcoordsH$objetoC

   if (LABELS) {
      PcoordsFI = matrix(toPolar(PcoordsC[,1],PcoordsC[,2]),nc=2)
      PcoordsFI[,2] = 1 +.15
      PcoordsFI = matrix(toCartesian(PcoordsFI[,1],PcoordsFI[,2]),nc=2)
   }

   if (clustered == TRUE) {
      objetoC.bk <- objetoC
      newobject <- H1WavoidCluttering(objetoC, 2)
      objetoC  <- newobject$newobjeto
      etiq2 <- etiq2[newobject$uniques]

      if (sum(is.na(itemsCol[newobject$uniques])) == 0) {
         itemsCol <- itemsCol[newobject$uniques]
      }
      colores <- colores[newobject$uniques]
      objetoP <- matrix(objetoP[newobject$uniques,], nc = 2)
   
      pchs = newobject$uniques %in% which(newobject$clusters == 0)
      pchs[pchs > 0] <- 19
      pchs[pchs == 0] <- 23
      Pminus = (newobject$uniques %in% which(newobject$clusters != 0))/2
   } else {
      Pminus = 0
      pchs <- 19
   }

   plot(circulo(0,0,1, circleCol, PLOT = FALSE),cex=.5,ylim=c(-1.15,1.15),xlim=c(-1.15,1.15),
                  ann=FALSE, axes=F,type='l', col = circleCol)

   points(objetoC, pch = pchs, bg = colores, col = colores, cex = 1.5 - (objetoP[,2] - Pminus) )

   text(objetoC[,1], objetoC[,2], labels = etiq2, cex = cgnsphrFont - objetoP[,2],
        pos = 3, col = itemsCol)

   abline(h = cx, col = 'grey', lty = 'dashed')
   abline(v = cy, col = 'grey', lty = 'dashed')

   points(PcoordsC,cex = 2, col = poisCircleCol)

   lines(Pcoords_objetoC, col = linesCol)

   segments(Pcoords_objetoC[nrow(Pcoords_objetoC),1],Pcoords_objetoC[nrow(Pcoords_objetoC),2],
            Pcoords_objetoC[1,1],Pcoords_objetoC[1,2], col = linesCol)

   if (LABELS) {
      text(PcoordsFI[,1],PcoordsFI[,2],toupper(etiq),cex=.75, col = poisTextCol)
   }

   if (clustered == TRUE) {
      objetoC <- objetoC.bk
   }

   if (selected != 1) {
      circulin(0,0, POI@circRadio, objeto = objetoC)   
   }

   if (!exists('POI.env')){
      POI.env <<- new.env()
   }
   poiCOPY = POI
   poiCOPY@objeto <- objeto
   poiCOPY@objetoC <- objetoC
   poiCOPY@newPcoords <- newPcoords
   poiCOPY@Pcoords <- Pcoords
   assign('POI',poiCOPY , envir = POI.env)

   }
 )


 setMethod("POIPlot" ,
            signature = "POIGraph",
            function(POI){
   par(bg=POI@plotCol, mar = c(0.1,0.1,0.1,0.1), family = POI@itemsFamily)
   if (exists('POI.env')) {
      if (exists('POI', envir = POI.env)) {
        POI <- get('POI', envir = POI.env)
      }
   }
   selected = POI@selected
   objeto = POI@objeto
   newcoords = POI@newcoords
   newcoords_1 = POI@newcoords_1
   NC = length(POI@wordsInQuery)
   cx=0
   cy=0
   r=1
   if(length(POI@docs) == 0) {
     etiq2 <- matrix(seq(1:nrow(objeto)))
   } else {
     etiq2 = POI@docs[,1]
   }
   fishEYE = TRUE
   M = POI@M
   poisTextCol = POI@poisTextCol
   colores = POI@colores[as.numeric(POI@docs[,2])]
   poisCircleCol = POI@poisCircleCol
   linesCol = POI@linesCol
   itemsCol = POI@itemsCol
   circleCol = POI@circleCol
   LABELS =  POI@LABELS
   Pcoords = POI@Pcoords
   newPcoords = POI@newPcoords
   cgnsphrFont = POI@cgnsphrFont
   EDGES = POI@EDGES
   newcoords_par = newcoords
   newcoords = matrix(rep( c(newcoords,newcoords_1),
                      nrow(objeto)),nc=2,byrow=TRUE)

   objeto = objeto+newcoords
   objetoH = toHiperbolico(objeto, M)
   objetoC = objetoH$objetoC
   objetoP = objetoH$objetoP
   plot(circulo(0,0,1, circleCol, PLOT = FALSE),cex=.5,ylim=c(-1.15,1.15),xlim=c(-1.15,1.15),
                  ann=FALSE, axes=F,type='l', col = circleCol)
   points(objetoC, pch=19, col = colores, cex = 1.5 - objetoP[,2])
   if (length(EDGES) != 0) {
     EDGES.lines = cbind(objetoC[ EDGES[ ,1] ,], objetoC[ EDGES[ ,2] ,])
     segments(EDGES.lines[,1], EDGES.lines[,2], EDGES.lines[,3], EDGES.lines[,4], col = 'red')
   }
   text(objetoC[,1], objetoC[,2], labels = etiq2, cex = cgnsphrFont - objetoP[,2],
        pos = 3, col = itemsCol)
   abline(h = cx, col = 'grey', lty = 'dashed')
   abline(v = cy, col = 'grey', lty = 'dashed')
   if (selected != 1) {
      circulin(0,0, POI@circRadio, objeto = objetoC)
   }
   if (!exists('POI.env')){
      POI.env <<- new.env()
   }
   poiCOPY = POI
   poiCOPY@objeto <- objeto
   poiCOPY@objetoC <- objetoC
   assign('POI',poiCOPY , envir = POI.env)
   }
 )
