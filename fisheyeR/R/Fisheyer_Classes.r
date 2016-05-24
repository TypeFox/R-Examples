# PACKAGE fisheyeR
# Classes

 setClass(
    Class = 'POI',
    representation(
      matrizSim = 'matrix',
  		cos.query.docs = 'vector',
  	  wordsInQuery = 'ANY',
 			docs = 'matrix',
      objeto = 'matrix',
      objetoC = 'matrix',
      Pcoords = 'matrix',
      PcoordsFI = 'matrix',
      newPcoords = 'matrix',
      newcoords = 'numeric' ,
      newcoords_1 = 'numeric',
      M = 'numeric',
      poisTextCol = 'character' ,
      colores = 'vector' ,
      poisCircleCol = 'character' ,
      linesCol = 'character',
      itemsCol = 'character',
      LABELS =  'logical',
      vscale = 'numeric',
      hscale = 'numeric',
      circleCol = 'character',
      plotCol = 'character',
      itemsFamily = 'character',
      pal = 'character',
      selected = 'numeric' ,
      circRadio = 'numeric' ,
      IncVscale = 'numeric',
      cgnsphrFont = 'numeric',
      xClick_old = 'numeric',
      yClick_old = 'numeric',
      wordsInQueryFull = 'character',
      clustered = 'logical'
      ),
     prototype(cos.query.docs = 0, 
               colores = 0,
               newcoords = 0,
               newcoords_1 = 0,
               M = 3,
               vscale = 1.25 ,
               hscale = 1.25 ,
               circleCol = 'white' ,
               itemsCol = 'white',
               poisTextCol =  '#fff5ee',
               poisCircleCol = '#fff5ee',
               linesCol = 'white',
               plotCol = 'black',
               itemsFamily = 'sans',
               pal = 'topo' ,
               selected = 1 ,
               circRadio = 0.25  ,
               IncVscale = 0.005  ,
               cgnsphrFont = 1.01,
               LABELS = TRUE ,
               clustered = TRUE
               )
 )

 setClass(
    Class = 'POIGraph',
    representation(EDGES = 'matrix'),
    prototype(LABELS = FALSE),
    contains = 'POI')

 setClass(
    Class = 'multiPOI',    
    contains = 'POI')

 setClass(
    Class = 'mPOIAnd',    
    contains = 'multiPOI')
    
 setClass(
    Class = 'mPOIOr',    
    contains = 'multiPOI') 

    