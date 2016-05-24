bordear <-
function(wb,hoja,fila=1,columna=1,ancho=1,alto=1,
                 estilo=createStyle(border="topleftbottomright",
                                   borderStyle = "double",borderColour = "blue4")){
  argb=function(x) paste0("#",substr(x,3,8))
  if(!is.null(estilo$borderTop)){
    es=createStyle(border="top",borderColour = argb(estilo$borderTopColour),
                   borderStyle =estilo$borderTop)
    addStyle(wb,hoja,es,rows = fila ,cols=columna+(1:ancho)-1,
             gridExpand = TRUE,stack = TRUE)
  }
  if(!is.null(estilo$borderBottom)){
    es=createStyle(border="bottom",borderColour = argb(estilo$borderBottomColour),
                   borderStyle =estilo$borderBottom)
    addStyle(wb,hoja,es,rows = fila+alto-1 ,cols=columna+(1:ancho)-1,
             gridExpand = TRUE,stack = TRUE)
  }
  if(!is.null(estilo$borderLeft)){
    es=createStyle(border="left",borderColour = argb(estilo$borderLeftColour),
                   borderStyle =estilo$borderLeft)
    addStyle(wb,hoja,es,rows = fila+(1:alto)-1,cols=columna,
             gridExpand = TRUE,stack = TRUE)
  }
  if(!is.null(estilo$borderRight)){
    es=createStyle(border="right",borderColour = argb(estilo$borderRightColour),
                   borderStyle =estilo$borderRight)
    addStyle(wb,hoja,es,rows = fila+(1:alto)-1 ,cols=columna+ancho-1,
             gridExpand = TRUE,stack = TRUE)
  }
}
