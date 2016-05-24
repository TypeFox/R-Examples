`printGEOinfo` <-
function(MAP, kstroke)
  {


df = data.frame(id=kstroke, nam=MAP$STROKES$nam[kstroke],
  index=MAP$STROKES$index[kstroke],
  num=MAP$STROKES$num[kstroke],
  col=MAP$STROKES$col[kstroke],
  style=MAP$STROKES$style[kstroke],
   code=MAP$STROKES$code[kstroke]
  )

  print(df)
  

  }

