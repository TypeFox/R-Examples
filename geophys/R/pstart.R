pstart <-
function(xscale=30, expand=1.2)
  {
    if(missing(xscale)) xscale=30
    if(missing(expand)) expand=1.2

    
   ##  print(paste(sep=" ", "xscale=", xscale))

    
    plot(c(-expand*xscale, expand*xscale), c(-expand*xscale, expand*xscale), type='n' ,  axes=FALSE, ann=FALSE, asp=1)
  }

