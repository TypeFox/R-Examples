`gaddtix` <-
function(side=3, pos=0,   tck=0.005, at=c(0,1), labels=NULL, col=2, addline=FALSE, ...)
  {
 ##X##   add tick marks to plot   
##X## ###     addtix(side=3, pos=y3+dy,   tck=0.005, at=ttics, labels=FALSE, col=2 )
    if(missing(side)) {  side = 3 }
    if(missing(pos)) {  pos = 0 }
    if(missing(tck)) {  tck = 0 }
    if(missing(at)) {  at = 0 }
    if(missing(labels)) {  labels = NULL }
    if(missing(col)) {  col = "black" }
    if(missing(addline) )  {  addline = FALSE } 
    
    n = length(at)
    if(side==1)
      {
        if(addline==TRUE) lines( c(at[1], at[n]), c(pos,pos) , col=col, ...)
        segments(at,pos, at, pos+tck , col=col, xpd=TRUE)
        if(!is.null(labels))
          {
           ##  mtext(labels, side=1, line=.5, at=at)

            text(at, rep(pos, length(at)), labels=labels, pos=1, xpd=TRUE)
            
          }
      }
    
    if(side==3)
      {
       if(addline==TRUE)  lines( c(at[1], at[n]), c(pos,pos) , col=col, ...)
        segments(at, pos, at, pos-tck , col=col, xpd=TRUE)

       if(!is.null(labels))
         {
        ###   mtext(labels, side=3, line=.5, at=at)

             text(at, rep(pos, length(at)), labels=labels, pos=3, xpd=TRUE)

           
         }
       
      }
    
    
    if(side==2)
      {
        if(addline==TRUE) lines(  c(pos,pos) , c(at[1], at[n]), col=col, ...)
        segments(pos,at,  pos+tck ,at, col=col, xpd=TRUE)


        if(!is.null(labels))
          {
           ###  mtext(labels, side=2, line=.5, at=at)

            text( rep(pos, length(at)), at, labels=labels, pos=2, xpd=TRUE)
            
          }
        
        
      }
    
    if(side==4)
      {
       if(addline==TRUE)  lines(  c(pos,pos) , c(at[1], at[n]), col=col, ...)
        segments( pos, at, pos-tck , at,  col=col, xpd=TRUE)

       if(!is.null(labels))
         {
          ###  mtext(labels, side=4, line=.5, at=at)

           text( rep(pos, length(at)), at, labels=labels, pos=4, xpd=TRUE)
         }
        
       
      }
    


    

  }

