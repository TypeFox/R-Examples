`Markup` <-
  function(MM=list(), sel=1, cex=1, ...)
  {

    if(missing(cex)) { cex=1 }
    if(missing(sel)) { sel = 1:length(MM) }

encoding="ISOLatin1"
myscreenfams = c("serif", "sans")
mypsfams = c("NimbusRom", "URWHelvetica")

    
   ###   ps.options(encoding=encoding, family=family)

    ps.options(encoding=encoding, family=mypsfams[1])

   
    for(i  in  sel)
      {
        if(is.null(MM[[i]]$CEX))
          { charex = cex }
        else
          {
            charex = MM[[i]]$CEX
          }

        ARRTF = MM[[i]]$ARR
        if(identical(ARRTF , "F"))ARRTF=FALSE
        if(identical(ARRTF , "T"))ARRTF=TRUE
        if(identical(ARRTF , "f"))ARRTF=FALSE
        if(identical(ARRTF , "t"))ARRTF=TRUE

        ROTTF = MM[[i]]$ROT

        if(identical(ROTTF , "F"))ROTTF=FALSE
        if(identical(ROTTF , "T"))ROTTF=TRUE
        if(identical(ROTTF , "f"))ROTTF=FALSE
        if(identical(ROTTF , "t"))ROTTF=TRUE

        
        
        if(MM[[i]]$pos<=0) { POS = NULL;  ADJ = c(0,0)} else { POS =MM[[i]]$pos ;  ADJ = NULL  }
        
        if(ARRTF==TRUE)
          {
            arrows(MM[[i]]$x1,MM[[i]]$y1 , MM[[i]]$x2,MM[[i]]$y2, ...)
          }
        
        if(ROTTF==TRUE)
          {
            if(MM[[i]]$pos==0)
              {
                MM[[i]]$adj = c(0,0)
                text(MM[[i]]$x1,MM[[i]]$y1, labels=MM[[i]]$lab, adj=ADJ,  srt=MM[[i]]$angdeg, cex=charex, font=MM[[i]]$font)
              }
            else
              {
                text(MM[[i]]$x1,MM[[i]]$y1, labels=MM[[i]]$lab,   pos=POS, srt=MM[[i]]$angdeg, cex=charex, font=MM[[i]]$font)
              }
          }
        else
          {
            
            text(MM[[i]]$x1,MM[[i]]$y1, labels=MM[[i]]$lab, adj=ADJ,   pos=POS, cex=charex, font=MM[[i]]$font)
          }
        

        
      }
    
  }

