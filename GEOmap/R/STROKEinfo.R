STROKEinfo<-function(map, w=1, h=NULL)
  {
    ###  source("STROKEinfo.R")
    ############   example:  STROKEinfo(haiti.map, h=c("nam", "col") )
    ############   example:  STROKEinfo(haiti.map, w=5:10, h=c("nam", "col") )
    ############   example:  STROKEinfo(haiti.map, w="Car", h=c("nam", "col") )
    
    if(missing(w))
      {

        w = 1:length(map$STROKES$num)
        
      }
    NAMS = names(map$STROKES)
    if(missing(h))
      {
        h = 1:length(NAMS)
        cnum = h
      }
    else
      {

        if(is.character(h))
          {
            cnum  = match(h, NAMS)
          }
        else
          {
            cnum = h

          }

      }


     if(is.character(w))
          {
            rnum  = grep(w, map$STROKES$nam)
          }
        else
          {
            rnum = w

          }

    
    D1 = data.frame(map$STROKES)

    if(is.numeric(rnum))
      {
        return(D1[rnum,cnum])

      }

    
  }

