`parseFN2STA` <-
function(fn, ista, icomp, sep="\\.", dir=0)
  {
    ###  get the station name and the component name from the file name
    ### function default assumes that the station name and the component name
      ###  are the last  items on the finle name seperated by a period
    ## first get the file name from the full path name

    ##  ista and icomp are computed from the end of the file name

    ##  if dir>0 then the counting is from the beginning starting with 1

    if(missing(ista)) { ista = 1 }
    if(missing(icomp)) { icomp = 0 } 
    if(missing(sep)) { sep="\\." } 


    fn1 = basename(fn)
    
    f2 = unlist(strsplit(fn1, split=sep))
    
    if(dir<1)
      {
        sta = f2[length(f2)-ista]
        comp = f2[length(f2)-icomp]
      }
    else
      {
        sta = f2[ista]
        comp = f2[icomp]
      }
    
    return(list(sta=sta, comp=comp) )

    
  }

