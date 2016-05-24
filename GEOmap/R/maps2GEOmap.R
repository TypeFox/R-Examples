maps2GEOmap<-function(zz, wx=1, mapnam="temp" )
{
  if(missing(mapnam)) mapnam="temp"
  if(missing(wx)) wx = c(which(is.na(zz$x))) 

  
  #########   convert maps information from a maps data base to a GEOmap data
  Kx = vector()
  Ky = vector()
  Knum = vector()
  Knam = vector()
  Kindex = vector()
  Kstyle = vector()
  Kcode = vector()
  Kcol = vector()

  a1 = 1
  
   start = 0

  if(length(wx)<1)
    {
ibegs = c(1)
iends = c(length(zz$x))
}
  else
  {
ibegs = c(1, wx+1)
iends = c(wx-1, length(zz$x))
  }
  
  n = length(ibegs)
  for(i in 1:n)
    {
      i1 = ibegs[i]
      i2 = iends[i]
     ## print(c(i1, i2))
      el = length(i1:i2)
      Kx = c(Kx, zz$x[i1:i2])
      Ky = c(Ky, zz$y[i1:i2])

      
      Knum[i] = el
      Kindex[i] =  start
      start = start + el
      
      Knam[i] = paste(mapnam, i, sep="_")
      Kstyle[i] = 2
      Kcode[i] =  'o'
      Kcol[i] =  'black'
      
    }

    NEWMAP = list(STROKES = list(nam =Knam , num =Knum , index =Kindex ,
                  col = Kcol, style = Kstyle, code = Kcode, LAT1 = NULL, LAT2 = NULL,
                  LON1 = NULL, LON2 = NULL), POINTS = list(lat = Ky,
                                               lon = Kx))

  NEWMAP =  boundGEOmap(NEWMAP )

  return(NEWMAP )
  

}

