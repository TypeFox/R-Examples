HiJACK<-function(lps, sta, vel)
{


  BSI = vector(mode="list")

  for(ipi in 1:length(lps))
    {
      g1 = RSEIS::getpfile(lps[ipi], stafile = NULL)

      MA = match(g1$STAS$name, sta$name)

      g1$STAS$lat = sta$lat[MA]
      g1$STAS$lon = sta$lon[MA]
      g1$STAS$z = sta$z[MA]

      w1 = which(!is.na(g1$STAS$lat) & !is.na(g1$STAS$lon) )
      Ldat  = LDATlist(g1, w1)

      BSI[[ipi]]   = BLACKJACK(Ldat, vel)


    }
################################

  ZEYE = vector(length=length(sta$name), mode="list")
  XEYE = vector(length=length(sta$name), mode="list")
  YEYE = vector(length=length(sta$name), mode="list")

  names(ZEYE) = sta$name
  names(XEYE) = sta$name
  names(YEYE) = sta$name

  for(k in 1:length(BSI))
    {
      g = BSI[[k]]
      SI = g$SI
      
      MA = match( rownames(SI), names(ZEYE) )
      for(i in 1:length(MA))
        {
          YEYE[[MA[i]]] = c(YEYE[[MA[i]]],SI[i,1] ) 
          XEYE[[MA[i]]] = c(XEYE[[MA[i]]],SI[i,2] ) 
          ZEYE[[MA[i]]] = c(ZEYE[[MA[i]]],SI[i,3] ) 

        }

      
    }

  
  return(list(X = XEYE, Y = YEYE, Z =ZEYE, sta=sta ) )   

}



