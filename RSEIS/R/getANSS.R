getANSS<-function(fn, skip=2)
  {
    ###   convert an event catalog file downloaded from the ANSS site
    ###  convert data to a list

   ###     http://www.ncedc.org/anss/catalog-search.html

    ###   skip two first lines

    if(file.exists(fn))
      {
        
        send1 = scan(file=fn, what="", sep="\n", skip=skip)
      }
    else
      {
        print(paste( "Cannot find file:", fn) )
        return(NULL)
        
      }

    

    labs = unlist( strsplit( send1[1], split=" ") )
    labs = labs[labs!=""]

    send1 = send1[3:length(send1)]
    s1 = strsplit( send1, split=" ")

    anss = list(yr=0, dom=0, mo=0, hr=0, mi=0, sec=0, lat=0, lon=0, z=0, mag=0)

    for(i in 1:length(s1))
      {

        h = s1[[i]]
        h = h[h!=""]
        
                                        #dtim = h[1]
        d0 = c(h[1], h[2])

        d1  = as.numeric( unlist(strsplit(d0[1], split="/")) )
        t1 = as.numeric( unlist(strsplit(d0[2], split=":")) )
        lat = as.numeric(h[3])
        lon = as.numeric(h[4])
        z = as.numeric(h[5])
        mag = as.numeric(h[6])


        anss$yr[i] = d1[1]
        anss$mo[i] = d1[2]
        anss$dom[i] = d1[3]
        anss$hr[i] = t1[1]
        anss$mi=t1[2]
        anss$sec[i] = t1[3]
        anss$lat[i] = lat
        anss$lon[i] = lon
        anss$mag[i] = mag
        anss$z[i] = z
      }
    return(anss)	
  }


