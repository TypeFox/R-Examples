getIRIS<-function(fn, skip=0)
  {


###   convert an event catalog file downloaded from the IRIS site
    ###  convert data to a list

    if(file.exists(fn))
      {
        
        send1 = scan(file=fn, what="", sep="\n", skip=skip)
      }
    else
      {
        print(paste( "Cannot find file:", fn) )
        return(NULL)
        
      }

   
    s1 = strsplit( send1, split="\t")

    jim = list(yr=0, dom=0, mo=0, hr=0, mi=0, sec=0, lat=0, lon=0, z=0, mag=0)

    for(i in 1:length(s1))
      {

        h = s1[[i]]
        dtim = h[1]
        d0 = unlist(strsplit(split=" ", dtim))

        d1  = as.numeric( unlist(strsplit(d0[1], split="/")) )
        t1 = as.numeric( unlist(strsplit(d0[2], split=":")) )
        lat = as.numeric(h[2])
        lon = as.numeric(h[3])
        z = as.numeric(h[4])
        mag = as.numeric(h[5])


        jim$yr[i] = d1[1]
        jim$mo[i] = d1[2]
        jim$dom[i] = d1[3]
        jim$hr[i] = t1[1]
        jim$mi=t1[2]
        jim$sec[i] = t1[3]
        jim$lat[i] = lat
        jim$lon[i] = lon
        jim$mag[i] = mag
        jim$z[i] = z
      }


    return(jim)	



  }

