clusterWPX<-function(twpx, tol=200, PLOT=FALSE)
  {
####  Given a set of pix,
####  separate and save them as a list of
####  new wpx dataframes that
####  can be stored

####     if all the pix are within tol, just return

####jj  = RSEIS::swig(GH, sel=GH$COMPS=="V", PADDLAB="YPIX" ); twpx = jj$g$WPX
    ##  clusterWPX(twpx)
    


    nona = !is.na(twpx$tag)
    
    twpx = twpx[nona,]
    
    A1T = Qrangedatetime(twpx)
    s1 = RSEIS::secdifL(A1T$min,  twpx)
    

    D1 = dist(s1)

    if(all(D1<tol)) return(list(MEM=twpx))

########   force points that are less than the tolerance
########   distance to be zero 
    D1[D1<tol] = 0

    
###    H2 = hclust(D1)
    H2 = hclust(D1,  method = "ward" )

    b = boxplot(H2$height, plot=FALSE)


    if(length(b$out)<1)
      {
        jheight = max(H2$height)
      }
    else
      {
        jheight = mean( c( min( b$out) , b$stats[5,])) 
      }
    
    hcut =  cutree(H2, h=jheight )
    

    
    if(PLOT)
      {
        plot(H2)

      }
    
    icut = unique(hcut)

    LI = length(icut)
    KL = vector(mode="list", length=LI)

    for(i in 1:LI )
      {
        j = icut[i]

        wj = which(hcut==j)
        
        KL[[i]] = twpx[wj,]

##### saveWPX(clusx, destdir="." )

      }

    return(KL)
    
  }
