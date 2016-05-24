`choosfilt` <-
  function(thefilts=thefilts, ncol=5)
  {
    if(missing(ncol)) ncol = 5
###  choosfilt()
    
    if(missing(thefilts))
      {
        thefilts =
          list(flo=
               c(0.02, 0.02, 0.02, 0.02, 0.02,   0.02,
                 0.02, 0.02, 0.02,  0.02, 0.02,  0.02,
                 0.02,
                 1/2, 1/50,1/100, 1/100,1,1,
                 0.2, 15, 5, 2,1,
                 100),
               fhi=
               c(1/10, 1/6, 1/5, 1/4, 1/3, 1/2,
                 0.2,  0.5, 1.0,  2.0, 3.0,  4.0,
                 7.0,
                 8, 1/2.0,1/5.0,1/10.0,10,5,
                 7.0, 100, 100, 100,10,
                 100),
               type =
               c("LP","LP", "LP", "LP", "LP", "LP",
                 "LP","LP", "LP", "LP", "LP", "LP",
                 "LP",
                 "BP", "BP","BP","BP","BP","BP",
                 "HP", "HP","HP", "HP","HP",
                 "None"))
      }


    if(is.null(thefilts))
      {
        thefilts =
          list(flo=
               c(0.02, 0.02, 0.02, 0.02, 0.02,   0.02,
                 0.02, 0.02, 0.02,  0.02, 0.02,  0.02,
                 0.02,
                 1/2, 1/50,1/100, 1/100,1,1,
                 0.2, 15, 5, 2,1,
                 100),
               fhi=
               c(1/10, 1/6, 1/5, 1/4, 1/3, 1/2,
                 0.2,  0.5, 1.0,  2.0, 3.0,  4.0,
                 7.0,
                 8, 1/2.0,1/5.0,1/10.0,10,5,
                 7.0, 100, 100, 100,10,
                 100),
               type =
               c("LP","LP", "LP", "LP", "LP", "LP",
                 "LP","LP", "LP", "LP", "LP", "LP",
                 "LP",
                 "BP", "BP","BP","BP","BP","BP",
                 "HP", "HP","HP", "HP","HP",
                 "None"))



      }


    gfl = c("flo", "fhi", "type")

    
    nfl = names(thefilts)

   m =  match(gfl, nfl )

    if(any(is.na(m)))
      {

        print("ERROR: Wrong Input: check the filter list for flo, fhi and type" )
 
        GIVE = list(type="None",  fl=100,fh=100,ON=FALSE,  proto="BU" )

  
    return(GIVE)


      }


if( length(thefilts$type) != length(thefilts$flo ) &   length(thefilts$type)  != length(thefilts$fhi ) )
  {
    print(paste(sep=' ',"problem with filter definition", length(thefilts$type),length(thefilts$flo ), length(thefilts$fhi )   ))
   print("ERROR: Wrong Input: check the filter list for flo, fhi and type" )
 
        GIVE = list(type="None",  fl=100,fh=100,ON=FALSE,  proto="BU" )

  
    return(GIVE)


  }
    
###  print(data.frame(cbind(thefilts$flo, thefilts$fhi,   1/thefilts$flo, 1/thefilts$fhi,  thefilts$type)))
    
###    namcols = c("springgreen2",   "plum2" ,         "cyan3"  ,        "darkgoldenrod2")
###    match(namcols, colors())
   ###  print(Z)
    colpals = matrix(c(0,238,118,238,174,238,0,205,205,238,173,14)/255, nrow=3)

    colrgb = rgb(colpals[1,],colpals[2,], colpals[3,])
    
    N = length( thefilts$flo)

    utyp  = unique(thefilts$type)
    thecols = rep(colrgb[1], N)

   
       mcol =    match(thefilts$type, utyp )

     thecols =  colrgb[mcol]
          
   ### get(getOption("device"))()
    dev.new()
    
   ### X11()
    

    cols = topo.colors(1.5*N)
   
    nrow = round((N/ncol)+.5)
    
    dx = 1/ncol
    dy =  1/nrow
 plot(c(0,1), c(0,1), type='n', axes=FALSE, xlab='', ylab='')


mtyp = match(thefilts$type,  utyp)
          
lolab = paste(sep=' ', thefilts$flo, "Hz")
          
lolab[thefilts$flo<1] = paste(sep=' ', 1/thefilts$flo[thefilts$flo<1], "s")

hilab = paste(sep=' ', thefilts$fhi, "Hz")
hilab[thefilts$fhi<1] = paste(sep=' ', 1/thefilts$fhi[thefilts$fhi<1], "s") 

 B =  RPMG::itoxyz(1:N, ncol, nrow, 1)
    x = (B$ix-1)*dx
        y = (B$iy-1)*dy
        rect(x , y , x+dx, y+dy, lty=1, col=thecols )

 lab = paste(sep='\n',thefilts$type )
 lab[thefilts$type=="LP"] = paste(sep='\n',thefilts$type[thefilts$type=="LP"],  hilab[thefilts$type=="LP"] )
 lab[thefilts$type=="HP"] = paste(sep='\n',thefilts$type[thefilts$type=="HP"],  lolab[thefilts$type=="HP"] )

  lab[thefilts$type=="BP"] =paste(sep='\n',thefilts$type[thefilts$type=="BP"],  lolab[thefilts$type=="BP"],  hilab[thefilts$type=="BP"] )
     text(x+dx/2, y+dy/2, lab)      
   
    z = locator(n=1, type='p')

    if(length(z$x)<1)
      {
       
        dev.off(dev.cur())
        return(NULL)
        
      }


    
    ii = 1+floor(z$x/dx)
    jj = 1+floor(z$y/dy)
    B =  ii+(jj-1)*(ncol)

    i = B[length(B)]

    GIVE = list(type=thefilts$type[i],  fl=thefilts$flo[i],fh=thefilts$fhi[i],ON=FALSE,  proto="BU" )

    
    print(paste(sep=" ", "DATA is FILTERED", paste(collapse=" ", unlist(GIVE))))

    dev.off(dev.cur())
    
    return(GIVE)

    
########## data.frame(cbind( formatC(thefilts$flo, digits =5, width=6, flag=" "),
##########     formatC(thefilts$fhi, digits =5, width=6, flag=" ") ,
##########     formatC(1/thefilts$flo, digits =5, width=6, flag=" "),
##########     formatC(1/thefilts$fhi, digits =5, width=6, flag=" ") ))


  #####   ans=readline("which filter do you want? ")

   #####  ians = match(ans, 1:length(thefilts$flo))
    
   #####  return(list(ON=FALSE, fl=thefilts$flo[ians], fh=thefilts$fhi[ians], type=thefilts$type[ians], proto="BU"))

  }

