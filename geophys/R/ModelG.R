ModelG<-function(Data, labs=c("Done"), obs=NULL, ZCOLS=RPMG::pastel.colors(24, seed=2 ) )
  {

    if(missing(obs)) { obs=NULL }
    ###  Data is the input model

    
    ####### 2.5D  gravity Modeling
    
     if(missing(ZCOLS))
       {
         ZCOLS = RPMG::pastel.colors(24, seed=2 )
       }
   
    #####  
   if(missing(labs)) {
     labs = c("DONE", "REFRESH", "SAVE" , "POLY", "ERASE" , "DELETE", "SHOWpts"  ,"COLOR"  , "MVpt", "MVpoly",  "Rho", "MODEL", "SELpol" )
   }


    whichSELpoly<-function(zloc, pts, TOL)
      {
        npt = length(zloc$x) - 1
        if(npt>0)
          {
            
            SELPOLY = NULL
            for(i in 1:npt)
              {
                x1 = zloc$x[i]
                y1 = zloc$y[i]
                 sp1 = which( sqrt(  (pts$x-x1)^2 + (pts$y-y1)^2)  < TOL )
                if(length(sp1)>0)
                  {
                    SELPOLY = c(SELPOLY, sp1)
                  }
              }
            
            SELPOLY = unique(SELPOLY)
            
            print(SELPOLY) 
            
          }
      }


    
centroid<-function(p)
  {
    n <- length(p$x)
    x1 <- p$x
    i2 <- c(n, 1:(n - 1))
    x2 <- p$x[i2]
    y1 <- p$y
    y2 <- p$y[i2]
    a <- x1 * y2 - x2 * y1
    s <- sum(a) * 3
    if (s == 0) 
      c(mean(x1), mean(y1))
    else c(sum((x1 + x2) * a)/s, sum((y1 + y2) * a)/s)
  }


    create.mod<-function()
      {
        EX = c(0,140000, 0, 140000)
        ZEE = c(0, 0, 100000, 100000)

        Data = vector(mode="list")
        Data$mod = vector(mode="list")
        Data$cens = list(x=NA, y=NA, onoff=NA )
        Data$n = 0
        Data$xmin=EX[1]
        Data$xmax=EX[2]
        Data$zmin=-ZEE[1]
        Data$zmax=-ZEE[3]
        
        return(Data)
      }

delete.poly<-function(Data, w1 )
  {
    NewDATA = create.mod()
    NewDATA$xmin=Data$xmin
    NewDATA$xmax= Data$xmax
   NewDATA$zmin =  Data$zmin
   NewDATA$zmax =  Data$zmax

    g =  1:Data$n
    g = g[-w1]
    NewDATA$n = length(g)
    NewDATA$mod = vector(length=NewDATA$n,    mode="list")

    
    for(k in 1:NewDATA$n )
      {
        NewDATA$mod[[k]] = Data$mod[[g[k]]]
        NewDATA$cens$x[k] = Data$cens$x[g[k]]
        NewDATA$cens$y[k] = Data$cens$y[g[k]]
      }
    return(NewDATA)

  }
    

   DefaultDen = 0.2

   if(is.null(Data$mod))
     {
      Data =  create.mod()
   
 }


   button.color = grey(0.6)
    button.color2 = 'red'

   colabs = rep(button.color, length=length(labs))
   pchlabs = rep(0,length(labs))
    
    wb = which("SELpol" == labs)
    colabs[wb] = button.color2 
   
   NLABS = length(labs)
   NOLAB = NLABS +1000  ## some large number
   ##################
   ##################
   ##################
    TOL = 0.02*(Data$xmax-Data$xmin)
    SELPOLY = NULL
    showpoints = FALSE
    showcolor = FALSE

    tcolON  = rgb(1,.85,.85)
    tcolOFF = rgb(.85,.85, 1)
    
    
##################
    replot<-function(Data)
      {
        plot(c(Data$xmin, Data$xmax), c(Data$zmin, Data$zmax),
             type="n", xlab="X-m", ylab="Depth, m")
        
        grid()
        
        if(Data$n>0)
          {
            for(i in 1:Data$n)
              {
                acol=NA
                pcol=ZCOLS[i]
                if(showcolor) { acol=ZCOLS[i] }
                polygon(c(Data$mod[[i]]$x), c(Data$mod[[i]]$y  ) , col=acol )
                if(showpoints) { points(c(Data$mod[[i]]$x), c(Data$mod[[i]]$y  ) , pch=21   , bg=pcol, fg='black' )  }
                cenP =  centroid(Data$mod[[i]])
                tolcol = tcolOFF
                if(!is.null(SELPOLY))
                  {
                    if(i %in% SELPOLY)
                      {
                        tolcol = tcolON
                        points(c(Data$mod[[i]]$x), c(Data$mod[[i]]$y  ) , pch=21   , bg=pcol, fg='black', cex=1.4 ) 
                      }
                  }
                rect(cenP[1]-TOL, cenP[2]-TOL, cenP[1]+TOL, cenP[2]+TOL, col=tolcol )
                text( cenP[1], cenP[2], paste(i,":", format(Data$mod[[i]]$rho)  , sep="") , col='black' )
              }
            if(!is.null(SELPOLY))
              {
                
                for(k in 1:length(SELPOLY))
                  {
                    i = SELPOLY[k]
                    points(c(Data$mod[[i]]$x), c(Data$mod[[i]]$y  ) , pch=21   , bg=acol, fg='black' )
                  }
                
              }
            
              }
            
        buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
        return(list(buttons=buttons))
      }
   ##################
   ##################
   ##################
   
    ZOLTAN = replot(Data)
   
    
    iloc = locator(1, type='p')
    zloc = iloc

    Nclick = length(iloc$x)
    if(is.null(zloc$x)) { return(NULL) }
    K =  RPMG::whichbutt(zloc , ZOLTAN$buttons)
    sloc = zloc
   
 
    while(TRUE)
      {
        ############   button actions

        ###########   quit and break loop
        if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
          {
            title("DONE", cex.main=3)
            break;
          }

        ###########   refresh the screen
        if(K[Nclick] == match("REFRESH", labs, nomatch = NOLAB))
          {
            zloc = list(x=NULL, y=NULL)
            ZOLTAN =replot(Data)
            buttons = ZOLTAN$buttons
          }


 if(K[Nclick] == match("SHOWpts", labs, nomatch = NOLAB))
          {
            showpoints = !showpoints
            zloc = list(x=NULL, y=NULL)
            ZOLTAN =replot(Data)
            buttons = ZOLTAN$buttons
          }

    if(K[Nclick] == match("COLOR", labs, nomatch = NOLAB))
          {
            showcolor = !showcolor
            zloc = list(x=NULL, y=NULL)
            ZOLTAN =replot(Data)
            buttons = ZOLTAN$buttons
          }

        

        ###########   Insert a polygon
        if(K[Nclick] == match("POLY", labs, nomatch = NOLAB))
          {
            npt = length(zloc$x) - 1
            Data$n = Data$n + 1
           #  pol = list(x=zloc$x[1:npt] , y = zloc$y[1:npt], cenx=NA, ceny=NA, rho=DefaultDen, active=1 )
            cenP =  centroid(list(x=zloc$x[1:npt] , y = zloc$y[1:npt])   )
            
            pol = list(x=zloc$x[1:npt] , y = zloc$y[1:npt], cenx=cenP[1], ceny=cenP[2], rho=DefaultDen, active=1 )
            cenP =  centroid(pol)
            if(Data$n==1)
              {
                Data$cens$x =  cenP[1]
                Data$cens$y = cenP[2]
              }
            else
              {
                Data$cens$x = c(Data$cens$x, cenP[1])
                Data$cens$y = c(Data$cens$y, cenP[2])
              }

            
            Data$mod[[Data$n]] = pol
            SELPOLY = Data$n
            zloc = list(x=NULL, y=NULL)
            ZOLTAN =replot(Data)
            buttons = ZOLTAN$buttons
          }
###########   establish a density contrast
        if(K[Nclick] == match("SELpol", labs, nomatch = NOLAB))
          {
          SELPOLY =  whichSELpoly(zloc, Data$cen, TOL)
           
        zloc = list(x=NULL, y=NULL)
        ZOLTAN =replot(Data)
        buttons = ZOLTAN$buttons
        
      }
        if(K[Nclick] == match("MVpt", labs, nomatch = NOLAB))
          {
            
            npt = length(zloc$x) - 1
            if(npt>1)
              {

                Ipol = SELPOLY[1]
                POL = Data$mod[[Ipol]]
                
                 i = npt-1
                  
                    x1 = zloc$x[i]
                    y1 = zloc$y[i]
                 sp1 = which( sqrt(  (POL$x-x1)^2 + (POL$y-y1)^2)  < TOL*.6 )
                    if(length(sp1)>0)
                      {
                      ##  points( POL$x[sp1],  POL$y[sp1], cex=2, col='blue')
                        POL$x[sp1] = zloc$x[npt]
                        POL$y[sp1] = zloc$y[npt]
                       ##  points( POL$x[sp1],  POL$y[sp1], cex=1, col='red')

                        cenP =  centroid(POL)
                       
                            Data$cens$x[Ipol] =  cenP[1]
                            Data$cens$y[Ipol] = cenP[2]
                          
                        Data$mod[[Ipol]]=POL


                        
                        
                      }
                  
                    
           
                
              }
        zloc = list(x=NULL, y=NULL)
        ZOLTAN =replot(Data)
        buttons = ZOLTAN$buttons
        
      }
        
        ###########   establish a density contrast
        if(K[Nclick] == match("Rho", labs, nomatch = NOLAB))
          {
            SELPOLY =  whichSELpoly(zloc, Data$cen, TOL)
            w1=SELPOLY[1]
            if(w1>0)
              {
                
                w1=SELPOLY[1]
                cat(paste("chose:", w1), sep="\n" )
                rho = readline(prompt = "Type in the density (cgm/cm^3) : ")
                myrho = as.numeric(rho)
                print(myrho)
                Data$mod[[w1]]$rho=myrho

                
              }
            
        
        zloc = list(x=NULL, y=NULL)
        ZOLTAN =replot(Data)
        buttons = ZOLTAN$buttons
      }
        ###########   establish a density contrast
        if(K[Nclick] == match("DELETE", labs, nomatch = NOLAB))
          {
            ####  delete a polygon
           
            w1=SELPOLY[1]

            
          Data = delete.poly(Data, w1)
        
        zloc = list(x=NULL, y=NULL)
        ZOLTAN =replot(Data)
        buttons = ZOLTAN$buttons
      }

        ###########   establish a density contrast
        if(K[Nclick] == match("MVpoly", labs, nomatch = NOLAB))
          {
            w1=SELPOLY[1]
            print(zloc)
            
            if(w1>=1)
              {
            ###  need 2 points
            npt = length(zloc$x)
            if(npt>2)
              {
                dvec = list(x=  zloc$x[2] -zloc$x[1]  , y=  zloc$y[2] -zloc$y[1]   )
                print(dvec)
                Data$mod[[w1]]$x = Data$mod[[w1]]$x+dvec$x
                Data$mod[[w1]]$y = Data$mod[[w1]]$y+dvec$y
                Data$cens$x[w1] = Data$cens$x[w1]+dvec$x
                Data$cens$y[w1] = Data$cens$y[w1]+dvec$y

                Data$mod[[w1]]$cenx = Data$cens$x[w1]
                Data$mod[[w1]]$ceny = Data$cens$y[w1]
              }

          }
            
        zloc = list(x=NULL, y=NULL)
        ZOLTAN =replot(Data)
        buttons = ZOLTAN$buttons
      }



##############  plot and calculate the gravity model
       if(K[Nclick] == match("MODEL", labs, nomatch = NOLAB))
          {
            IDcur = dev.cur()
            dev.new()
            BMOD(Data, nstn=100, obs=obs)
            dev.set(IDcur)
            zloc = list(x=NULL, y=NULL)
            ZOLTAN =replot(Data)
            buttons = ZOLTAN$buttons
          }
         if(K[Nclick] == match("SAVE", labs, nomatch = NOLAB))
          {
            save(file="Gsavemodel.Rdata", Data)
            zloc = list(x=NULL, y=NULL)
            ZOLTAN =replot(Data)
            buttons = ZOLTAN$buttons
          }
        iloc = locator(1,type='p')
##### print(iloc)
        zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
        Nclick = length(iloc$x)
        if(is.null(zloc$x)) { return(sloc) }
        K =  RPMG::whichbutt(iloc , ZOLTAN$buttons)
##### print(K)   
      }

   return(Data)
    
  }

