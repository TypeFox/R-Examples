`detail.pick` <-
function(y, ex, dt, TIT="")
  {
    if(missing(TIT)) { TIT=NULL }
   labs = c("DONE", "PROJ", "XING", "YMIN", "YMAX", "SAVED", "NONE" )
   colabs =  rep(1,length(labs))
   pchlabs = rep(1,length(labs))
   NSEL = 1
   N = 0
   NLABS = length(labs)
   NOLAB = NLABS +1000
   KSAVE = NULL
    xsave = NULL
    ysave = NULL
    pwink = 0.01*diff(range(ex))
    pcol=rgb(1,.5, 0)
    plot(ex, y, type='n', col=1)
    abline(h=0, col=1)
    points(ex, y, col=rgb(0.75,0.75,0.8) )
    lines(ex, y, col=1)
    
    title(main=TIT)
    
   ###    title(main="click middle mouse to end session")
    buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
    ###  zloc = plocator(col=rgb(1,0, 0), NUM=TRUE , YN=NSEL, style=-1)
   zloc = locator(1, type='p', col=pcol)

   
    
    Nclick = length(zloc$x)
    if(is.null(zloc$x)) { return(NULL) }
    K = RPMG::whichbutt(zloc ,buttons)
   sloc = zloc
   while(Nclick>0)
     {
        ##  abline(v=zloc$x[1], col=2)
        xsave = c(xsave, zloc$x)
        ysave = c(ysave, zloc$y)
        N = N+1
      ###  print(paste("CLICKER ", N, zloc$x, zloc$y , xsave[N], ysave[N], length(xsave)  ))
       if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
         {
           N = N-1
           xsave = xsave[1:N]
           ysave = ysave[1:N]

           break;
         }
       if(K[Nclick] == match("Postscript", labs, nomatch = NOLAB))
          {
          }


        
      if(K[Nclick] == match("XING", labs, nomatch = NOLAB))
          {
            ###  find the zero crossing

             ###print(paste("xing ", N, xsave[N], ysave[N]))
            ### print(paste("should be  ",xsave[N-1], ysave[N-1]))
             ### print(xsave)

            N = N-1
            xsave = xsave[1:N]
            ysave = ysave[1:N]

            
           LX = xsave[N]
           LY =  ysave[N]
           rim = findInterval(LX, ex)
           nflag = seq(from=(rim-5), to=rim+5, by=1)
           lex = ex[nflag]
           lwhy = y[nflag]
           sy = sign(lwhy[1])
           ww = which(sign(lwhy) !=sy) 
           
           x1 = lex[ww[1]-1]
           y1 =  lwhy[ww[1]-1]
           x2 = lex[ ww[1] ]
           y2 =  lwhy[ww[1]]

           m = (y2-y1)/(x2-x1)
           b = y2-m*x2

           xingx = -b/m
           xingy = 0
           points(c(x1,x2,xingx), c(y1, y2, xingy), col=2)

           
           xsave[N] = xingx
           ysave[N] = xingy
            text(xsave[N], ysave[N], labels= N, pos=3)
      

          }
   
      if(K[Nclick] == match("PROJ", labs, nomatch = NOLAB))
          {
            ###  find the zero crossing

             ###print(paste("xing ", N, xsave[N], ysave[N]))
            ### print(paste("should be  ",xsave[N-1], ysave[N-1]))
             ### print(xsave)

            n = length(xsave)
           ##  xsave = xsave[1:N]
            ##  ysave = ysave[1:N]

           x1 = xsave[n-1]
            
           y1 = ysave[n-1]
           x2 = xsave[n-2]
           y2 =  ysave[n-2]

           m = (y2-y1)/(x2-x1)
           b = y2-m*x2
          ###  print(paste(sep=' ',"slope and intercept", m,b))
          ###  print(paste(sep=' ', n, x1,y1, x2,y2))
           xingx = -b/m
           xingy = 0
           points(c(x1,x2,xingx), c(y1, y2, xingy), col=2)

            N = N-2
            xsave = xsave[1:N]
            ysave = ysave[1:N]
            
            xsave[N] = xingx
            ysave[N] = xingy
            text(xsave[N], ysave[N], labels= N, pos=3)
      

          }

        
      if(K[Nclick] == match("YMAX", labs, nomatch = NOLAB))
          {
            ##  find the local max


            N = N-1
            xsave = xsave[1:N]
            ysave = ysave[1:N]


            
           ### print(paste("MAX ", N, xsave[N], ysave[N]))
            LX = xsave[N]
            LY =  ysave[N]
            ax = LX
            
            flag = ex > (ax-pwink) & ex < (ax+pwink)
            w1 = which(flag)[1]-1
            
            rim = which.max(y[flag])
            abline(v=ex[w1+rim], col=4)

           
           xsave[N] = ex[w1+rim]
           ysave[N] = y[w1+rim]
             text(xsave[N], ysave[N], labels= N, pos=3)
 

            points(xsave[N],  ysave[N] , col=3, pch=7) 
        
           ###  points(xsave[N],  ysave[N] , col=3, pch=7) 
           
            
          }
       if(K[Nclick] == match("YMIN", labs, nomatch = NOLAB))
          {


             N = N-1
            xsave = xsave[1:N]
            ysave = ysave[1:N]

           ### print(paste("MIN ", N, xsave[N], ysave[N]))
            LX = xsave[N]
            LY =  ysave[N]
            ax = LX
            flag = ex > (ax-pwink) & ex < (ax+pwink)
            w1 = which(flag)[1]-1

             points( ex[flag]  ,  y[flag], col=5, pch=7) 
            rim = which.min(y[flag])

            
            abline(v=ex[w1+rim], col=4)
            
                  
            xsave[N] = ex[w1+rim]
            ysave[N] = y[w1+rim]
             text(xsave[N], ysave[N], labels= N, pos=3)
            
            points(xsave[N],  ysave[N] , col=3, pch=7) 
           
            
          }
         if(K[Nclick] == match("NONE", labs, nomatch = NOLAB))
          {
            N = 0
            KSAVE = NULL
            xsave = NULL
            ysave = NULL
            plot(ex, y, type='n', col=1)
            abline(h=0, col=1)
             points(ex, y, col=rgb(0.75,0.75,0.8) )
            lines(ex, y, col=1)
            title(main=TIT)
###    title(main="click middle mouse to end session")
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
          }
        if(K[Nclick] == match("POINTS", labs, nomatch = NOLAB))
          {
            if(N>1)
              {
                N = N-1
                xsave = xsave[1:N]
                ysave = ysave[1:N]
              }
            else
              {
                N = 0
                xsave = NULL
                ysave = NULL
              }
            points(ex, y, col=rgb(0.8,0.8,0.8) )
            
          }
       if(K[Nclick] == match("SAVED", labs, nomatch = NOLAB))
          {
            
            points(xsave, ysave, col=rgb(0.5,1, 0.5) )
            text(xsave, ysave, labels=1:length(xsave), pos=1)
          }
      
      
       
       #### buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
       ####  zloc = plocator(col=rgb(1,0, 0), NUM=TRUE , YN=NSEL, style=-1)
        zloc = locator(1, type='p', col=pcol)
       Nclick = length(zloc$x)
       if(is.null(zloc$x)) { return(sloc) }
       K =  RPMG::whichbutt(zloc ,buttons)
      ###  print(paste("Button",K, Nclick))
     }
   
   

   KSAVE = list(x=xsave, y=ysave)
   
   return(KSAVE)
 }

