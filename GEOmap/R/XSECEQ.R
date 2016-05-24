XSECEQ<-function(MAP, EQ , XSECS=NULL, labs=c("DONE","REFRESH", "XSEC", "MSEC"), width=10, kmaxes=TRUE , pch=".", demo=FALSE)
  {
    if(missing(demo)) { demo = FALSE }
    if(missing(width)) { width=10 }
    if(missing(pch)) { pch="." }
    if(missing(XSECS)) { XSECS=NULL}
    if(missing(kmaxes)) { kmaxes=TRUE}
    if(missing(labs)) { labs=c("DONE","REFRESH", "XSEC", "MSEC",  "KLEAN", "KMAXES", "COLZ", "CONT", "width", "PS" ) }


    polycol = 'blue'
    linecol='black'
    ZCOLOR.FLAG = FALSE
#######  
    
##### source("/home/lees/XSECEQ.R")

    ###  labs=c("DONE","REFRESH", "XSEC", "MSEC",  "KMAXES", "CONT", "width", "DIST", "PS" )
    
  ### XSECEQ(  MAP, EQ , labs, demo=FALSE  )   
    #####  
  
      
    iseclab  = 0
    secmat = NULL
    ncol = 100
    TPALS = c("rainbow", "topo", "terrain")
    colabs = rep(1, length=length(labs))
    pchlabs = rep(0,length(labs))
    FUN = match.fun(TPALS[1])
    pal = FUN(ncol)

    WINWIDTH = 8
    CONT.FLAG = FALSE
    XSEC.FLAG = FALSE
    PS.FLAG =  FALSE
    

    XYLIM = attr(MAP, "XYLIM")
    PROJ  = attr(MAP, "PROJ")
    

    GRIDcol=1


    ur = par("usr")
    urx = expandbound(ur[1:2], -.05)
    ury = expandbound(ur[3:4], -.05)

    
    
    Abound = XY.GLOB(urx , ury, PROJ)
    
    PLAT =  pretty( Abound$lat)
    PLAT = c(min(Abound$lat),  PLAT[PLAT>min(Abound$lat) & PLAT<max(Abound$lat)],max(Abound$lat))
    PLON  = pretty(Abound$lon)

    if(kmaxes)
      {
        kmaxis = TRUE
      }
    else
      {
        
        kmaxis = FALSE
      }
    

    
    ##########  map function:
    plot(MAP[[1]]$x, MAP[[1]]$y, xlim=XYLIM$x, ylim= XYLIM$y,   type='n', asp=1, ann=FALSE, axes=FALSE )


    if(kmaxis==TRUE)
      {
        axis(1)
        axis(2)
        title(xlab="km", ylab="km")
        box()
        
      }
    else
      {
       # addLLXY(PLAT, PLON,  GRIDcol=GRIDcol, LABS=0, BORDER=0 , PROJ=PROJ )
        print("sqrtix")
        sqrTICXY(XYLIM, proj=PROJ)
        
      }
    
    for(i in 2:length(MAP)) {
      
      plotGEOmapXY(MAP[[i]], LIM= c(XYLIM$lon[1], XYLIM$lat[1], XYLIM$lon[2], XYLIM$lat[2])  , PROJ=PROJ, add=TRUE, shiftlon=0)
      
    }
    for(i in 1:length(EQ)) {
      
      if(ZCOLOR.FLAG)
        {
          points(EQ[[i]]$x, EQ[[i]]$y, col=EQ[[i]]$COL, pch=EQ[[i]]$pch, cex=EQ[[i]]$cex)
        }
      else
        {
          points(EQ[[i]]$x, EQ[[i]]$y, col=EQ[[i]]$col, pch=EQ[[i]]$pch, cex=EQ[[i]]$cex)
        }
      
     
      
    }
##########
   
   SW = list()
    
  if(!is.null(XSECS))
    {
      XSEC.FLAG = TRUE
      dev1 = dev.cur()
      
      for(i in 1:length(XSECS))
        {
          iseclab =  iseclab +1
          SW[[iseclab]] = XSECS[[i]]
          secmat = rbind(secmat, attr(SW[[iseclab]], "xsec" ))
          polygon(SW[[iseclab]]$InvBox, border=polycol)
          seclabs = attr(SW[[iseclab]],"LAB")
          segments(secmat[iseclab,1],secmat[iseclab,2],secmat[iseclab,3],secmat[iseclab,4], col=linecol)

          ang = (180/pi)*atan2(secmat[iseclab,4]-secmat[iseclab,2], secmat[iseclab,3]-  secmat[iseclab,1] )
          
          text(secmat[iseclab,1],secmat[iseclab,2], labels= seclabs, pos=3, srt=ang)
          text(secmat[iseclab,3],secmat[iseclab,4], labels= paste(sep="",seclabs, "'") , pos=3, srt=ang)

          
          
        }
      for(i in 1:length(SW))
        {
          rr = range(SW[[i]]$r)
          rdep = range(SW[[i]]$depth)
          dr = abs(diff(rr))
          ddep = abs(diff(rdep))
          ##  rat = diff(rr)/diff(rdep)

          if(dr>0 & ddep>0)
            {
              wid = WINWIDTH
              hei = wid*ddep/dr
             ##  print(paste(i, wid, hei, dr, ddep))
              dev.new(width=wid, height=hei)
               seclabs = attr(SW[[i]],"LAB")
              XSECwin( SW[[i]] , i, xLAB = seclabs,   demo=TRUE  )
            }
        }
      dev.set(dev1)

      
    }

    if(demo==TRUE)  return(NULL)

    
    cdev = dev.cur()
    
    buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
    NLABS = length(labs)
    NOLAB = NLABS +1000  ## some large number
    
    
    iloc = locator(1, type='p')
    zloc = iloc
    
    Nclick = length(iloc$x)
    if(is.null(zloc$x)) { return(NULL) }
    K =  RPMG::whichbutt(zloc , buttons)
    sloc = zloc
    


    while(TRUE)
      {
        ############   button actions

        ###########   quit and break loop
        if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
          {


             buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.8), length(labs)), pch=rep("NULL", length(labs)))
            title("Return to Calling Program")
            
            break;
          }
        if(K[Nclick] == match("DIST", labs, nomatch = NOLAB))
          {

            zen = length(zloc$x)
            if(zen > 3)
              {
            LEN = sqrt( (zloc$x[zen-1] - zloc$x[zen-2])^2+(zloc$y[zen-1] - zloc$y[zen-2])^2)
            print(paste("Length=", LEN, "km"))
          }
            
          }

        ###########   refresh the screen
        if(K[Nclick] == match("REFRESH", labs, nomatch = NOLAB))
          {
            zloc = list(x=NULL, y=NULL)
          }

        ###########   refresh the screen
        if(K[Nclick] == match("Next", labs, nomatch = NOLAB))
          {
            dev.set(dev.next())
            
            zloc = list(x=NULL, y=NULL)
          }
       
        if(K[Nclick] == match("COLZ", labs, nomatch = NOLAB))
          {
            ZCOLOR.FLAG = !ZCOLOR.FLAG
            zloc = list(x=NULL, y=NULL)
          }
################  add contours
        if(K[Nclick] == match("CONT", labs, nomatch = NOLAB))
          {
            CONT.FLAG = !CONT.FLAG
            zloc = list(x=NULL, y=NULL)
          }
 ################  
        if(K[Nclick] == match("KMAXES", labs, nomatch = NOLAB))
          {
            kmaxis = !kmaxis
            zloc = list(x=NULL, y=NULL)
          }

 ################  
       



        
        if(K[Nclick] == match("KLEAN", labs, nomatch = NOLAB))
          {
           kdev = dev.cur()
           devall = dev.list()
           
           for(i in 1:length(devall))
             {
               if(devall[i]==kdev) next()
               dev.off(which = devall[i])
               
             }
           SW = list()
           secmat = NULL
           XSEC.FLAG = FALSE
           iseclab = 0
            zloc = list(x=NULL, y=NULL)
          }

        ################  make postscript file
        if(K[Nclick] == match("PS", labs, nomatch = NOLAB))
          {
            PS.FLAG = !PS.FLAG
            zloc = list(x=NULL, y=NULL)
          }
        if(K[Nclick] == match("width", labs, nomatch = NOLAB))
          {
            print(paste(sep=" ", "Current width=", width))
            answer = readline(prompt = "type in the width (km): ")

            width = as.numeric(answer)
            
            if(!is.numeric(width)) width=10

            
            if(width<=0) width=10

            print(paste(sep=" ", "New width=", width))
            
            zloc = list(x=NULL, y=NULL)
          }

################  cut cross sections
        if(K[Nclick] == match("XSEC", labs, nomatch = NOLAB))
          {
            n = length(zloc$x)
            x1 = zloc$x[n-2]
            y1 = zloc$y[n-2]
            x2 = zloc$x[n-1]
            y2 = zloc$y[n-1]

            XSEC.FLAG = TRUE
            print(c(x1,y1, x2,y2))
            iseclab = iseclab + 1
            LAB = LETTERS[iseclab]
            secmat = rbind(secmat, c(x1, y1, x2,  y2))

            
           ##### GETXprofile(jx, jy, Data, myloc=list(x=c(x1, x2), y=c(y1, y2)), LAB=LAB, PLOT=TRUE)
            L = list(x=c(x1, x2) ,  y=c(y1, y2))
            
            j = 1
            SW[[iseclab]] = list()
             SW[[iseclab]] = eqswath(EQ[[j]]$x, EQ[[j]]$y, EQ[[j]]$z, L, width=width, PROJ=PROJ)

            polygon(SW$InvBox, border="blue")

            attr(SW[[iseclab]], "LAB" ) <- LAB
            attr(SW[[iseclab]], "NUM" ) <- iseclab
            attr(SW[[iseclab]], "xsec" ) <- c(x1,y1,x2, y2  )
            
            SW[[iseclab]]$proj = PROJ
            
           #### get(getOption("device"))()


            rr = range(SW[[iseclab]]$r)
            rdep = range(SW[[iseclab]]$depth)
            dr = abs(diff(rr))
            ddep = abs(diff(rdep))
            
            if(dr>0 & ddep>0)
              {
                wid = WINWIDTH
                hei = wid*ddep/dr
                if(hei<4) hei=4
              }
            else
              {
                wid = WINWIDTH
                hei = 4
              }
            
            dev.new(width=wid, height=hei)
            
            xlabs=c("DONE","REFRESH", "PS" )
            XSECwin( SW[[iseclab]] , iseclab, LAB , xlabs, demo=FALSE  )   
            
####   plot(SW[[iseclab]]$r , -SW[[iseclab]]$depth,  main=paste( iseclab, LAB) , xlab="km", ylab="Depth", asp=1)

            
        
            dev.set(cdev)

            
            zloc = list(x=NULL, y=NULL)

           

            
          }
################  cut multiple cross sections
        if(K[Nclick] == match("MSEC", labs, nomatch = NOLAB))
          {

            XSEC.FLAG = TRUE

            
            n = length(zloc$x)

            for(M in seq(from=1, to=n-2, by=2))
              {
                x1 = zloc$x[M]
                y1 = zloc$y[M]
                x2 = zloc$x[M+1]
                y2 = zloc$y[M+1]
                
                
                print(c(x1,y1, x2,y2))
                iseclab = iseclab + 1
                LAB = LETTERS[iseclab]
                secmat = rbind(secmat, c(x1, y1, x2,  y2))
                
                
##### GETXprofile(jx, jy, Data, myloc=list(x=c(x1, x2), y=c(y1, y2)), LAB=LAB, PLOT=TRUE)
                L = list(x=c(x1, x2) ,  y=c(y1, y2))
                
                j = 1
                SW[[iseclab]] = list()
                SW[[iseclab]] = eqswath(EQ[[j]]$x, EQ[[j]]$y, EQ[[j]]$z, L, width=width, PROJ=PROJ)
                
                polygon(SW$InvBox, border="blue")
                
                attr(SW[[iseclab]], "LAB" ) <- LAB
                attr(SW[[iseclab]], "NUM" ) <- iseclab
                attr(SW[[iseclab]], "xsec" ) <- c(x1,y1,x2, y2  )

                SW[[iseclab]]$proj = PROJ
                
              ####   get(getOption("device"))()



                rr = range(SW[[iseclab]]$r)
                rdep = range(SW[[iseclab]]$depth)
                dr = abs(diff(rr))
                ddep = abs(diff(rdep))
                
                if(dr>0 & ddep>0)
                  {
                    wid = WINWIDTH
                    hei = wid*ddep/dr
                  }
                else
                  {
                    wid = WINWIDTH
                    hei = 4
                  }

    
                  if(hei<4)  hei = 4

                
                
                dev.new(width=wid, height=hei)
                
                
                xlabs=c("DONE","REFRESH", "PS" )
                XSECwin( SW[[iseclab]] , iseclab, LAB , xlabs, demo=TRUE  )   
                
          ####   plot(SW[[iseclab]]$r , -SW[[iseclab]]$depth,  main=paste( iseclab, LAB) , xlab="km", ylab="Depth", asp=1)
                
                
                
                dev.set(cdev)
              }

            
            zloc = list(x=NULL, y=NULL)
            
            
            
            
          }

        if(K[Nclick] > 0)
          {
            if(PS.FLAG) {
              P = round(par('pin'), digits=2);
              psname = RPMG::local.file("XMAP", "eps")

              postscript(file=psname  , width=P[1], height=P[2],
                         paper = "special", horizontal=FALSE, onefile=TRUE,print.it=FALSE)
            }


##########  map function:
            plot(MAP[[1]]$x, MAP[[1]]$y, xlim=XYLIM$x, ylim= XYLIM$y,   type='n', asp=1, ann=FALSE, axes=FALSE )
            if(kmaxis==TRUE)
              {
                axis(1)
                axis(2)
                title(xlab="km", ylab="km")
                box()
                
              }
            else
              {
                u = par("usr")
                
                ## print("sqrtix")


                ##  XYLIM
                sqrTICXY(list(x=u[1:2], y=u[3:4] ) , proj=PROJ)
        
                ##  addLLXY(PLAT, PLON,  GRIDcol=GRIDcol, LABS=TRUE, BORDER=0 , PROJ=PROJ )
                
              }
            
            for(i in 2:length(MAP)) {
              
              plotGEOmapXY(MAP[[i]], LIM= c(XYLIM$lon[1], XYLIM$lat[1], XYLIM$lon[2], XYLIM$lat[2])  , PROJ=PROJ, add=TRUE, shiftlon=0)
              
            }


    ##########

            
            for(i in 1:length(EQ)) {

              if(ZCOLOR.FLAG)
                {
                  points(EQ[[i]]$x, EQ[[i]]$y, col=EQ[[i]]$COL, pch=EQ[[i]]$pch, cex=EQ[[i]]$cex)
                }
              else
                {
                  points(EQ[[i]]$x, EQ[[i]]$y, col=EQ[[i]]$col, pch=EQ[[i]]$pch, cex=EQ[[i]]$cex)
                }
              
            }
    ##########

          
            
            #####    if(CONT.FLAG) contour(x=jx, y=jy, Data, add=TRUE)
            if(XSEC.FLAG) 
              {

                
                segments(secmat[,1],secmat[,2],secmat[,3],secmat[,4])
                seclabs = LETTERS[1:iseclab]
                
                for(k in 1:length(SW))
                  {
                    polygon(SW[[k]]$InvBox, border=polycol)
                   
                    
                    seclabs = attr(SW[[k]],"LAB")
                    segments(secmat[k,1],secmat[k,2],secmat[k,3],secmat[k,4], col=linecol)
                    
                    ang = (180/pi)*atan2(secmat[k,4]-secmat[k,2], secmat[k,3]-  secmat[k,1] )
                    
                    text(secmat[k,1],secmat[k,2], labels= seclabs, pos=3, srt=ang)
                    text(secmat[k,3],secmat[k,4], labels= paste(sep="",seclabs, "'") , pos=3, srt=ang)
          
                   
                  }

                
               

               

              }

            
           
            
            if(PS.FLAG) {
              dev.off();
              cat(paste(sep=" ", "the postscript file is in: ", psname), sep="\n")
              PS.FLAG =  FALSE
            }

            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)


          }
        else
          {
###  in case the plot was resized with asp=1, need to replot the buttons
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
          }


        iloc = locator(1,type='p')
##### print(iloc)
        zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
        Nclick = length(iloc$x)
        if(is.null(zloc$x)) { return(sloc) }
        K =  RPMG::whichbutt(iloc , buttons)
##### print(K)   
      }

    invisible(SW)
    
  }

