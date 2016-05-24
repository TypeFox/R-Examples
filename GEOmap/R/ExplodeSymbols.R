ExplodeSymbols<-function(XY, fsiz=1, STARTXY=NULL, MAP=NULL )
  {
####  find a good way to explode symbols (focal mechanisms)
    ####  source("ExplodeSymbols.R")
    
    
    if(missing(MAP)) { MAP = NULL }
    if(missing(STARTXY))
      {
        STARTXY=NULL 
      }

    if(is.null(STARTXY))
      {
        NEWXY = XY
        STARTXY   = XY
      }
    else
      {
        NEWXY =STARTXY

      }


    
    draw.circ<-function (x, y, r, ...) 
      {
        CI = RPMG::circle(1)
        for (i in 1:length(x)) {
          Cx = x[i] + r * CI$x
          Cy = y[i] + r * CI$y
          lines(c(Cx, Cx[1]), c(Cy, Cy[1]), type = "l", ...)
        }
      }


    EXPreplot<-function()
      {
        rx = range(c(NEWXY$x, XY$x))
        ry = range(c(NEWXY$y, XY$y))

        drx = diff(rx)
        dry = diff(ry)
        
        rx = c(rx[1]-XPCT*drx, rx[2]+XPCT*drx)
        ry = c(ry[1]-XPCT*dry, ry[2]+XPCT*dry)
        
        plot(rx , ry , type='n', asp=1, xlab="km", ylab="km")

        if(!is.null(MAP) ) { plotGEOmap(MAP, LIM = c(rx[1],ry[1] , rx[2], ry[2]) , add = TRUE ) }

        points(XY$x, XY$y, pch=3, col='blue')
        points(NEWXY$x, NEWXY$y, pch=6, col='red')
        segments(XY$x, XY$y, NEWXY$x, NEWXY$y)
        u = par("usr")
        
        focsiz = fsiz* (u[2]-u[1])
        draw.circ(NEWXY$x, NEWXY$y, focsiz)

        invisible(focsiz)
        
      }

   
    XPCT = 0
    focsiz =EXPreplot()

    labs = c("DONE", "QUIT", "SEL","CIRC", "RECT2","LINE","NOVER", "HAND", "RECT",  "Xpand", "Zoom", "REFRESH", "RESTORE", "DOC" ) 
    colabs = rep("blue", length = length(labs))
    colabs[labs=="LINE"] = "green"
    colabs[labs=="RECT2"] = "red"
    colabs[labs=="CIRC"] = "red"
    
    pchlabs = rep(0, length(labs))
    buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
    NLABS = length(labs)
    NOLAB = NLABS + 1000


    SEL = NULL

    iloc = locator(1, type = "p")
    zloc = iloc
    Nclick = length(iloc$x)
    if (is.null(zloc$x)) {
      return(NULL)
    }
    K = RPMG::whichbutt(zloc, buttons)
    sloc = zloc
    while (TRUE) {

##############################################################      
      if (K[Nclick] == match("DONE", labs, nomatch = NOLAB)) {
        break
      }
##############################################################      
      if (K[Nclick] == match("QUIT", labs, nomatch = NOLAB)) {
        ##  print("Pressed Quit")
         break
        zloc = list(x = NULL, y = NULL)
      }
##############################################################      
      if (K[Nclick] == match("HAND", labs, nomatch = NOLAB)) {
        print("Pressed Hand")
        NP = length(zloc$x)-1
        if(NP<2)
          {
            print("Must have at least 2 clicks on screen to move a symbol")
            focsiz = EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
            print(zloc)
            zloc = list(x = NULL, y = NULL)
            next

          }
        kloc = list(x=zloc$x[1:NP], y=zloc$y[1:NP ])
        print(kloc)
        for(k in seq(from=1, to = length(kloc$x), by=2) )
          {
            kpoint = which.min(  (kloc$x[k]-NEWXY$x)^2  +  (kloc$y[k]-NEWXY$y)^2)
            
            NEWXY$x[kpoint] = kloc$x[k+1]
            NEWXY$y[kpoint] = kloc$y[k+1]

          }

        print(NEWXY)
        focsiz =EXPreplot()
        buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
        print(zloc)
        zloc = list(x = NULL, y = NULL)
      }

##############################################################
      if (K[Nclick] == match("SEL", labs, nomatch = NOLAB)) {
####  select points in rect
        NP = length(zloc$x)-1
        if(NP<2)
          {
            print("Must have at least 2 clicks on screen to select a set of symbols")
            focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
            print(zloc)
            zloc = list(x = NULL, y = NULL)
            next
            
          }
        else
          {

###  take last two clicks before clicking the RECT button
            p1 = list(x=c(zloc$x[(NP-1):NP]), y=c(zloc$y[(NP-1):NP]))
            
            SELX = which(NEWXY$x>min(p1$x) & NEWXY$x<max(p1$x) & NEWXY$y>min(p1$y) & NEWXY$y<max(p1$y))

            SEL = SELX

            rect(min(p1$x), min(p1$y), max(p1$x), max(p1$y), border="red", lty=2)
            
            cat("Selection is done, choose a button to apply to selection.", sep="\n")

           ##  zloc = list(x = NULL, y = NULL)

          }
      }
 ##############################################################

     if (K[Nclick] == match("NOVER", labs, nomatch = NOLAB)) {
       

        NP = length(zloc$x)-1
        if(NP<2)
          {
            print("Must have at least 2 clicks on screen to select a set of symbols")
            focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
            print(zloc)
            zloc = list(x = NULL, y = NULL)
            next
            
          }
        else
          {

  
            bx = NEWXY$x
            by = NEWXY$y

          #######NXY = NoOverlap(x,y,focsiz, cenx=PL$x, ceny=PL$y)
            NXY = NoOverlap(bx,by,focsiz, SEL=SEL,  OLDx=XY$x, OLDy=XY$y)
            
            NEWXY$x[SEL] =NXY$x[SEL]
            NEWXY$y[SEL] =NXY$y[SEL]
          }

        focsiz =EXPreplot()
        buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
        print(zloc)
        zloc = list(x = NULL, y = NULL)
      }
      

      
      if (K[Nclick] == match("LINE", labs, nomatch = NOLAB)) {
        print("Pressed Line")

        NP = length(zloc$x)-1
        if(NP<2)
          {
            print("Must have at least 2 clicks on screen to select a set of symbols")
            focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
            print(zloc)
            zloc = list(x = NULL, y = NULL)
            next
            
          }
        else
          {

  
            bx = NEWXY$x[SEL]
            by = NEWXY$y[SEL]
           ##   bx = runif(10, -10, 10)
           ##  by = runif(10, -10, 10)

          ##   plot(c(-18, 18) , c(-18, 18), type='n', asp=1)
          ##   points(bx, by)

           ##    zloc = locator(4, type='p', col='red', pch=3)
            
          ##   NP = length(zloc$x)-1 
            p1 = list(x=c(zloc$x[(NP-2):NP]), y=c(zloc$y[(NP-2):NP]))


          ##   segments(p1$x[1], p1$y[1], p1$x[2], p1$y[2], col='red', lty=2, lwd=2)
            
            
            
          
            mfxy = list(x=p1$x[1],y=p1$y[1])

            ##   vector from P1 to P2
            V1 = list(x=p1$x[2]-p1$x[1], y=p1$y[2]-p1$y[1])
            
            V3 = list(x=p1$x[3]-p1$x[1], y=p1$y[3]-p1$y[1])
            absV1= sqrt(V1$x^2+V1$y^2)
            absV3= sqrt(V3$x^2+V3$y^2)
            
            a3  = sqrt(V3$x^2 + V3$y^2)
            c3 =  (V1$x*V3$x + V1$y*V3$y)/absV1

            V1vec = list(x=V1$x/absV1, y=V1$y/absV1)

            V3vec = list(x=V3$x/absV3, y=V3$y/absV3)
            
         ##    arrows(p1$x[1],p1$y[1], p1$x[1]+absV1*V1vec$x,p1$y[1]+absV1*V1vec$y)
         ##    arrows(p1$x[1],p1$y[1], p1$x[1]+absV3*V3vec$x,p1$y[1]+absV3*V3vec$y)

         ##    arrows(p1$x[1],p1$y[1], p1$x[1]+c3*V1vec$x,p1$y[1]+c3*V1vec$y, col='blue' )
            
            

            width1 = sqrt(  (p1$x[3]- (p1$x[1]+c3*V1vec$x))^2 + (p1$y[3]- (p1$y[1]+c3*V1vec$y))^2)

            perpvec = list(x=-V1vec$y, y=V1vec$x)



         ##    arrows(p1$x[1]+c3*V1vec$x,p1$y[1]+c3*V1vec$y    ,
## p1$x[1]+c3*V1vec$x+width1*perpvec$x,    p1$y[1]+c3*V1vec$y+width1*perpvec$y , col='purple'  )

            

       ##      cornsX = c(p1$x[1]+width1*perpvec$x, p1$x[2]+width1*perpvec$x, p1$x[2]-width1*perpvec$x,p1$x[1]-width1*perpvec$x,p1$x[1]+width1*perpvec$x )

        ##     cornsY = c(p1$y[1]+width1*perpvec$y, p1$y[2]+width1*perpvec$y, p1$y[2]-width1*perpvec$y,p1$y[1]-width1*perpvec$y,p1$y[1]+width1*perpvec$y )

        ##     lines(cornsX, cornsY)

            #############


B1 = list(x=bx-p1$x[1], y=by-p1$y[1])

            
            B3 =  (V1$x*B1$x + V1$y*B1$y)/absV1

## points(p1$x[1]+B3*V1vec$x,p1$y[1]+B3*V1vec$y)

##   plot(c(-18, 18) , c(-18, 18), type='n', asp=1)
##             points(bx, by)
## segments(p1$x[1], p1$y[1], p1$x[2], p1$y[2], col='red', lty=2, lwd=2)
            
##             arrows(p1$x[1]+B3*V1vec$x, p1$y[1]+B3*V1vec$y   ,  bx,by, col="brown")


            

widthB  = sqrt(  (bx- (p1$x[1]+B3*V1vec$x))^2 + (by- (p1$y[1]+B3*V1vec$y))^2)            
  
            g1vec = list(x=(bx-(p1$x[1]+B3*V1vec$x))/widthB , y=(by-(p1$y[1]+B3*V1vec$y))/widthB)
            
    ##            arrows(p1$x[1]+B3*V1vec$x,    p1$y[1]+B3*V1vec$y   ,           p1$x[1]+B3*V1vec$x+ widthB*g1vec$x , p1$y[1]+B3*V1vec$y  + widthB*g1vec$y   , col=grey(.8) )

           
##    arrows(bx, by, bx+ (width1-widthB)*g1vec$x   , by+  (width1-widthB)*g1vec$y)

            
       ##       NEWXY$x[SEL]=bx+ (width1-widthB)*g1vec$x
       ##       NEWXY$y[SEL]=by+  (width1-widthB)*g1vec$y)
          
            NEWXY$x[SEL]=bx+ (width1)*g1vec$x
            NEWXY$y[SEL]=by+  (width1)*g1vec$y
  
          }
        
         focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
        print(zloc)
        zloc = list(x = NULL, y = NULL)
      }

##############################################################
      if (K[Nclick] == match("Xpand", labs, nomatch = NOLAB)) {

        
        XPCT = XPCT + 0.1
        
         focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
       
        zloc = list(x = NULL, y = NULL)
      }
##############################################################    
       if (K[Nclick] == match("Zoom", labs, nomatch = NOLAB)) {

         XPCT = XPCT - 0.1

         if(XPCT<0) { XPCT = 0 }

         
         
         focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
       
        zloc = list(x = NULL, y = NULL)
      }
     
##############################################################      
      if (K[Nclick] == match("REFRESH", labs, nomatch = NOLAB)) {
       
         focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
       
        zloc = list(x = NULL, y = NULL)
      }
##############################################################      
      if (K[Nclick] == match("RESTORE", labs, nomatch = NOLAB)) {

        NEWXY = STARTXY
         focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
       
        zloc = list(x = NULL, y = NULL)
      }
##############################################################
      if (K[Nclick] == match("CIRC", labs, nomatch = NOLAB)) {


####  select points in rect

        NP = length(zloc$x)-1
        if(NP<2)
          {
            print("Must have at least 2 clicks on screen to select center and radius of circle")
            focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
            print(zloc)
            zloc = list(x = NULL, y = NULL)
            next
            
          }
        else
          {
            p1 = list(x=c(zloc$x[(NP-1):NP]), y=c(zloc$y[(NP-1):NP]))
            
            RAD = sqrt(  (p1$x[1]-p1$x[2])^2 +(p1$y[1]-p1$y[2])^2)

            bx = NEWXY$x[SEL]
            by = NEWXY$y[SEL]
            mfxy = list(x=p1$x[1],y=p1$y[1])
            
            dis1 = sqrt( (bx-mfxy$x)^2 + (by-mfxy$y)^2)
            
            
            dixplo  = RAD
            DX = (bx-mfxy$x)/dis1
            DY = (by-mfxy$y)/dis1
            
            PX = p1$x[1]+ dixplo* DX
            PY = p1$y[1]+ dixplo* DY
            
            NEWXY$x[SEL] = PX
            NEWXY$y[SEL] = PY 
            
          }
        focsiz =EXPreplot()
        buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
        zloc = list(x = NULL, y = NULL)
      }

##############################################################
      if (K[Nclick] == match("RECT2", labs, nomatch = NOLAB)) {


####  select points in rect

        NP = length(zloc$x)-1
        if(NP<2)
          {
            print("Must have at least 2 clicks on screen to select center and size of box")
            focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
            print(zloc)
            zloc = list(x = NULL, y = NULL)
            next
            
          }
        else
          {
            ##########  zloc = locator(2); NP = 2
            
            p1 = list(x=c(zloc$x[(NP-1):NP]), y=c(zloc$y[(NP-1):NP]))

            ###  rect(p1$x[1], p1$y[1], p1$x[2], p1$y[2])

            
            
            diagon = sqrt(  (p1$x[1]-p1$x[2])^2 +(p1$y[1]-p1$y[2])^2)


            theta = atan2((p1$y[1]-p1$y[2]), (p1$x[1]-p1$x[2]))

            digx = diagon*cos(theta)
            digy = diagon*sin(theta)

#### plot(epos$x, epos$y)
#### points(NEWXY$x[SEL] , NEWXY$y[SEL] , col='red')
#### bx = NEWXY$x[SEL]
#### by = NEWXY$y[SEL]

            
            rect(p1$x[1]-digx, p1$y[1]-digy, p1$x[1]+digx, p1$y[1]+digy)

            

            bx = NEWXY$x[SEL]
            by = NEWXY$y[SEL]


            rekt = list(x=c(p1$x[1]-digx, p1$x[1]+digx, p1$x[1]+digx, p1$x[1]-digx),
              y=c(p1$y[1]-digy, p1$y[1]-digy, p1$y[1]+digy,p1$y[1]+digy ))

            pp = list(x1=rep(p1$x[1], length(bx))  , y1=rep(p1$y[1], length(bx))    , x2=bx  , y2=by  )

            edges = rekt2line(rekt, pp )


          ####   points(edges, pch=6, col='brown')
          ####   segments(bx, by, edges$x, edges$y)

            

           
            NEWXY$x[SEL] = edges$x
            NEWXY$y[SEL] = edges$y
            
          }
        focsiz =EXPreplot()
        buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
        zloc = list(x = NULL, y = NULL)
      }
      


      
###  source("ExplodeSymbols.R")
###  NEWLOCS = ExplodeSymbols(SMXY, 0.03)
##############################################################
      if (K[Nclick] == match("DOC", labs, nomatch = NOLAB)) {


        exdoc = c(
          "############################",
          "############################",
          "This is an interactive program.  You must click on the screen",
          "and on the buttons to get this working -",
          "  the program will not work in batch mode or run as a script",
          "You click in the active screen area and then press a button",
          "on top (or bottom) - the button takes your clicks and does something",
          "Here are some hints:","\n",
          "HAND:",
          "If you want to move only one symbol (focal mech) click near it and then click",
          "where you want it to go.  Then click the HAND button",
          "You may click several at once...but for each click oin a symbol",
          "there has to be a click somewhere to relocate it.",
          "(i.e. there must be an even number of clicks on the screen before hitting the",
          " HAND button)","\n",

          "SEL:",
          "If you want to explode several symbols at once, first select them:",
          "  click lower left, then upper right of rectangle enclosing the selection.",
          "Once a selection is made it remains active until another selection is made",
          "so you can keep changing the radius and center for different explosions",
          "Then click CIRC.","\n",

          "RECT:",
          "Choose a rectangle (lower left and upper right), then click RECT for an explosion",
          "\n",
          "RECT2:","After selecting, choose a center and a distance.",
          "symbols will be moved to a rectangular perimeter defined by the two points",
          "\n",

          "CIRC:",
          "After selection, click once for the circle center, and a second time for the radius, then click CIRC",
          "\n",
          "LINE:\n","After selection, three clicks: 1,2 for line, 3 for distance from line",
          "############################",
          "############################",
          "\n",
          "\n","\n"

          )


        for(kdoc in 1:length(exdoc)){ cat(exdoc[kdoc], sep="\n") }


        
        focsiz =EXPreplot()
        buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
        zloc = list(x = NULL, y = NULL)
      }
##############################################################      
      if (K[Nclick] == match("RECT", labs, nomatch = NOLAB)) {


####  select points in rect

        NP = length(zloc$x)-1
        if(NP<2)
          {
            print("Must have at least 2 clicks on screen to move a symbol")
            focsiz =EXPreplot()
            buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
            print(zloc)
            zloc = list(x = NULL, y = NULL)
            next
            
          }
        else
          {

###  take last two clicks before clicking the RECT button
            p1 = list(x=c(zloc$x[(NP-1):NP]), y=c(zloc$y[(NP-1):NP]))
            
            SELX = which(NEWXY$x>min(p1$x) & NEWXY$x<max(p1$x) & NEWXY$y>min(p1$y) & NEWXY$y<max(p1$y))

            bx = NEWXY$x[SELX]
            by = NEWXY$y[SELX]
            mfxy = list(x=mean(bx),y=mean(by ))
            
            dis1 = sqrt( (bx-mfxy$x)^2 + (by-mfxy$y)^2)

            mult = sqrt(  (p1$x[1]-p1$x[2])^2 +(p1$y[1]-p1$y[2])^2) 
            
            dixplo  = mult
            DX = (bx-mfxy$x)/dis1
            DY = (by-mfxy$y)/dis1
            
            PX = bx+ dixplo* DX
            PY = by+ dixplo* DY
            
            NEWXY$x[SELX] = PX
            NEWXY$y[SELX] = PY 
            
            
          }

        focsiz =EXPreplot()
        buttons = RPMG::rowBUTTONS(labs, col = colabs, pch = pchlabs)
        zloc = list(x = NULL, y = NULL)
      }
##############################################################      
      iloc = locator(1, type = "p")
      zloc = list(x = c(zloc$x, iloc$x), y = c(zloc$y, iloc$y))
      Nclick = length(iloc$x)
      if (is.null(zloc$x)) {
        return(sloc)
      }
      K = RPMG::whichbutt(iloc, buttons)
    }

    return(NEWXY)

  }
