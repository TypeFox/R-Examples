stress<-function(PPs=matrix(ncol=4, nrow=3) ,  Rview=c(-130, -50) ,  xscale=100, Stensor=matrix(ncol=3, nrow=3)  )
  {


    if(missing(xscale)) xscale = 30
    if(missing(Stensor)) { Stensor = NULL }
    if(missing(Rview))
      {
        Rview=c(-130, -50)

      }
    if(missing(PPs))
      {

        P1 = xscale*c(0, 1,1,0)
        P2 = xscale*c(1, 0,1,0)
        P3 = xscale*c(1, 1,0,0)
        
        PPs = rbind(P1, P2, P3)
        
        
      }


    
    
    ampvec<-function(a){return(sqrt(sum(a*a))) }

    TauNumber = 0
    TauOnLine = list()
    ES = NULL
    colt = 'black'
    cex = 1

    devCUBE  =  NULL
    devMOHR = NULL
    devAUX    =  NULL
    
    opendevs = dev.list()
    if(is.null(opendevs)  )
      {
        dev.new()
        devCUBE  = dev.cur()
        opendevs = dev.list()
      }
    else
      {

        devCUBE  = opendevs[1]

      }

    if(length(opendevs)>1)
      {
        devMOHR   =  dev.next()
        dev.set(devCUBE)   
      }
    else
      {
        dev.new()
        devMOHR   = dev.cur()
        opendevs = dev.list()
      }

    
    if(length(opendevs)>2)
      {      
        devAUX    =    devMOHR+1
      }
    
    
    typefaces  =  c("serif","sans serif", "script",
      "gothic english", "serif symbol" , "sans serif symbol")

    fontindeces = c("plain", "italic", "bold", "bold italic", "cyrillic")

    typeface = typefaces[1]
    fontindex = fontindeces[1]

    vfont = c(typeface, fontindex)

    pdcols = c("black", "darkmagenta", "forestgreen", "blueviolet",
      "tan3", "lightseagreen", "deeppink", "cyan3", "bisque3",
      "magenta1", "lightsalmon3", "darkcyan", "darkslateblue",
      "chocolate4", "goldenrod4", "mediumseagreen", "DeepSkyBlue1", "DarkSeaGreen" )

    if(is.vector(Rview))
      {
        ANG = Rview
        Rview  =    RFOC::ROTZ(ANG[1]) %*% RFOC::ROTX(ANG[2]) 

      }
    

    if(!is.null(Stensor))
      {
####  get eigenvalue decomposition of stress tensor
        if(is.matrix(Stensor))
          {
            ES = eigen(Stensor)
          }
        if(is.vector(Stensor))
          {
            if(length(Stensor)==3)
              {
                Stensor = diag(Stensor)  
              }
            if(length(Stensor)==6)
              {
                nstress = Stensor[1:6]
                MAT = matrix(ncol=3, nrow=3)
                
                MAT[ lower.tri(MAT, diag=TRUE)] = nstress
                MAT[ upper.tri(MAT) ] = MAT[lower.tri(MAT)]
                
                Stensor = MAT
                
                
              }
            
            ES = eigen(Stensor)
          }


        
        MAXSTRESS = TRUE
        MAXSHEAR = TRUE

        if(is.null(devMOHR) )
          {
            dev.new()
            devMOHR = dev.cur()
           ##  DoMohr(Stensor)
            dev.set(dev.prev())

          }
        else
          {            
            devMOHR = dev.set(devMOHR )
         ##    DoMohr(Stensor)
            dev.set(devCUBE)
          }

        
      }
    else
      {
        MAXSTRESS = FALSE
        MAXSHEAR = FALSE

      }

   

    BOX <-matrix(c(0,0,0,0,
                   0, 1, 0,0,
                   0, 1, 1,0,
                   0, 0, 1,0,
                   1,0,0,0,
                   1, 1, 0,0,
                   1, 1, 1,0,
                   1, 0, 1,0), ncol=4, byrow=TRUE)

    BOX = xscale*BOX
    
    AX = matrix(c(0,0,0,0,
      1, 0, 0,0,
      0, 0, 0,0,
      0, 1, 0,0,
      0,0,0,0,
      0, 0, 1,0), ncol=4, byrow=TRUE)

    AX = 1.5*xscale*AX

    

    Rax =  AX %*% Rview

    Rbox =   BOX %*% Rview

    Rp = PPs  %*% Rview
    
    headlen =xscale* .3/6
    len =xscale* .7/6
    basethick =xscale* 0.05/2
    headlip =xscale* .02/2
    aglyph = Z3Darrow(len = len , basethick =basethick , headlen =headlen , headlip=headlip )
    
    axcol = 'black'
    boxcol = 'blue'
    planecol = 'brown'
    
    rim = 0.05*xscale

    SHOWNORM = TRUE
    
    labs = c( "QUIT","DONE", "NORM", "MXnrm","MXshr", "LINE", "AXES", "STens", "TAU", "Replot", "Refresh" )


    
    colabs = rep(1, length=length(labs))
    pchlabs = rep(NA,length(labs))
    
    NLABS = length(labs)
    NOLAB = NLABS +1000  ## some large number


   
    gvars = list(
      Stensor=Stensor,
      ES=ES,
      xscale=xscale,
      rim = rim,
      Rax=Rax,
      Rbox=Rbox,
      axcol=axcol ,
      boxcol=boxcol,
      Rp=Rp,
      planecol=planecol,
      PPs=PPs,
      Rview=Rview,
      aglyph=aglyph,
      SHOWNORM=SHOWNORM,
      cex=cex,
      vfont = vfont,
      colt=colt,
      MAXSHEAR=MAXSHEAR,
      MAXSTRESS=MAXSTRESS,
      TauNumber=TauNumber,
      devMOHR  = devMOHR, 
      devCUBE  = devCUBE,
      devAUX   = devAUX 
      
      
      )

    plotnormvec<-function(B, PPs, xscale, Rview, aglyph)
      {
        L = list(x1 = PPs[3, 1], y1 = PPs[3, 2], z1 = PPs[3, 3], 
          x2 = PPs[3, 1] + xscale * B[1]/5, y2 = PPs[3, 2] + xscale * 
          B[2]/5, z2 = PPs[3, 3] + xscale * B[3]/5)
          RFOC::BOXarrows3D(L$x1, L$y1, L$z1, L$x2, L$y2, L$z2, aglyph = aglyph, 
            Rview = Rview, col = "green")
      }
    
    fresheng<-function(g)
      {
        NN = NORMvec(g$PPs, g$xscale, g$Rview,aglyph=g$aglyph, add=g$SHOWNORM)
        sigNORMmax = NN[1]^2*g$ES$values[1] + NN[2]^2 * g$ES$values[2]  +NN[3]^2 *  g$ES$values[3]
         tauSHEARmax  = NN[1]^2*NN[2]^2*(g$ES$values[1]-g$ES$values[2])^2 +
              NN[2]^2 *NN[3]^2 * (g$ES$values[2]-g$ES$values[3])^2  +
                NN[3]^2 *NN[1]^2 *  (g$ES$values[1] -g$ES$values[3])^2

            tauSHEARmax = sqrt(tauSHEARmax)
        g$NN = NN
        g$sigNORMmax = sigNORMmax
        g$tauSHEARmax = tauSHEARmax
        

        return(g)
      }

    replot<-function(g)
      {
        ##
      ####  print(paste("check xscale=",g$xscale))

      ####  print(g$devCUBE)
        dev.set(g$devCUBE)
      ##  print("plot 1")
        pstart(xscale=g$xscale)
        PLOTbox(g$Rax, g$Rbox, axcol=g$axcol , boxcol=g$boxcol )
        PLOTplane(g$Rp, planecol=g$planecol )
        NN = g$NN
        if(g$SHOWNORM)  plotnormvec(NN, g$PPs, g$xscale, g$Rview, g$aglyph)

        if(!is.null(g$Stensor))
          {
            ###########  this plots the legend in the upper right
            mohrleg(g$ES)
          }


        u = par("usr")
        posTOPLEFT = u[4]-0.1*(u[4]-u[3])
        shigh =  strheight("TEST" , units = "user", cex = g$cex, vfont = g$vfont)
        
        if(g$MAXSTRESS)
          {
           
            s1 = substitute(sigma[max] == sig, list(sig= g$sigNORMmax)  )
            text(u[1], posTOPLEFT , labels=s1, vfont=g$vfont, cex=g$cex, col=g$colt, xpd=TRUE, pos=4)
            
          }
        if(g$MAXSHEAR)
          {
            
            tlab = substitute(tau[max] == sig, list(sig=g$tauSHEARmax)  )
            
            text(u[1], posTOPLEFT-1.5*shigh , labels=tlab,
                 vfont=g$vfont, cex=g$cex, col=g$colt, xpd=TRUE, pos=4)
            
          }

        if(g$SHOWNORM){
          u = par("usr")
          pos1 = (u[3]+u[4])/2
          s1 = substitute(n[1] == sig, list(sig=NN[1])  )
          shigh =  strheight(s1, units = "user", cex = g$cex, vfont = g$vfont)

          text(u[1], pos1, labels=s1,
               vfont=g$vfont, cex=g$cex, col=g$colt, xpd=TRUE, pos=4)

          s2 = substitute(n[2] == sig, list(sig=NN[2])  )
          text(u[1], pos1-2*shigh, labels=s2,
               vfont=g$vfont, cex=g$cex, col=g$colt, xpd=TRUE, pos=4)

          s3 = substitute(n[3] == sig, list(sig=NN[3])  )
          text(u[1], pos1-4*shigh, labels=s3,
               vfont=g$vfont, cex=g$cex, col=g$colt, xpd=TRUE, pos=4)
          
        }
      
        if(g$MAXSTRESS & g$MAXSHEAR)
          {

             ##print("plot 2")
            dev.set(  g$devMOHR )
            DoMohr(g$Stensor, axis=c(1,4) )
            points(g$sigNORMmax, g$tauSHEARmax, pch=23, col='red', bg='yellow' )


            if(g$TauNumber>0)
              {
                for(ktau in 1:g$TauNumber)
                  {
                    temptau = g$TauOnLine[[ktau]]$tau
                    points(g$sigNORMmax, temptau, pch=21, col='blue', bg='gold', cex=.8 )

                  }
              }

            dev.set(g$devCUBE)
            
          }
      }

    
 

    zloc  = list(x=NULL, y=NULL)
    sloc = zloc
##   print("start")
    gvars$SHOWNORM = FALSE
    gvars = fresheng(gvars)
 ##      print("after")
      gvars$SHOWNORM = TRUE
    replot(gvars )
  
    gvars$devCUBE = dev.cur()
    
    
    u = par("usr")
    posTOPLEFT = u[4]-0.1*(u[4]-u[3])
    charheight  =  strheight("TEST" , units = "user", cex = cex, vfont = vfont)
    
    
    buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)

    zloc  = list(x=NULL, y=NULL)
    
    while(TRUE)
      {
        iloc = locator(1,type='p')
        K =  RPMG::whichbutt(iloc , buttons)
##### print(iloc)
        zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
        Nclick = length(iloc$x)

        if(Nclick<1)
          {
            return(NULL)
            
          }
        else
          {
        ####################################  change the plane if the point is near a node
        dis =  sqrt( (iloc$x-gvars$Rp[,1])^2 +  (iloc$y-gvars$Rp[,2])^2  )

        if(all(dis>gvars$rim))
          {
            m1 = 0

          }
        else
          {
            m1 = which.min(dis)
            
          }
        if(m1>0)
          {
            if(m1==1)
              {
                 basepoint = 3
                legpoints = c(7,4,2)
                
                
              
            }
            if(m1==2)
              {
                basepoint = 8
                legpoints = c(7, 4, 5)
                
              }
            if(m1==3)
              {
                
                basepoint = 6


                legpoints = c(7, 5, 2)

              }

            points(Rbox[c(legpoints,basepoint),1], Rbox[c(legpoints,basepoint),2], pch=c(rep(22, 3), 24), col=c(rep('red', 3), 'green') , bg=grey(.85) )
            
            Lp = locator(1)
            
            gvars$PPs = REplane(m1, Lp, gvars$PPs, gvars$Rbox, gvars$Rview, gvars$xscale)
            gvars$Rp = gvars$PPs  %*% gvars$Rview
            gvars = fresheng(gvars)
            replot(gvars )

              buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
            m1 = 0
            next
            
          }
          
            
          }
         ########################################################### 
########################################################### 
        ###########################################################  exit
        if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
          {

            buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.85), times=length(labs) ), pch=pchlabs)
            title(main="RETURN to MAIN")
            
                zloc  = list(x=NULL, y=NULL)
    
            return()
            break;
          }
########################################################### 
########################################################### 
########################################################### 
        if(K[Nclick] == match("QUIT", labs, nomatch = NOLAB))
          {

             buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.85), times=length(labs) ), pch=pchlabs)
             title(main="RETURN to MAIN")

             
                zloc  = list(x=NULL, y=NULL)
    
            return()
            break;
          }
########################################################### 
########################################################### 
########################################################### 

        
        if(K[Nclick] == match("Refresh", labs, nomatch = NOLAB))
          {
           ##  PPs = MOVEpt(PPs, Rbox, Rview, xscale ) 
         ##  print("pressed Refresh")
           TauNumber = 0
           TauOnLine = list()

           replot(gvars )

           buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
          ## print("pressed refresh")
           m1 = 0
               zloc  = list(x=NULL, y=NULL)
           next
          }


        
        if(K[Nclick] == match("Replot", labs, nomatch = NOLAB))
          {
           ##  PPs = MOVEpt(PPs, Rbox, Rview, xscale )
              replot(gvars)
            ##  NN = NORMvec(PPs, xscale, Rview)
           buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
           m1 = 0
      
         ##  print("pressed replot")
                zloc  = list(x=NULL, y=NULL)
      
           next
          }
########################################################### 
########################################################### 
########################################################### 
        
        if(K[Nclick] == match("LINE", labs, nomatch = NOLAB))
          {

            if(is.null(gvars$ES))
              {

                cat("Cannot estimate Shear Stress without a Stress Tensor input", sep="\n")
                cat("click on STens to input tensor first", sep="\n")
                

              }
            else
              {
            npoints = length(zloc$x)

            if(npoints>2)
              {
                
            ##    DUMPLOC(P1)
            ##    DUMPLOC(P2)

                kn = (npoints-1)
                
                gvars$P2 = list(x=zloc$x[kn] , y=zloc$y[kn])
                gvars$P1 = list(x=zloc$x[kn-1] , y=zloc$y[kn-1])

               TauNumber = TauNumber+1

                
               ## L = list(x=zloc$x[1:(npoints-1)] , y=zloc$y[1:(npoints-1)])
             
                NewTauOUT  =   tauline(gvars$Rp, gvars$P1, gvars$P2, gvars$Rview, gvars$ES, gvars$NN)

                NewTau=NewTauOUT$tau
              ##   replot()
                TauOnLine[[TauNumber]] = NewTauOUT
                acol = pdcols[TauNumber]
                
                arrows(gvars$P1$x, gvars$P1$y, gvars$P2$x, gvars$P2$y, length = 0.1, col=acol)
                text(gvars$P1$x, gvars$P1$y, labels=TauNumber, pos=2, cex=.8, col=acol)
                
               ##  buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
               ##  m1 = 0

                
                u = par("usr")
                ypos = u[3]+0.8*(u[4]-u[3]) - TauNumber*charheight*1.4
               
                  text(u[2],ypos , labels=substitute(tau[t] == sig, list(sig=NewTau, t=TauNumber)  ),
                 vfont=vfont, cex=cex, col=acol, xpd=TRUE, pos=2)
            
                cat(paste(sep=" ", "tau=",NewTauOUT$tau, "vec="), sep=" ")
                cat(paste(collapse=" ",NewTauOUT$vec), sep="\n\n")


                

                
              }
           
          }

             zloc  = list(x=NULL, y=NULL)
            next
           ##  PPs = MOVEpt(PPs, Rbox, Rview, xscale ) 
           ### print("pressed LINE")
           
          }
########################################################### 
########################################################### 
########################################################### 
       if(K[Nclick] == match("NORM", labs, nomatch = NOLAB))
          {
           ##  PPs = MOVEpt(PPs, Rbox, Rview, xscale ) 
            gvars$SHOWNORM = !gvars$SHOWNORM
            
             replot(gvars )

             buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
           m1 = 0
             zloc  = list(x=NULL, y=NULL)
            next
           
          }
########################################################### 
########################################################### 
########################################################### 

        
       if(K[Nclick] == match("MXnrm", labs, nomatch = NOLAB))
          {
           ##  PPs = MOVEpt(PPs, Rbox, Rview, xscale ) 

####  calculate the max normal stress

            gvars$MAXSTRESS = !gvars$MAXSTRESS

            if(is.null(gvars$Stensor)) {
              print("No Stress Tensor")
              gvars$MAXSTRESS =FALSE }
            
            replot(gvars)
            
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
            m1 = 0
             zloc  = list(x=NULL, y=NULL)
            next
            
          }

########################################################### 
########################################################### 
########################################################### 
        
     if(K[Nclick] == match("MXshr", labs, nomatch = NOLAB))
          {
           ##  PPs = MOVEpt(PPs, Rbox, Rview, xscale ) 

####  calculate the max normal stress

            gvars$MAXSHEAR = !gvars$MAXSHEAR

            if(is.null(Stensor)) {
              print("No Stress Tensor")
              gvars$MAXSHEAR =FALSE }
            
            replot(gvars)
            
            buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
            m1 = 0
             zloc  = list(x=NULL, y=NULL)
            next
            
          }
########################################################### 
########################################################### 
########################################################### 
     if(K[Nclick] == match("TAU", labs, nomatch = NOLAB))
          {


            print(gvars$PPs)
            
            BV = TriangleCenter(gvars$PPs[1,1:3],gvars$PPs[2,1:3], gvars$PPs[3,1:3] )
            CIRCview =   BV$Cinscribed  %*% gvars$Rview

            lines(CIRCview[,1], CIRCview[,2], col='purple')

            cview =    BV$Center %*%  gvars$Rview

            points(cview[1,1], cview[1,2])


            acols = rainbow(length(CIRCview[,1]))
            
            segments(  cview[1,1], cview[1,2]  , CIRCview[,1], CIRCview[,2], col=acols)
            
            vees = sweep(BV$Cinscribed , 2, BV$Center , "-")
            rees = sweep(vees , 1, apply(vees, 1,ampvec) , "/")

            TAUline =  gvars$NN[1]*rees[,1]*gvars$ES$values[1] +
              gvars$NN[2]*rees[,2] * gvars$ES$values[2]  +
                gvars$NN[3]*rees[,3] *  gvars$ES$values[3]

            #####   get an angle here, relative to downdip direction unless
            ########    normal is pointing up, then use X direction

            whTAU = which.max(TAUline)
           
           ###  ampTAU = ampvec(BV$Center[,1:3]-BV$Cinscribed[whTAU, 1:3])

            ####   get horizontal: cross product of normal(NN) and vertical vector:

            hozvec = AXB.prod(gvars$NN, c(0,0,1))
            Ahoz = ampvec(hozvec)
            if(Ahoz!=0)
              {
                hozvec =  hozvec/Ahoz
              }
            else
              {
                hozvec = c(0,1,0)
              }
            
            dotTAU = (hozvec[1] *rees[,1]  +  hozvec[2] *rees[,2]  + hozvec[3] *rees[,3])

            angTAU = 180*acos(dotTAU)/pi

            H1N1 = c(BV$Center[1:3]+ BV$r*hozvec, 0)
            
            VH1N1 =    H1N1  %*% Rview

             arrows(  cview[1,1], cview[1,2]  ,VH1N1[ 1], VH1N1[2], col='black', length=0.1, lty=2 )

             arrows(  cview[1,1], cview[1,2]  , CIRCview[whTAU,1],
                    CIRCview[whTAU,2], col='black', length=0.1, lwd=1.2 )

            
            umine = par("usr")

            jnum = 6
            ybase = umine[3]+charheight*1
            xbase = umine[1]+0.05*(umine[2]-umine[1])
            xdx = (umine[2]-umine[1])/(jnum+1)
            yh = 0
            for(j1 in 1:length(TAUline))
              {
                yh = yh + charheight*1.2
                if( (j1 %% jnum) == 0 ) yh = 0
                ypos = ybase + yh
                
                xpos = xbase+floor(j1/jnum)*xdx

                ftau = as.numeric( format(TAUline[j1], digits=3) )
                
                text(xpos, ypos , labels=substitute(tau[t] == sig, list(sig=ftau , t=j1)  ),
                     vfont=vfont, cex=cex, col=acols[j1], xpd=TRUE, pos=4)
                
                cat(paste(sep=" ", "tau=",ftau , "vec="), sep=" ")
                cat(paste(collapse=" ",format(rees[j1,], digits=5)), sep="\n\n")
                
          }


            if( is.null(gvars$devAUX) )
              {
                dev.new()
                gvars$devAUX = dev.cur()
              }
            else
              { 
                dev.set(gvars$devAUX )
              }
                plot(angTAU, TAUline, col=acols)
                abline(h=gvars$tauSHEARmax, col=grey(.85))
                abline(v=angTAU[whTAU ], lty=2, col=grey(.85))
                       
                dev.set(gvars$devCUBE)
              

            ## replot()
            
          ##   buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
           m1 = 0
             zloc  = list(x=NULL, y=NULL)
            next
           
          }



########################################################### 
########################################################### 
###########################################################
        
        if(K[Nclick] == match("STens", labs, nomatch = NOLAB))
          {   ####### STens


            cat(c("################################",
                  "Input stress tensor:  six or three values",
                  "s11 s12 s13 s22 s23 s33, as in:",
                  "s11 s12 s13",
                  "s12 s22 s23",
                  "s13 s23 s33",
                  "example: 15 0 0 10 0 5",
                  "or example: 15 10 5"),
                sep="\n")
            

            ANUMX = readline(prompt=">: ")

            nstress = as.numeric(unlist(strsplit(ANUMX, split=" " )))

            if(length(nstress)<3)
              {
                print("ERROR: need six distinct values separated by a blank")
                break
              }
            else
              {

                if(length(nstress)==3)
                  {
                    gvars$Stensor = diag(nstress)

                  }
                else
                  {
####   incase they input more than 6
                    nstress = nstress[1:6]
                    
                    MAT = matrix(ncol=3, nrow=3)

                    
                    MAT[ lower.tri(MAT, diag=TRUE)] = nstress
                    MAT[ upper.tri(MAT) ] = MAT[lower.tri(MAT)]

                    gvars$Stensor = MAT

                  }
                cat("INPUT TENSOR:", sep="\n")
                
                for(is in 1:3){ cat(paste(gvars$Stensor[is,], collapse=" "), sep="\n" ) }
                cat("", sep="\n")

                gvars$ES = eigen(gvars$Stensor)
                gvars$MAXSTRESS = TRUE
                gvars$MAXSHEAR = TRUE

                gvars = fresheng(gvars)

                if(is.null(gvars$devMOHR)  )
                  {
                    dev.new()
                    gvars$devMOHR = dev.cur()

                  }
                else
                  {    
                     dev.set( gvars$devMOHR)
                   }
                    DoMohr(Stensor)
                    dev.set(dev.prev())

                
                replot(gvars)
                
                buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
                m1 = 0

              }

             zloc  = list(x=NULL, y=NULL)
            
            next


          }  

########################################################### 
########################################################### 
###########################################################

     
        if(is.null(zloc$x)) { return(sloc) }
        K =  RPMG::whichbutt(iloc , buttons)


 
      }#########   end while for clicking






  }  ##########   end function definition
