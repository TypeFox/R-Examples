`PMOT.drive` <-
function(temp,  dt, pmolabs=c("Vertical", "North", "East"), STAMP="", baz=0 )
  {

    if(missing(STAMP)) { STAMP = " " }
    if(missing(baz)) { baz=0 }
    if(missing(pmolabs)) {pmolabs=c("Vertical", "North", "East")  }
    
    TPALS = c("rainbow", "topo.colors", "terrain.colors", "JGRAY", "tomo.colors")
    APALS = c("rainbow", "topo", "terrain", "JGRAY", "tomo")
    ADDBUTS = c("More" )
  
    rotlabs=c("Vertical", "Radial", "Transvers")
    vnelabs=c("Vertical", "North", "East")
    
    labs = c("DONE", "Angles", "PTS", "PRINC", "LOCS", "Postscript", "ROTATE", APALS, ADDBUTS )
    BLABS = labs
    NLABS = length(labs)
    NOLAB = NLABS +1000


    
###  FUN = match.fun(TPALS[1])
     pal = RPMG::Gcols(plow=0, phi=0,  N=100, pal=TPALS[1])
    scale.def = 0
   
    colabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
    pchlabs = c(rep(1,2) , rep(2, length(APALS) ), rep(4,length(ADDBUTS) ))
 
    gridon = FALSE
 
  
    NSEL = 1

    ROTATEFLAG = FALSE

    ddim = dim(temp)
    txlen = (ddim[1]-1)*dt

    
    atemp = temp

    if(baz!=0)
      {
        rbaz = grotseis(baz, flip=FALSE)
        btemp  = atemp  %*%  rbaz
      }
    else
      {
        
         btemp = atemp
      }
    

  ###  X11()
###
    ilocstyle = -1
    ######################################################  global variables
    global.vars = list(AMAT=atemp,
      BMAT=btemp,
      ROTATEFLAG = FALSE,
      dt=dt,
      BLABS=pmolabs,
      pal=pal,
      ilocstyle = ilocstyle,
      iloccol = rgb(1,0.6, 0.6),
      ilocnum = 1,
      PTS = FALSE,
      ptscol="brown",
      ptssize = 1,
      ptspch = 3, stamp=STAMP)


    ######################################################
    ######################################################  replot
    ######################################################
    DOreplot<-function()
    {

      if(global.vars$ROT)
        {
          showmat = global.vars$BMAT
          global.vars$pmolabs = rotlabs
        }
      else
        {
          showmat = global.vars$AMAT
          global.vars$pmolabs = vnelabs
        }
     
      sx = complex.hodo(showmat,
        dt=global.vars$dt,
        labs=global.vars$pmolabs,
        COL=global.vars$pal,
        STAMP=global.vars$stamp)


      if(global.vars$PTS)
        {
          addpoints.hodo(showmat, global.vars$dt, sx,  pch=global.vars$ptspch, col=global.vars$ptscol)
        
        }
      
      invisible(sx)
    }
    ######################################################
    ######################################################
    ######################################################


    
    YN = DOreplot()
    ##  sx = complex.hodo(temp, dt=dt, labs=pmolabs, COL=pal, STAMP=STAMP)

    global.vars$buttoncex = 0.8
    
    buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs, cex=global.vars$buttoncex )
    
    global.vars$MAINdev = dev.cur()

   ##   buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
      zloc = list(x=NULL, y=NULL)
    sloc = zloc

    
    while(TRUE)
      {
        iloc = RPMG::ilocator(global.vars$ilocnum ,COL=global.vars$iloccol ,NUM=FALSE , style=global.vars$ilocstyle )
        
        Nclick = length(iloc$x)
           
          if(Nclick>0)
            {
              zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
              zenclick = length(zloc$x)
              K =  RPMG::whichbutt(iloc ,buttons)
              sloc = zloc
            }
          else
            {
              Nclick = 0
              K = 0
              buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)))
              title("Return to Calling Program")
              break;
            }
     


        
        if(K[Nclick] == match("DONE", BLABS, nomatch = NOLAB))
          {
            buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)), cex=global.vars$buttoncex)
            title("Return to Calling Program")

            
            break;
          }
       
        if(K[Nclick] == match("REPLOT", BLABS, nomatch = NOLAB))
        {
          YN = DOreplot()
          ## print("clicked REPLOT")
          buttons = RPMG::rowBUTTONS(global.vars$BLABS, col=global.vars$colabs, pch=global.vars$pchlabs, cex=global.vars$buttoncex)
          Nclick = 0
          K = 0
          zloc = list(x=NULL, y=NULL)
        
          next;
        }
      

        
      ####################   POSTSCRIPT  ##################
      
        if(K[Nclick] == match("Postscript", BLABS, nomatch = NOLAB))
        {

          print("Start postscript plot.ts")
          plfname = RPMG::local.file("pmot","eps")
          jdev = dev.cur()
          RPMG::jpostscript("pmot")
          YN = DOreplot()
          
           print("Done creating postscript")
          dev.off()
          dev.set(jdev)
          zloc = list(x=NULL, y=NULL) 
        }

        if(K[Nclick] == match("ROTATE", BLABS, nomatch = NOLAB))
          {
            global.vars$ROTATEFLAG = !global.vars$ROTATEFLAG
           
            YN = DOreplot()
          ###  sx = complex.hodo(temp, dt=dt, labs=pmolabs, COL=pal, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
            zloc = list(x=NULL, y=NULL) 
            

          }


        if(K[Nclick]==match("Angles", BLABS, nomatch = NOLAB) & (zenclick-1)>=2 )
          {
            
            
            print(zloc$x[1:(zenclick-1)])
            LN = length(zloc$x[1:(zenclick-1)])
            LN = 2*(floor(LN/2))
            ###  use only pairs of clicks
            
            sn1 = seq(from=1, to=LN-1, by=2)
            sn2 = sn1+1
            segments(zloc$x[sn1], zloc$y[sn1], zloc$x[sn2], zloc$y[sn2], col="black")
            print(STAMP)
            a1 = 180*atan2(zloc$y[sn2]-zloc$y[sn1], zloc$x[sn2]-zloc$x[sn1])/pi
            print(a1)
            zloc = list(x=NULL, y=NULL) 
          
          }

        ###########################################################
        if(K[Nclick]==match("PRINC", BLABS, nomatch = NOLAB))
          {

            
            print(zloc$x)
            LN = length(zloc$x[1:(zenclick-1)])
            if(LN>=2)
              {
                Timex =  zloc$x[1:(zenclick-1)]
                TEES = RPMG::RESCALE(Timex , 0,  txlen , 0, 3)
                NEES = round( RPMG::RESCALE(Timex , 0,  ddim[1] , 0, 3) )

                print(NEES)

                
                tsel = NEES[1]:NEES[2]
                NOR = temp[ tsel, 2 ]
                EAS  = temp[tsel, 3  ]
                VER  = temp[tsel, 1  ]
                len = length(NOR)

                mEAS = EAS-mean(EAS)
                mNOR = NOR - mean(NOR)

                covtem = var(cbind(mEAS, mNOR)  )
                
                eg=eigen(covtem, symmetric = TRUE )

                efact = sqrt(eg$values)
                
                dev.new()

               ### par(mfrow=c(2,1))
                
                plot(mEAS, mNOR ,asp=1,  type='l')
                
                vecs = list(x=mEAS[2:len]-mEAS[1:(len-1)],
                  y = mNOR[2:len]-mNOR[1:(len-1)])
                arrows(0,0 ,efact[1]*eg$vectors[1,1] , efact[1]*eg$vectors[2,1], col='red')

               princdir =  90-180*atan2(eg$vectors[2,1], eg$vectors[1,1])/pi
               princdir =  RPMG::fmod(princdir , 180)
                
                slop = paste(global.vars$stamp , format(princdir))
                cat(slop, sep="\n")
                title(main=slop )

                
               ###    arrows(0,0 ,efact[2]*eg$vectors[2,1] , efact[2]*eg$vectors[2,2], col='blue')

                

              ###  plot(vecs$x, vecs$y, type='n', asp=1)
                
               ### arrows(0, 0, vecs$x, vecs$y, length=0.1)
                
               ### plot(EAS-mean(EAS) , NOR-mean(NOR),
               ###      type='n', asp=1)
                
               ### arrows(0, 0, EAS-mean(EAS), NOR-mean(NOR) , length=0.1)
                
                
                locator(1)
                dev.off()
                dev.set(global.vars$MAINdev)

                
              }
            
            
            
          }

###########################################################


        
        
        if(K[Nclick]==match("LOCS", BLABS, nomatch = NOLAB))
          {
            print(zloc$x)
            LN = length(zloc$x[1:(zenclick-1)])
            sn1 = seq(from=1, to=LN, by=1)
            print(STAMP)

            plt1 = 1+floor(zloc$x[sn1])
           
            print(zloc$x[sn1])
            print(plt1)

            ids = idpoints.hodo(temp, YN, zloc$x[sn1], zloc$y[sn1])
            
            if(length(ids)>=2)
              {
                print(ids)

                t1 = ids[1]*dt
                t2t1 = dt*(ids[2]-ids[1])
                
                addpoints.hodo(temp, dt, YN,  flag=c(ids[1]:ids[2]) , pch=3, col="brown")
                
                femp = temp[c(ids[1]:ids[2]), ]
                covtem = var(femp)
                eg=eigen(covtem, symmetric = TRUE )
                
                arrows(.5, .5, .5+0.3*eg$vectors[3,1] , .5+0.3*eg$vectors[2,1])
                
                
                alpha=180*atan2(eg$vectors[2,1], eg$vectors[3,1])/pi
####az=90-alpha
                
                az = alpha
                inci=180*atan2(eg$vectors[1,1], sqrt(eg$vectors[2,1]^2+eg$vectors[3,1]^2))/pi
                rateig  = 1 - ((eg$values[2]+eg$values[3])/(2*eg$values[1]))
                
                plop = paste(sep="  ", STAMP, "Az=", format.default(az, width=6,digits=4, trim=FALSE),
                  "Inc=", format.default(inci, width=5, digits=4, trim=FALSE),
                  "Rat=", format.default(rateig, width=5,digits=4, trim=FALSE),
                  "T1=", format.default(t1, width=5,digits=4, trim=FALSE),
                  "T2-T1=", format.default(t2t1, width=5,digits=4, trim=FALSE)
                  )
                print(plop)
                text(0, 1.05, labels=plop, adj=0)
              }
             zloc = list(x=NULL, y=NULL) 
          }

        if(K[Nclick]==match("PTS", BLABS, nomatch = NOLAB))
          {
           global.vars$PTS = !global.vars$PTS
           
            YN = DOreplot()
          ###  sx = complex.hodo(temp, dt=dt, labs=pmolabs, COL=pal, STAMP=STAMP)
            buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
            zloc = list(x=NULL, y=NULL) 
            

          }


           
        if(K[Nclick]==match("More", BLABS, nomatch = NOLAB))
          {
            YN = DOreplot()
             buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
            zloc = list(x=NULL, y=NULL) 
            
           
          }
        if( length(which(K[Nclick] == match(APALS, labs, nomatch = NOLAB)))>0 )
          {
            J = match(labs[K[Nclick]] ,  APALS   )
            
            ##FUN = match.fun(TPALS[J])
            ##  pal = FUN(NCOL)

            
            global.vars$pal = RPMG::Gcols(plow=0, phi=0,  N=100, pal=TPALS[J])
            
            YN = DOreplot()
            buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
             zloc = list(x=NULL, y=NULL) 
          }





        

      }

    print("DONE with PMOT")
    
    
  }

