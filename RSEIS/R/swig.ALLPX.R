`swig.ALLPX` <-
function(t0, STNS, COMPS, YPX, PHASE=NULL, POLS=TRUE,  FILL=FALSE, FORCE=TRUE, cex=cex, srt=srt)
  {
    ############  plot picks on a swig window
    
    if(missing(FILL)) { FILL=FALSE }
    if(missing(FORCE)) { FORCE=TRUE }
     if(missing(PHASE)) {  PHASE=NULL }
     if(missing(POLS)) {  POLS=TRUE }
   
    if(missing(cex))    { cex=1 }
    if(missing(srt))    { srt=0 }

#########  if FORCE==true then if there is no match, plot picks anyway on first trace
   if(length(STNS)<=0) { return(0) }

    ####   STNS is a list of showing stations
    ####   COMPS = vector of associated components
     ####   YPX =  list of picks

   ####    print(YPX)

  ####   print(paste("PHASE=",PHASE))
  ####   print(paste( "POLS=", POLS,  "FILL=",FILL, "FORCE=",FORCE, "cex=",cex, "sr=",srt))
    
      du = 1/length(STNS)
    
    uphase = unique(YPX$phase)


####    print(paste("input parameters  PHASE=", PHASE, "POLS=", POLS,  "FILL=", FILL,
####                "FORCE=",FORCE,"cex=", cex, "srt=", srt, collapse=" "))

    

    if(is.null(PHASE))
      {
       ## print("PHASE is NULL")
        ###########  match the pick with the component only
       ####   m1 = match(STNS , YPX$name)

       #### uphase = unique( YPX$phase)
        
       #### m1 = match( YPX$name, STNS )
       #### m2 = match(STNS, YPX$name)
        
      ####  wmd = m1
        
      ####  m1 = which(!is.na(wmd))
      ####  wmd = m1
        
       #### mpicks = m2[YPX$onoff[m2]>=0]
        wmd = vector()
        mpicks = vector()

        for(jj in 1:length(STNS))
          {
            m1 = match( YPX$name, STNS[jj] )
            m2 = which(!is.na(m1))
            mpicks = c(mpicks, m2)
            wmd  = c(wmd, rep(jj, length( m2)))

          }

        moff = YPX$onoff[ mpicks]>=0

        mpicks  = mpicks[moff]
        wmd =     wmd[moff]
        qres = YPX$res[ mpicks]
        qerr = YPX$err[mpicks ]
        LABS = YPX$phase[mpicks]
        
###   
        ypixA = (length(STNS)-wmd)*du

        

      }
    else
      {


        tag1 = paste(sep=".",YPX$name, YPX$comp)
        tag2 = paste(sep=".",STNS, COMPS)

        
          wmd = vector()
          mpicks = vector()

          for(jj in 1:length(STNS))
            {
              m1 = match( tag1, tag2[jj] )
              m2 = which(!is.na(m1))
              mpicks = c(mpicks, m2)
              wmd  = c(wmd, rep(jj, length( m2)))

            }

          moff = YPX$onoff[ mpicks]>=0

          
          mpicks  = mpicks[moff]
          wmd =     wmd[moff]
          qres = YPX$res[ mpicks]
          qerr = YPX$err[mpicks ]
          LABS = YPX$phase[mpicks]
          
###   
          ypixA = (length(STNS)-wmd)*du

      }

    if(length(mpicks)<1) return()




    
    pols = YPX$pol[ mpicks]
    
    ypixB =  ypixA+du

    
    
    
    x1 = YRsecdif(   t0$jd, t0$hr, t0$mi, t0$sec, YPX$jd[ mpicks], YPX$hr[mpicks ],YPX$mi[mpicks ], YPX$sec[mpicks ], yr1=t0$yr, yr2=YPX$yr[ mpicks] )

   ## print(x1)

    
#############  set the color of the pix 
    if(is.null(YPX$col))
      {
        pcol=rep("springgreen4", length(YPX$sec[mpicks]))
        phas = YPX$phase[mpicks ]
        pcol[phas=="P"] = "violetred2"
        pcol[phas=="S"] = "deepskyblue4"
       
        
      }
    else
      {
        pcol=YPX$col[mpicks]


      }

  ###   print(paste(sep=" ", "PLOT.ALLPX", pcol))


    #############  plot main pick



    
    if(FILL==TRUE)
      {
#############  draw the pick down the whole window
        ypixA = rep(0, times=length(ypixA))
         ypixB = rep(1, times=length(ypixB))
        
      }

    
      segments(x1, ypixA, x1, ypixB, col=pcol)
        
      ####  text(x1, ypixB, labels=LABS, col=pcol, pos=4, cex=cex, srt=srt)


#################  introduction of jittered labels


######   print(cbind(ypixB, LABS))
    noLAB = FALSE
    
    noLAB = all(YPX$flg[!is.na(YPX$flg)] =="N")


    
   #### print(mpicks)
   #### print(YPX)
    

    if(!noLAB)
      {

        jlab = YPX$flg[mpicks]
        wlab = which(jlab!="N")
        
       ####   print(wlab)
        
        nx1 =x1[wlab]
        nLABS=LABS[wlab]
        nypixB=ypixB[wlab]
        npcol=pcol[wlab]
        
        tx1 =  split(nx1,nypixB)
        tLABS = split(nLABS,nypixB)
        typixB  = split(nypixB,nypixB)
        tpcol  = split(npcol,nypixB)

          ###  print(c(length(wlab), length(tx1)))
          ### print(tx1)

        
##############  typixB are the categories
        for(ipixB in 1:length(typixB))
          {
            gx1 =   tx1[[ipixB]]
            gypixB = typixB[[ipixB]]
            gLABS = tLABS[[ipixB]]
            gpcol = tpcol[[ipixB]]
            
            RMAT =   RPMG::textrect( gx1, gypixB,
              gLABS , textcol=gpcol   ,xpd=TRUE, add=FALSE, font=1, cex=cex, brd = 0.03 )

            newjitx = jitter.lab(RMAT[,1]  , RMAT[,3]-RMAT[,1])
            newy = (gypixB)-newjitx*(RMAT[,4]-RMAT[,2])        

            
            RMAT =   RPMG::textrect(gx1 , newy,gLABS , textcol=gpcol,
              xpd=TRUE, add=TRUE, font=1, cex=cex, brd = 0.03 )
            
          }
      }
    
    

         #############    if a duration is given, mark it with a bar and a hook
   
    segments( x1, ypixB, x1+qres, ypixB, col=pcol  )
    segments( x1+qres, ypixB, x1+qres, ypixB-.2*du, col=pcol  )
    
    
    
#############    if an err is given, mark it with a small error bat
    
    if(!any(is.na(qerr) ))
      {
           

           
             segments( x1-qerr, ypixB-.1*du, x1+qerr, ypixB-.1*du, col=pcol  )
             segments( x1+qerr, ypixB-.05*du, x1+qerr, ypixB-.15*du, col=pcol  )
             segments( x1-qerr, ypixB-.05*du, x1-qerr, ypixB-.15*du, col=pcol  )

             
         }

    if(POLS)
      {
        
        points(x1[pols=="U"], ypixB[pols=="U"]-.25*du, pch=2, col=pcol[pols=="U"], cex=cex)
        points(x1[pols=="D"], ypixB[pols=="D"]-.25*du, pch=6, col=pcol[pols=="D"], cex=cex)
        

      }
       
        
      
    
  }

