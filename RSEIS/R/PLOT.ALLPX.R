`PLOT.ALLPX` <-
function(t0, STNS, COMPS, YPX, PHASE=NULL, POLS=TRUE,  FILL=FALSE, FORCE=TRUE, cex=cex, srt=srt)
  {
    
    ##     YPX[[NPX]] =  list( pick=c(NH$info$yr[ipick], NH$info$jd[ipick], NH$info$hr[ipick], NH$info$mi[ipick], asec), kind="SHAPE", sta= STNS[ypick], comp=COMPS[ypick], PULS=PULS)

    ### match(JH$COMPS, "V")
    ###  match(JH$STN, "MAS")&match(JH$COMPS, "V")

    if(missing(FILL)) { FILL=FALSE }
    if(missing(FORCE)) { FORCE=TRUE }
     if(missing(PHASE)) {  PHASE=NULL }
     if(missing(POLS)) {  POLS=TRUE }
   
    if(missing(cex))    { cex=1 }
    if(missing(srt))    { srt=0 }

#########  if FORCE==true then if there is no match, plot picks anyway on first trace
   if(length(STNS)<=0) { return(0) }
    
    du = 1/length(STNS)
    
 ###  print(YPX)
  ### print(paste(sep=' ', "PLOT.WPX NPX=", length(YPX), "du=", du))

    ### print(cbind(STNS, COMPS))

    uphase = unique(YPX$phase)

    if(is.null(PHASE))
      {
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
        
###   M = RPMG::meshgrid(1:length(STNS), 1:length(COMPS) )
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
          
###   M = RPMG::meshgrid(1:length(STNS), 1:length(COMPS) )
          ypixA = (length(STNS)-wmd)*du

      }

    if(length(mpicks)<1) return()

    pols = YPX$pol[ mpicks]
    
    ypixB =  ypixA+du

    
    
    
    x1 = secdif(   t0$jd, t0$hr, t0$mi, t0$sec, YPX$jd[ mpicks], YPX$hr[mpicks ],YPX$mi[mpicks ], YPX$sec[mpicks ])
    
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
        segments(x1, 0, x1, 1, col=grey(.8) )
      }
    else
      {
        
        segments(x1, ypixA, x1, ypixB, col=pcol)
        
      ####  text(x1, ypixB, labels=LABS, col=pcol, pos=4, cex=cex, srt=srt)


#################  introduction of jittered labels


        ######   print(cbind(ypixB, LABS))


        tx1 =  split(x1,ypixB)
        tLABS = split(LABS,ypixB)
        typixB  = split(ypixB,ypixB)
        tpcol  = split(pcol,ypixB)


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

