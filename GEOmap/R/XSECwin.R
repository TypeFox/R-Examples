XSECwin<-function(SW, iseclab=1, xLAB="A" , labs=c("DONE","REFRESH", "PS" ), width=10, demo=FALSE)
  {
    if(missing(demo)) { demo = FALSE }
    if(missing(width)) { width=10 }
    if(missing(labs)) { labs=c("DONE","REFRESH", "PS" ) }
    if(missing(xLAB)) { xLAB="A" }
    if(missing(iseclab))  {  iseclab=1 }
    


    
    
#######  
    
##### source("/home/lees/XSECwin.R")

    ###  labs=c("DONE","REFRESH", "XSEC",  "CONT", "width", "PS" )
    
  ### XSECwin( SW , iseclab, xLAB , labs, demo=FALSE  )   
    #####  
    
    iseclab  = 0
    secmat = NULL
    ncol = 100
    TPALS = c("rainbow", "topo", "terrain")
    colabs = rep(1, length=length(labs))
    pchlabs = rep(0,length(labs))
    FUN = match.fun(TPALS[1])
    pal = FUN(ncol)

    ##########  xsec plotting  function:
    YNreplot<-function(global.vars)
    {

      r = global.vars$SW$r
      depth = global.vars$SW$depth
      iseclab = global.vars$iseclab
      xLAB = global.vars$xLAB
      plot(r , -depth,  main=paste( iseclab, xLAB) , xlab="km", ylab="Depth", asp=1)
      a1 = range(r, ns.rm=TRUE)
      mtext(xLAB, side=3, at=a1[1])
      mtext(paste(sep="", xLAB, "'") , side=3, at=a1[2])
      
    }

  
   
    NLABS = length(labs)
    NOLAB = NLABS +1000  ## some large number
    
    
    BLABS = labs

    ilocstyle = -1
    CONT.FLAG = FALSE
    XSEC.FLAG = FALSE
    PS.FLAG =  FALSE

    WIN  = LASTwin= NULL
#####################################     global variables
  global.vars = list(

    SW =SW,
    iseclab=iseclab,
    xLAB = xLAB,
       CONT.FLAG = FALSE,
    XSEC.FLAG = FALSE,
    PS.FLAG =  FALSE,

   ilocstyle = ilocstyle,
    iloccol = rgb(1,0.6, 0.6),
    ilocnum = 1,
    MAINdev=NULL,
    tempbuttons = NULL,

    BLABS = BLABS ,
    NLABS = length(BLABS),
    NOLAB = NOLAB,
    
    WIN =WIN,
    LASTwin = LASTwin,
    KLICK = NULL,
    thebuts = FALSE
    )
##################################### 
      

  YN = YNreplot(global.vars)
  
    if(demo==TRUE)
      {
        return(NULL)

      }
       
    
 global.vars$MAINdev = dev.cur()
  

  global.vars$buttoncex = 0.8
  
  buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs, cex=global.vars$buttoncex )

  global.vars$MAINdev = dev.cur()


###  Get.Screens(2)
  dev.set(global.vars$MAINdev)
 
  
  u = par("usr")
  sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
  zloc =list(x=NULL, y=NULL)
 
  zenclick = length(zloc$x)


  global.vars$BLABS = BLABS
  global.vars$buttons = buttons
  global.vars$zloc = zloc
  global.vars$sloc = sloc
  global.vars$zenclick = zenclick
  global.vars$action="donothing"
  OLDglobal.vars = global.vars
   global.vars$OLDglobal.vars = global.vars
 
 iloc = RPMG::ilocator(global.vars$ilocnum ,COL=global.vars$iloccol ,NUM=FALSE , YN=length(global.vars$sel), style=global.vars$ilocstyle )
       
 Nclick = length(iloc$x)
   
    K =  RPMG::whichbutt(zloc , buttons)
    sloc = zloc

    
    while(TRUE)
      {
############   button actions
        iloc = RPMG::ilocator(global.vars$ilocnum ,COL=global.vars$iloccol ,NUM=FALSE , YN=length(global.vars$sel), style=global.vars$ilocstyle )
        
        
        Nclick = length(iloc$x)

###########   quit and break loop

       ##  print(Nclick)
        

        if(Nclick>0)
          {
#######  add last click to list of clicks, continue 
            zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
            global.vars$zenclick = length(zloc$x)
            K =  RPMG::whichbutt(iloc ,buttons)
            sloc = zloc
            
            
            if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
              {


                buttons = RPMG::rowBUTTONS(labs, col=rep(grey(.8), length(labs)), pch=rep("NULL", length(labs)))
                title("Return to Calling Program")
                
                break;
              }

###########   refresh the screen
            if(K[Nclick] == match("REFRESH", labs, nomatch = NOLAB))
              {
                YNreplot(global.vars)
                buttons = RPMG::rowBUTTONS(labs, col=colabs, pch=pchlabs)
                zloc = list(x=NULL, y=NULL)
                next
              }
###########   
            if(K[Nclick] == match("Next", labs, nomatch = NOLAB))
              {
                dev.set(dev.next())
                
                zloc = list(x=NULL, y=NULL)
                next
              }

###########   
            if(K[Nclick] == match("PS", labs, nomatch = NOLAB))
              {
                
                PS.FLAG = !PS.FLAG

                psname = RPMG::local.file(paste("xsec", xLAB, sep=""), "eps")
                
                P = round(par('din'), digits=2); 
                postscript(file= psname , width=P[1], height=P[2],
                           paper = "special", horizontal=FALSE, onefile=TRUE,print.it=FALSE)
                
                
                YNreplot(global.vars)
                dev.off();
                
                dev.set(global.vars$MAINdev)
                
                zloc = list(x=NULL, y=NULL)

                next
              }


          }
        

        ## print("here at end of while loop")  
      }
    

  }
