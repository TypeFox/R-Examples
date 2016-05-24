`colwheel` <-
function(v=1, BACK="black")
  {
    
    
    if(missing(v)) {  v=1 }
    if(missing(BACK)) {  BACK="white" }


    opar = par(no.readonly=TRUE)

    
    require(RPMG)
    
  ######  labs = c("DONE", "REFRESH", "SHOW", "CLEAR", "0.5", "0.6", "0.7", "0.8", "0.9", "1"  )
    labs = c("DONE", "REFRESH", "SHOW", "CLEAR", "BACK.BW", "BACK", "CHANGEV")


    changev<-function(g, TE="Click in Bar for change of v")
  {

    pct1 = .25
    pct2 = .35
    
    u = par("usr")
    x = c(u[1] , u[2])
    y =  c(u[3] , u[4])

    curx = 

    x = c(x[1]+pct1*diff(x), x[2]-pct1*diff(x))
    y = c(y[1]+pct2*diff(y), y[2]-pct2*diff(y))

    rect(x[1], y[1], x[2], y[2], col=rgb(1,.7, .7), border="blue", lwd=2)
    text(x[1], mean(y), labels=0, pos=4)
    text(x[2], mean(y), labels=1, pos=2)

     curx = x[1] + diff(x)*g$v
    text(mean(x), mean(y), TE)

    segments(curx, y[1] , curx, y[2] , col='red', lwd=2)
    
    L = locator(1, type="p" )

    new = (L$x-x[1])/diff(x)
    if(new>1) new = 1
    if(new<0) new = 0
    invisible(new)
  }


    

    
    NLABS = length(labs)
    NOLAB = NLABS +1000
    
    pchlabs = c(rep(1,length(labs)))
    if(identical(BACK, "black"))
      {
        cmain = rgb(.7,.7,.7)
        par(bg=rgb(0.1,0.1,0.1), fg=rgb(1,1,1), col.axis=rgb(1,1,1), col.lab=rgb(1,1,1), col.main=cmain , col.sub=rgb(1,1,1)   )
        colabs = c(rep('white',length(labs)))
        
        
      }
    else
      {
        cmain = rgb(1-.7, 1-.7, 1-.7)
        par(bg=rgb(1,1,1), fg=rgb(0,0,0),col.axis=rgb(0,0,0),col.lab=rgb(0,0,0), col.main=cmain ,col.sub=rgb(0,0,0)   )
        colabs = c(rep('black',length(labs)))
        
      }
  
    BIGMESH = VVwheel( v=v)

    RY = BIGMESH$RY
    rye = range(c(BIGMESH$RX  , BIGMESH$RX+diff(BIGMESH$RY)))
    littley = .2*abs(diff(BIGMESH$RY))


    gv = list(BIGMESH=BIGMESH, RY=RY   ,v=v, labs=labs, colabs=colabs, pchlabs=pchlabs)

    CWreplot<-function(gv)
      {
        VVwheel(BIGMESH=gv$BIGMESH, v=gv$v)
        buttons = rowBUTTONS(gv$labs, col=gv$colabs, pch=gv$pchlabs)
        return(buttons)
      }
    
    
    
    buttons = CWreplot(gv)
     zloc = list(x=NULL, y=NULL)
  ########  iloc = locator(1, col=grey(.5), type='p')
  ########  zloc = iloc
  ########  Nclick = length(iloc$x)
    
 ########   Nclick = length(zloc$x)
    
   ######## if(is.null(zloc$x)) { return(NULL) }
 ########   K = whichbutt(zloc ,buttons)
    
    sloc = zloc
    DF = NULL
    MAINdev = dev.cur()
    PALdev = NULL

    while(TRUE)
      {
######################################################################
#######   rect( BIGMESH$RX[1], rye[1], BIGMESH$RX[2], rye[2], border='white', col=NA)
#######    BIGMESH$RX
#######            rye = range(c(BIGMESH$RX  , BIGMESH$RX+diff(BIGMESH$RY)))

        iloc = locator(1,type='p')
        zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
        Nclick = length(iloc$x)
        
        if(is.null(iloc$x)) { par(opar);  return(zloc) }
        K =  whichbutt(iloc ,buttons)
        
    

        if(iloc$x>=BIGMESH$RX[1] & iloc$x<=BIGMESH$RX[2] & iloc$y>=rye[1] & iloc$y<=rye[2])
          {


            Wcols = wheelrgb(iloc, gv$v, gv$RY)
            rect(BIGMESH$RX[1],rye[1]-littley/2, BIGMESH$RX[1]+littley, rye[1]-littley, col=Wcols, border=NA, xpd=TRUE)
          }

        
        if(K[Nclick] == match("DONE", labs, nomatch = NOLAB))
          {
            nloc=length(zloc$x)-1
            zloc = list(x=zloc$x[1:nloc], y=zloc$y[1:nloc])
            thecols = wheelrgb(zloc, gv$v, gv$RY)
            ###   cat(file="",paste(sep="", nam, "newcols=c(", paste(thecols, collapse=","), ")") , fill=TRUE)
            cprint(thecols)
            
            DF = thecols
            zloc = list(x=NULL, y=NULL)
            break;
          }
 ######################################################################
       
        if(K[Nclick] == match("REFRESH", labs, nomatch = NOLAB))
          {
            buttons =CWreplot(gv)
            next
          }
######################################################################
   
        if(K[Nclick] == match("CHANGEV", labs, nomatch = NOLAB))
          {

            newv =  changev(gv, TE="Click in Bar for change of v")
            gv$v = newv

BIGMESH = VVwheel( v=v)

    RY = BIGMESH$RY
    rye = range(c(BIGMESH$RX  , BIGMESH$RX+diff(BIGMESH$RY)))
    littley = .2*abs(diff(BIGMESH$RY))

gv$BIGMESH=BIGMESH
gv$RY=RY
          
            buttons =CWreplot(gv)
            next
          }


        
######################################################################

        if(K[Nclick] == match("BACK.BW", labs, nomatch = NOLAB))
          {
            if( identical(BACK, "black") ) { BACK="white" } else{  BACK="black" }

            if(identical(BACK, "black"))
              {
                cmain = rgb(.7,.7,.7)
                par(bg=rgb(0.1,0.1,0.1), fg=rgb(1,1,1),col.axis=rgb(1,1,1),col.lab=rgb(1,1,1),col.main=cmain ,col.sub=rgb(1,1,1)   )
                gv$colabs = c(rep('white',length(labs)))
                gv$pchlabs = c(rep(1,length(labs)))
                
              }
            else
              {
                cmain = rgb(1-.7, 1-.7, 1-.7)
                par(bg=rgb(1,1,1), fg=rgb(0,0,0),col.axis=rgb(0,0,0),col.lab=rgb(0,0,0),col.main=cmain ,col.sub=rgb(0,0,0)   )
                gv$colabs = c(rep('black',length(labs)))
                gv$pchlabs = c(rep(1,length(labs)))
              }
            
            
           buttons = CWreplot(gv)
            next
          }
######################################################################
        if(K[Nclick] == match("BACK", labs, nomatch = NOLAB))
          {
            
            nloc=length(zloc$x)-1
            if(nloc<1) next
            zloc = list(x=zloc$x[nloc], y=zloc$y[nloc])
            
            thecols = wheelrgb(zloc, gv$v, gv$RY)
            
            aa = col2rgb(thecols)
            
            print(paste(sep=" ", thecols, paste(aa, collapse=" ") ))
            aa = col2rgb(thecols)

            bb = c(255-aa[1], 255-aa[2], 255-aa[3])
            
            aacomp = rgb(bb[1], bb[2], bb[3], maxColorValue=255)

            print(paste(sep=" ", aacomp, paste( bb, collapse=" ")))
            
            cmain = rgb(.7,.7,.7)
            BACK = thecols
            par(bg=thecols, fg=aacomp ,col.axis=rgb(1,1,1),col.lab=rgb(1,1,1),col.main=cmain ,col.sub=rgb(1,1,1)   )
            gv$colabs = c(rep(aacomp,length(labs)))
            gv$pchlabs = c(rep(1,length(labs)))
            
            buttons = CWreplot(gv)
            
            zloc = list(x=NULL, y=NULL)
            next
            
          }
######################################################################

        
        if(K[Nclick]==match("SHOW", labs, nomatch = NOLAB) )
          {
            nloc=length(zloc$x)-1
            zloc = list(x=zloc$x[1:nloc], y=zloc$y[1:nloc])
            thecols = wheelrgb(zloc, gv$v, gv$RY)
            print(thecols)

            ###   cat(file="",paste(sep="", nam, "newcols=c(", paste(thecols, collapse=","), ")") , fill=TRUE)
            cprint(thecols)
            
            DF = thecols

            print(PALdev)
            

            if( is.null(PALdev) )
              { dev.new();
                PALdev=dev.cur();
                print(PALdev) }
            else
              {
                dev.set(PALdev)
              }
            SHOWPAL(thecols, BACK=BACK, NAME=TRUE)
##### 
            dev.set( MAINdev)
            next
            
          }
######################################################################

        if(K[Nclick]==match("CLEAR", labs, nomatch = NOLAB) )
          {
            
             DF = NULL
            buttons = CWreplot(gv)
           
            zloc = list(x=NULL, y=NULL)
             next
             
          }

 ######################################################################
      


      }
    par(opar)
    
    if(!is.null(DF)) invisible(DF)
        
}

