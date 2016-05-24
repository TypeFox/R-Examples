SELBUT<-function(OPTS, onoff=1, ocols ="white", default="opt"  )
  {

##    source("/home/lees/NEWseis/BUTREPLOT.R")

    Nopts = length(OPTS)
  
    if(missing(onoff)) { onoff =  rep(1 , Nopts )  }
    if(missing(default)) { default=OPTS  }

    
  
    stdlab = c("OKAY", "CANCEL", "DEFAULT", "REVERT", "ASORT", "SEL", "ADD", "SUB", "NONE", "ALL" )

    BLABS = c(stdlab)
    NLABS = length(BLABS)
    NOLAB = NLABS +1000

    colabs = rep(1,length(BLABS))
    
    pchlabs = rep(4,length(BLABS))

    defaultcol = grey(.9)
    addcol = rgb(1, .9, .9)
    subcol = rgb(.9, .9, 1)


    ocols = rep(defaultcol , Nopts)
    ocols[onoff==1] = addcol

    ORIGonoff =onoff
    ORIGcols = ocols
    
    
    YN = BUTREPLOT(OPTS, cols=ocols)
    
    buttons = RPMG::rowBUTTONS(BLABS, col=colabs, pch=pchlabs)
    u = par("usr")
    sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
    zloc =list(x=NULL, y=NULL)
    gvars = list(BLABS=BLABS, opts=OPTS, onoff=onoff, ocols=ocols, default=default )
    
    gclick<-function(gvars, zloc, YN)
      {
        klick = gvars$zenclick-1
        if(klick<1) return(NULL)
        
        thex = zloc$x[1:klick]
        they = zloc$y[1:klick]

        flag = thex>=YN$rx[1] & thex<=YN$rx[2] & they>=YN$ry[1] & they<=YN$ry[2]
        if(all(!flag)) { return(NULL) }

        
        thex = thex[flag]
        they = they[flag]
        
        klick = length(thex)
        w = vector()
        if(klick>0)
          {
            

            for(i in 1:klick)
              {
                w[i] = which.min( (YN$M$x+YN$dx/2 -thex[i])^2 +  (YN$M$y+YN$dy/2-they[i])^2)
                
              }
          }
        else
          {
            return(NULL)

          }
        return(w)
      }
    
    
    while(TRUE)
      {
    ##    iloc = RPMG::ilocator(1, COL=rgb(1,0.6, 0.6), NUM=FALSE , YN=length(gvars$sel), style=-1)
         iloc = locator(1, type='p',  pch=21, cex=3,  col='red', bg='yellow')
        Nclick = length(iloc$x)
        
     if(Nclick>0)
        {
          #######  add last click to list of clicks, continue 
          zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
          gvars$zenclick = length(zloc$x)
          K =  RPMG::whichbutt(iloc ,buttons)
          sloc = zloc
          
          
        }
      else
        {
###  Right button was clicked
          Nclick = 0
###  zenclick=zenclick+1
###   print(zenclick)
          K = 0
          gvars$zenclick = length(zloc$x)
          if(gvars$zenclick<1)
            {

           #######  No left mouse click was executed - stop and return to main 
 
              buttons = RPMG::rowBUTTONS(gvars$BLABS, col=rep(grey(.8), length(gvars$BLABS)),
                pch=rep("NULL", length(gvars$BLABS)))
              title("Done, Return to Calling Program")
              
              return(list(but="None", zloc=0, pix=0))
            }
         
        }
      
        if(K[Nclick] == match("OKAY", gvars$BLABS, nomatch = NOLAB))
          {

            if(length(zloc$x)>1)
              {
            w = gclick(gvars, zloc, YN)
            if(length(w)>0)
              {
                gvars$onoff =  rep(0 , Nopts)
                gvars$ocols =  rep(defaultcol , Nopts)
                gvars$onoff[w] =  1
                gvars$ocols[w] = addcol
              }
          }
            
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=rep(grey(.8), length(gvars$BLABS)), pch=rep("NULL", length(gvars$BLABS)))
            title("Return to Calling Program")
          

            retvec = gvars$opts[gvars$onoff==1]

            ## print(retvec)
            return(retvec)

            break;
            
          }
        
        if(K[Nclick] == match("DEFAULT", gvars$BLABS, nomatch = NOLAB))
          {


            
            kOPTS = unique(c(gvars$default, gvars$opts))
            
            kNopts = length(kOPTS)
            konoff =  rep(0 , kNopts )

            konoff[match(kOPTS , gvars$default)]  = 1
            
            kocols = rep(defaultcol , kNopts)
            
            kocols[konoff==1] = addcol

            
            gvars = list(BLABS=BLABS, opts=kOPTS, onoff=konoff, ocols=kocols, default=gvars$default )

            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
            
            sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
            
            
            
            Nclick=1
            K = 0
            zloc = list(x=NULL, y=NULL)
            gvars$zenclick = 0
            next
            
                 }

        if(K[Nclick] == match("ASORT", gvars$BLABS, nomatch = NOLAB))
          {


            KORD =   order(gvars$opts)
            kOPTS =   gvars$opts[KORD]
            
            
            kNopts = length(kOPTS)
   

            konoff = gvars$onoff[KORD]
            
            kocols = gvars$ocols[KORD]
            
            gvars = list(BLABS=BLABS, opts=kOPTS, onoff=konoff, ocols=kocols, default=gvars$default )

            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
            
            sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
            
            
            
            Nclick=1
            K = 0
            zloc = list(x=NULL, y=NULL)
            gvars$zenclick = 0
            next
            
                 }



         
        if(K[Nclick] == match("CANCEL", gvars$BLABS, nomatch = NOLAB))
          {
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=rep(grey(.8), length(gvars$BLABS)), pch=rep("NULL", length(gvars$BLABS)))
            title("Return to Calling Program")
            blist = zloc
            invisible(NULL)
            break;
            
          }


         if(K[Nclick] == match("SEL", gvars$BLABS, nomatch = NOLAB))
          {
            
            w = gclick(gvars, zloc, YN)
            gvars$onoff =  rep(0 , Nopts)
             gvars$ocols =  rep(defaultcol , Nopts)
            gvars$onoff[w] =  1
            gvars$ocols[w] = addcol
            
            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
              
              sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
              
              Nclick=1
              K = 0
              zloc = list(x=NULL, y=NULL)
              gvars$zenclick = 0
              
         next;
            
          }
         
        if(K[Nclick] == match("ADD", gvars$BLABS, nomatch = NOLAB))
          {

           w = gclick(gvars, zloc, YN)

            gvars$onoff[w] =  1
            gvars$ocols[w] = addcol
            
            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
              
              sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
              
              Nclick=1
              K = 0
              zloc = list(x=NULL, y=NULL)
              gvars$zenclick = 0
              
         next;
            
          }
        if(K[Nclick] == match("SUB", gvars$BLABS, nomatch = NOLAB))
          {


            w = gclick(gvars, zloc, YN)
            
            gvars$onoff[w] =  0
            gvars$ocols[w] = defaultcol
            
            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
            
            sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
              
            
            
              Nclick=1
              K = 0
              zloc = list(x=NULL, y=NULL)
              gvars$zenclick = 0
               next;

            
          }
        if(K[Nclick] == match("ALL", gvars$BLABS, nomatch = NOLAB))
          {

            gvars$onoff =  rep(1 , Nopts)
            gvars$ocols =  rep(addcol , Nopts)


            
            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
            
            sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
            
            
            
            Nclick=1
            K = 0
            zloc = list(x=NULL, y=NULL)
            gvars$zenclick = 0
            next
            
          }

        if(K[Nclick] == match("NONE", gvars$BLABS, nomatch = NOLAB))
          {


           
            gvars$onoff =  rep(0 , Nopts)
            gvars$ocols =  rep(defaultcol , Nopts)


            
            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
            
            sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
            
            
            
            Nclick=1
            K = 0
            zloc = list(x=NULL, y=NULL)
            gvars$zenclick = 0
            next
            
          }

        if(K[Nclick] == match("REVERT", gvars$BLABS, nomatch = NOLAB))
          {


            gvars$onoff =  ORIGonoff
            gvars$ocols =  ORIGcols
            
           
            YN = BUTREPLOT(gvars$opts, cols=gvars$ocols)
            buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
            
            sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
            
            
            
            Nclick=1
            K = 0
            zloc = list(x=NULL, y=NULL)
            gvars$zenclick = 0
            next
            
          }



      }


    retvec = gvars$opts[gvars$onoff==1]
    return(retvec)

  }


##########   source("SELBUT.R");   SELBUT(OPTS)


