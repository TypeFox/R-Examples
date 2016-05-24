SELOPT<-function(OPTS, onoff=-1, ncol=5, ocols ="white", cex=1, default="opt"  )
  {

##    
    Nopts = length(OPTS)
    ORIGsel  = NULL
    
    if(missing(onoff)) { onoff =  rep(0 , Nopts )  }
    if(missing(default)) { default=OPTS  }
    if(missing(cex)) {  cex=1  }


    
    if(length(onoff)==1)
      {
        
        if(onoff==0) {
          onoff = rep(1, times=Nopts)
           ORIGsel  = 1:Nopts
        }

      }

    if(length(onoff)==1)
      {
        
        if(onoff<0) {
          onoff = rep(0, times=Nopts)
          ORIGsel  = NULL
        }

      }
    if(length(onoff)<Nopts)
      {
    
        aonoff = rep(0, times=Nopts)
        aonoff[onoff] = 1
        ORIGsel  =onoff    
        onoff = aonoff
      }

    

    
    defaultcol = grey(.9)
    addcol = rgb(1, .9, .9)
    subcol = rgb(.9, .9, 1)
    
    
    stdlab = c("OKAY", "CANCEL", "DEFAULT", "REVERT", "ASORT", "SEL", "ADD", "SUB", "NONE", "ALL", "PRINT" )
    
    BLABS = c(stdlab)
    NLABS = length(BLABS)
    NOLAB = NLABS +1000
    
    colabs = rep("black",length(BLABS))
    
    pchlabs = rep(4,length(BLABS))
    
    
    if(missing(ocols))
      {
        ocols = rep(defaultcol , Nopts)
        ocols[onoff==1] = addcol
      }

    
        
    ORIGonoff =onoff
    ORIGcols = ocols
    

    DOREPLOT<-function(gvars)
    {
      
      Y1 = OPTREPLOT(gvars$opts, cols=gvars$ocols, sel=which(gvars$onoff==1  ),
        ncol=gvars$ncol, cex=gvars$cex, slwd=3, blwd=3,  bcol="transparent", mpct=0.05 )
      
      return(Y1)
    }

    
    
    u = par("usr")
    sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
    zloc =list(x=NULL, y=NULL)
    
    gvars = list(BLABS=BLABS, opts=OPTS, onoff=onoff, ocols=ocols,
      default=default,  ncol=ncol , sel=NULL, cex=cex, worder =ORIGsel , norder=length(ORIGsel) )
    
    gclick<-function(gvars, zloc, YN, one=FALSE)
      {

        if(missing(one))
          {  one=FALSE }
        
        if(one) { klick = 1 } else{   klick = gvars$zenclick-1 }
        
        if(klick<1) return(NULL)
        
        thex = zloc$x[1:klick]
        they = zloc$y[1:klick]

        flag = thex>=YN$rx[1] & thex<=YN$rx[2] & they>=YN$ry[1] & they<=YN$ry[2]
        if(all(!flag)) { return(NULL) }

        
        thex = thex[flag]
        they = they[flag]
        
        klick = length(thex)
        if(klick<1)  return(NULL)
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


    YN = DOREPLOT(gvars)
    
    buttons = rowBUTTONS(BLABS, col=colabs, pch=pchlabs)


    
    
    while(TRUE)
      {
    ##    iloc = ilocator(1, COL=rgb(1,0.6, 0.6), NUM=FALSE , YN=length(gvars$sel), style=-1)
        iloc = locator(1, type='p',  pch=21, cex=3,  col='red', bg='yellow')
        Nclick = length(iloc$x)
        
        if(Nclick>0)
          {
#######  add last click to list of clicks, continue

                w = gclick(gvars, iloc, YN, one=TRUE)
                if(!is.null(w))
                  {
                    rect(YN$M$x[w], YN$M$y[w], YN$M$x[w]+YN$dx, YN$M$y[w]+YN$dy, border="black", lwd=3 )
                    
                    gvars$worder  = unique(c(gvars$worder, w) )
                    gvars$norder = length(gvars$worder)
                    
                  }

            zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
            gvars$zenclick = length(zloc$x)
            K =  whichbutt(iloc ,buttons)
            sloc = zloc
            
            
          }
        else
          {
###  Right button was clicked
            Nclick = 0
###  zenclick=zenclick+1
            print("right button clicked")
            K = 0
            gvars$zenclick = length(zloc$x)
            if(gvars$zenclick<1)
              {
                
#######  No left mouse click was executed - stop and return to main 
                
                buttons = rowBUTTONS(gvars$BLABS, col=rep(grey(.8), length(gvars$BLABS)),
                pch=rep("NULL", length(gvars$BLABS)))
                title("Done, Return to Calling Program")
                
                return(NULL)
                break
              }
            else
              {
                w = gclick(gvars, zloc, YN, one=TRUE)
                if(!is.null(w))
                  {
                    gvars$onoff[w] =  1
                    gvars$worder = c(gvars$worder,  w)
                    gvars$worder =  unique(gvars$worder)
                    gvars$norder = length(gvars$worder)
                    
                  }
                ##   gvars$ocols[w] = addcol
                
                
                
                print(gvars$worder)
                print(gvars$opts[gvars$worder])
                
                retvec = gvars$opts[gvars$worder]
                
                ## print(OKAY)
                return(retvec)
                break


              }
            break
            
        }
        
        if(K[Nclick] == match("OKAY", gvars$BLABS, nomatch = NOLAB))
          {
            
            if(length(zloc$x)>1)
              {
                w = gclick(gvars, zloc, YN)
                if(length(w)>0)
                  {
                    gvars$onoff =  rep(0 , Nopts)

                    gvars$worder  = unique(c(gvars$worder, w) )
                    gvars$norder = length(gvars$worder)
                    
                  ##  gvars$ocols =  rep(defaultcol , Nopts)
                 ##   gvars$onoff[w] =  1
                 ##   gvars$ocols[w] = addcol
                  }
              }
            
            buttons = rowBUTTONS(gvars$BLABS, col=rep(grey(.8), length(gvars$BLABS)), pch=rep("NULL", length(gvars$BLABS)))
            title("Return to Calling Program")
          

            ## print(gvars$worder)
           ##  print(gvars$opts[gvars$worder])
            
           ##  retvec = gvars$opts[gvars$onoff==1]
            retvec = gvars$opts[gvars$worder]
                
            ## print(OKAY)
            return(retvec)

            break;
            
          }
        
        if(K[Nclick] == match("DEFAULT", gvars$BLABS, nomatch = NOLAB))
          {


            
            kOPTS = unique(c(gvars$default, gvars$opts))
            
            kNopts = length(kOPTS)
            konoff =  rep(0 , kNopts )

            konoff[match(kOPTS , gvars$default)]  = 1
            kocols = ORIGcols
            ##  kocols = rep(defaultcol , kNopts)
            
           ##  kocols[konoff==1] = addcol

            
            gvars = list(BLABS=BLABS, opts=kOPTS, onoff=konoff, ocols=kocols,
              default=gvars$default, ncol=ncol, cex=cex, worder =ORIGsel , norder=length(ORIGsel) )

            YN = DOREPLOT(gvars)
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
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
            YN = DOREPLOT(gvars)
            
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
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
            buttons = rowBUTTONS(gvars$BLABS, col=rep(grey(.8), length(gvars$BLABS)), pch=rep("NULL", length(gvars$BLABS)))
            title("Return to Calling Program")
            blist = zloc
            return(NULL)
            break;
            
          }


         if(K[Nclick] == match("SEL", gvars$BLABS, nomatch = NOLAB))
          {
            
            w = gclick(gvars, zloc, YN)
            gvars$onoff =  rep(0 , Nopts)
          ##    gvars$ocols =  rep(defaultcol , Nopts)
            gvars$onoff[w] =  1
         ##    gvars$ocols[w] = addcol
            gvars$worder = w
            gvars$norder = length(gvars$worder)
            YN =DOREPLOT(gvars)
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
              
              sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
              
              Nclick=1
              K = 0
              zloc = list(x=NULL, y=NULL)
              gvars$zenclick = 0
              
         next;
            
          }

         if(K[Nclick] == match("PRINT", gvars$BLABS, nomatch = NOLAB))
          {
            print(gvars$worder)
            CHOICE = gvars$opts[gvars$worder]
           ##     print(gvars$worder)
  ## 
               b = as.character(paste(sQuote(CHOICE), collapse = ","))

            
             cat(file = "", paste(sep = "", "CHOICE" , "=c(", b, ")"), fill = TRUE)


            
              Nclick=1
              K = 0
              zloc = list(x=NULL, y=NULL)
              gvars$zenclick = 0
              
         next;
            
          }
         


        
        if(K[Nclick] == match("ADD", gvars$BLABS, nomatch = NOLAB))
          {

           w = gclick(gvars, zloc, YN)
           
           if(!is.null(w))
             {
               gvars$onoff[w] =  1
               gvars$worder = c(gvars$worder,  w)
               gvars$worder =  unique(gvars$worder)
               gvars$norder = length(gvars$worder)
               
             }
          ##   gvars$ocols[w] = addcol

            YN =DOREPLOT(gvars)
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
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
            if(!is.null(w))
              {
                gvars$onoff[w] =  0
                mw = match( w, gvars$worder )

              ##   print(w)
               ##  print(mw)
               ##  print(gvars$worder)
             ##    print(gvars$worder[ mw ] )
                
                if(length(mw)>0)
                  {
                    gvars$worder = gvars$worder[-mw]
                    gvars$norder = length(gvars$worder)
                  }

              }
          ##   gvars$ocols[w] = defaultcol
            
             YN =DOREPLOT(gvars)
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
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
          ##   gvars$ocols =  rep(addcol , Nopts)

               gvars$worder =  seq(from=1, to= Nopts)
               gvars$norder = length(gvars$worder)
            
             YN =DOREPLOT(gvars)
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
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
          ##   gvars$ocols =  rep(defaultcol , Nopts)

               gvars$worder =  vector()
               gvars$norder = 0
 
             YN =DOREPLOT(gvars)
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
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
            gvars$worder =  ORIGsel
            gvars$norder =  length(ORIGsel)
 
           
             YN =DOREPLOT(gvars)
            buttons = rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
            
            u = par("usr")
            
            sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
            
            Nclick=1
            K = 0
            zloc = list(x=NULL, y=NULL)
            gvars$zenclick = 0
            next
            
          }



      }

    print("BOTTOM of CODE")
    retvec = gvars$opts[gvars$onoff==1]
    return(retvec)

  }


##########   

