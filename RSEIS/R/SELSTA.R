SELSTA<-function(GH, sel=1, newdev=TRUE, STAY=FALSE)
  {

    if(missing(newdev)) {  newdev=TRUE   }
    if(missing(STAY)) {  STAY=FALSE   }
    if(missing(sel)) {  sel=1:length(GH$STNS)  }

litecolors = c( "peachpuff2",      "darkolivegreen2", "slategray1" ,     "lightgoldenrod1",
  "darkseagreen3",   "lavenderblush2" , "slategray2" ,     "thistle1"      , 
  "cadetblue2"  ,    "lemonchiffon3"  )

    
  somecolors = c("black", "darkmagenta", "forestgreen", "blueviolet",
    "tan3", "lightseagreen", "deeppink", "cyan3", "bisque3",
    "magenta1", "lightsalmon3", "darkcyan", "darkslateblue",
    "chocolate4", "goldenrod4", "mediumseagreen")


    stdlab = c("OKAY", "QUIT", "REVERT", "ALLSTA", "ALLCOMP" , "SEL", "ADD", "SUB", "NONE", "ALL" )

    BLABS = c(stdlab)
    NLABS = length(BLABS)
    NOLAB = NLABS +1000

    colabs = rep(somecolors[2],length(BLABS))
    
    pchlabs = rep(4,length(BLABS))

    defaultcol = grey(.9)
    addcol = rgb(1, .9, .9)
    subcol = rgb(.9, .9, 1)


   
    ustas = unique(GH$STNS)
    ucomps = unique(GH$COMPS)

    stacols1 =  litecolors[3]
    stacols2 =  litecolors[1]

    compscols1 = litecolors[3]
    compscols2 = litecolors[1]

    staORIG = unique(GH$STNS[sel])
    compORIG = unique(GH$COMPS[sel])

     msta =     ustas  %in% staORIG 
     mcomp =    ucomps %in%  compORIG 

    scols1 = rep(stacols1, length(ustas))
    ccols1 = rep(compscols1, length(ucomps))

    
    sonoff = rep(0, length(ustas))
    sonoff[msta] = 1
    
    conoff = rep(0, length(ucomps))
    conoff[mcomp] = 1

    
    
    gvars = list(zenclick=0, ustas=ustas, NS=length(ustas), ucomps=ucomps, NC=length(ucomps), sonoff=sonoff, conoff=conoff, BLABS=BLABS )



    stareplot<-function(gvars)
      {
        plot(c(0,1), c(0,1) , type='n', axes=FALSE, xlab='', ylab='')
        scols1 =  rep(stacols1, length(gvars$ustas))
        ccols1 = rep(compscols1, length(gvars$ucomps))

      ##  print(gvars$sonoff)
        ##   print(gvars$conoff)

        
        scols1[gvars$sonoff==0] = stacols2
        ccols1[gvars$conoff==0] = compscols2 
        
        YN1 = BUTREPLOT(gvars$ustas, cols=scols1, ylim=c(.2, 1), newplot=FALSE)
         YN2 = BUTREPLOT(gvars$ucomps, cols=ccols1, ylim=c(0, .15), newplot=FALSE)
        buttons = RPMG::rowBUTTONS(gvars$BLABS, col=colabs, pch=pchlabs)
      
        return(list(YS=YN1, YC=YN2, buttons=buttons))
        
      }
    
    
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
            if(any(dim(YN$M$x)==1))
              {
                vx = as.vector(YN$M$x)
                vy = as.vector(YN$M$y)
              }
               else
                 {
                   vx = YN$M$x
                   vy = YN$M$y
                 }
                
            for(i in 1:klick)
              {
                w[i] = which.min( (vx+YN$dx/2 -thex[i])^2 +  (vy+YN$dy/2-they[i])^2)
              }
          }
        else
          {
            return(NULL)
          }
        return(w)
      }
    
    

    u = par("usr")
    sloc = list(x=c(u[1],u[2]), y=c(u[3],u[4]))
    zloc =list(x=NULL, y=NULL)
 
    Y = stareplot(gvars)

    
    while(TRUE)
      { ####### while loop
    ##    iloc = RPMG::ilocator(1, COL=rgb(1,0.6, 0.6), NUM=FALSE , YN=length(gvars$sel), style=-1)
        iloc = locator(1, type='p',  pch=21, cex=3,  col='red', bg='yellow')
        Nclick = length(iloc$x)
        
        if(Nclick>0)
          {
#######  add last click to list of clicks, continue 
            zloc  = list(x=c(zloc$x,iloc$x), y=c(zloc$y, iloc$y))
            gvars$zenclick = length(zloc$x)
            K =  RPMG::whichbutt(iloc , Y$buttons)
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
                Aselcomps = gvars$ucomps[gvars$conoff==1]
                Aselstas  = gvars$ustas[gvars$sonoff==1]
                
               ##    print(Aselcomps)
               ##    print(Aselstas )
                
                
                selp = which( GH$COMPS %in% Aselcomps & GH$STNS %in% Aselstas )
    
     
                buttons = RPMG::rowBUTTONS(gvars$BLABS, col=rep(grey(.8), length(gvars$BLABS)),
                  pch=rep("NULL", length(gvars$BLABS)))
                title("Done, Return to Calling Program")
                
                return(selp)
              }
            
          }
         if(K[Nclick] == match("QUIT", BLABS, nomatch = NOLAB))
          {
            zloc =list(x=NULL, y=NULL)
            buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)))
            title("Return to Calling Program")
            return(sel)
            break;
            
          }
        if(K[Nclick] == match("OKAY", BLABS, nomatch = NOLAB))
          {

            if(length(zloc$x)>1)
              {
                gvars$zenclick = length(zloc$x)
                w1 = gclick(gvars, zloc, Y$YS)
                 ##  print(w1)
                
                w2 = gclick(gvars, zloc, Y$YC)
                  ## print(w2)
                if(length(w1)>0)
                  {
                    gvars$sonoff[w1] =  1
                  }
                if(length(w2)>0)
                  {
                    gvars$conoff[w2] =  1
                  }
              }
              zloc =list(x=NULL, y=NULL)
            buttons = RPMG::rowBUTTONS(BLABS, col=rep(grey(.8), length(BLABS)), pch=rep("NULL", length(BLABS)))
            title("Return to Calling Program")
            break;
            
          }

        
        if(K[Nclick] == match("SEL", BLABS, nomatch = NOLAB))
          {

            if(length(zloc$x)>1)
              {
                gvars$zenclick = length(zloc$x)
                w1 = gclick(gvars, zloc, Y$YS)
                w2 = gclick(gvars, zloc, Y$YC)
                
                gvars$sonoff = rep(0, gvars$NS)
                gvars$conoff = rep(0, gvars$NC)
                 ##   print(w1)
                if(length(w1)>0)
                  {
                    gvars$sonoff[w1] =  1
                  }
               ##    print(w2)
                if(length(w2)>0)
                  {
                    gvars$conoff[w2] =  1
                  }
              }
            
            Y = stareplot(gvars)
              zloc =list(x=NULL, y=NULL)
            next;
            
          }
        if(K[Nclick] == match("ADD", BLABS, nomatch = NOLAB))
          {

            if(length(zloc$x)>1)
              {
                gvars$zenclick = length(zloc$x)
                w1 = gclick(gvars, zloc, Y$YS)
                w2 = gclick(gvars, zloc, Y$YC)
               
                ##    print(w1)
                if(length(w1)>0)
                  {
                    gvars$sonoff[w1] =  1
                  }
             ##      print(w2)
                if(length(w2)>0)
                  {
                    gvars$conoff[w2] =  1
                  }
              }
            
            Y = stareplot(gvars)
              zloc =list(x=NULL, y=NULL)
            next;
            
          }
        if(K[Nclick] == match("SUB", BLABS, nomatch = NOLAB))
          {

            if(length(zloc$x)>1)
              {
                gvars$zenclick = length(zloc$x)
                w1 = gclick(gvars, zloc, Y$YS)
                w2 = gclick(gvars, zloc, Y$YC)
               
                 ##   print(w1)
                if(length(w1)>0)
                  {
                    gvars$sonoff[w1] =  0
                  }
                ##   print(w2)
                if(length(w2)>0)
                  {
                    gvars$conoff[w2] =  0
                  }
              }
            
            Y = stareplot(gvars)
              zloc =list(x=NULL, y=NULL)
            next;
            
          }



        
        if(K[Nclick] == match("REVERT", BLABS, nomatch = NOLAB))
          {
                    gvars$sonoff =  sonoff
                    gvars$conoff =   conoff
                     Y = stareplot(gvars)
                      zloc =list(x=NULL, y=NULL)
            next;
          }
    
        
        if(K[Nclick] == match("NONE", BLABS, nomatch = NOLAB))
          {
                    gvars$sonoff =  rep(0, gvars$NS)
                    gvars$conoff =   rep(0, gvars$NC)
                     Y = stareplot(gvars)
                      zloc =list(x=NULL, y=NULL)
            next;
          }
        if(K[Nclick] == match("ALL", BLABS, nomatch = NOLAB))
          {
                    gvars$sonoff =  rep(1, gvars$NS)
                    gvars$conoff =   rep(1, gvars$NC)
                     Y = stareplot(gvars)
                      zloc =list(x=NULL, y=NULL)
            next;
          }
            
          if(K[Nclick] == match("ALLSTA", BLABS, nomatch = NOLAB))
          {
                    gvars$sonoff =  rep(1, gvars$NS)
                     Y = stareplot(gvars)
                      zloc =list(x=NULL, y=NULL)
            next;
          }
          if(K[Nclick] == match("ALLCOMP", BLABS, nomatch = NOLAB))
          {
                    gvars$conoff =  rep(1, gvars$NC)
                     Y = stareplot(gvars)
                      zloc =list(x=NULL, y=NULL)
            next;
          }
          

        
        
      }
    

    

    Aselcomps = gvars$ucomps[gvars$conoff==1]
    Aselstas  = gvars$ustas[gvars$sonoff==1]



    
    ##   print(Aselcomps)
    ##   print(Aselstas )
    
    
    selp = which( GH$COMPS %in% Aselcomps & GH$STNS %in% Aselstas )
    
    if(length(selp)<1)
      {

        selp = 1:length(GH$COMPS)
        

      }


    
    return(selp)

  }


########### source("/home/lees/selpgen.R");

####  g = selpgen(GH, newdev=FALSE, STAY=TRUE)

### cat(paste(GH$STNS[g], GH$COMPS[g]), sep="\n")





###########  selpgen(GH, newdev=FALSE, STAY=TRUE)

