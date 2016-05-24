selpgen<-function(MH, newdev=TRUE, STAY=FALSE)
  {

    if(missing(newdev)) {  newdev=TRUE   }
    if(missing(STAY)) {  STAY=FALSE   }

    ustas = unique(MH$STNS)
    ucomps = unique(MH$COMPS)


    pastcols = RPMG::pastel.colors(10, seed = 1)


    ns = length(ustas)
    nc = length(ucomps)


    RCOL = c(ucomps, "ALL STAS", "ALL COMPS", "ALL")
 
    nrs  = length(RCOL)


   
    
    ncol = 5
    
    if(ns<=nrs){ ncol = 2 }


   
    
    N = ns
    nrow1 = round((N/(ncol-1)) + 0.5)



    NROW = max(c(nrow1, nrs))

    Ntot = NROW *  ncol
    
    
    blnkrs = 0
    if(NROW>nrs)
      {

        blnkrs = (NROW-nrs)
        RCOL = c(RCOL, rep(" ", times=blnkrs))
      }

    

    LEFTover =  NROW*(ncol-1)-ns
    
  
    if(LEFTover>0)
      {

        ustas = c(ustas, rep(" ", times=LEFTover))

      }
    
    PP = c(ustas, RCOL )

    NP = length(PP)
    
    
    cols = c( rep(pastcols[1], times=length(PP)))
    cols[PP=="ALL STAS"] = pastcols[2]
    cols[PP=="ALL COMPS"] = pastcols[3]
    cols[PP=="ALL"] = pastcols[4]
    cols[match(ucomps,  PP  ) ] = pastcols[5]
    cols[match(ustas,  PP  ) ] = pastcols[6]

      
      P2 = RPMG::chooser( PP, ncol=5, nsel=NA, cols =cols, newdev=newdev, STAY=STAY, main="" , pch=21, cex=3,  col='red', bg='blue' )

    selp = 0
    
    if(any("ALL"==P2))
      {
        selp= 1:length(MH$STNS)
        return(selp)
      }

    if(any("ALL STAS"==P2))
      {
        Aselstas = ustas
        
      }
    else
      {
        Aselstas =  P2[P2 %in%  MH$STNS]
        
         if(length(Aselstas)<1)  Aselstas = ucomps

      }

    

    if(any("ALL COMPS"==P2))
      {
        Aselcomps = ucomps
       
      }
    else
      {
        Aselcomps =  P2[P2 %in%  MH$COMPS]
        if(length(Aselcomps)<1)  Aselcomps = ucomps
      }

    

     
    selp = which( MH$COMPS %in% Aselcomps & MH$STNS %in% Aselstas )
    

    return(selp)

    }


########### source("/home/lees/selpgen.R");

####  g = selpgen(GH, newdev=FALSE, STAY=TRUE)

### cat(paste(GH$STNS[g], GH$COMPS[g]), sep="\n")





###########  selpgen(GH, newdev=FALSE, STAY=TRUE)

