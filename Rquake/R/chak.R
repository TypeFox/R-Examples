chak <-
function(DBnov , gstas, gcomps , fn, stas, vel, kind=2, Iendian=1, BIGLONG=FALSE)
  {
 
    buts = c("GPIX","PPIX", "PickWin", "fspread", "gMAP", "RQ" , "CONTPF")
   
    BUTLAB = c("REPLOT", "DONE",  
        "ZOOM.out", "ZOOM.in", "RESTORE", "Pinfo", 
         "SPEC", "SGRAM", "WLET", "FILT", "UNFILT", 
        "SCALE", "Postscript")

    
    load(fn)
    ####  for now swithc the Gpix to Ppix
    twpx$phase[twpx$phase=="G" |  twpx$phase=="Y" ] = "P"
   PF =  INITpickfile(stas=stas, src=NULL, WPX=twpx)

    PF$STAS$phase[ PF$STAS$phase=="Y" | PF$STAS$phase=="G" ] = "P"
   
   

  A1T = Qrangedatetime(twpx)
    at0 = A1T$min$jd+ A1T$min$hr/24+ A1T$min$mi/(24*60)+A1T$min$sec/(24*60*60)
    at1 = at0 -15/(24*3600) 
    at2 = at0 + 40/(24*3600)


   APAL = c("black", "darkmagenta", "forestgreen", "blueviolet", 
        "tan3", "lightseagreen", "deeppink", "cyan3", "bisque3", 
        "darkcyan", "darkred", "firebrick3", "rosybrown4", "royalblue3", 
        "darkblue", "red2", "violetred3", "springgreen3", "darkorange4", 
        "palevioletred3", "mediumpurple3", "tomato", "dodgerblue1", 
        "olivedrab4", "yellow4", "pink4")

     GH = RSEIS::Mine.seis(at1, at2, DBnov , gstas, gcomps, kind=2, Iendian=1, BIGLONG=FALSE)
###   add in the station information
    GH$sta = stas
    GH$pickfile = PF
    GH$vel = vel
    
    
###  assign colors to the traces
    
      
      stacols = list(name=stas$name, col=APAL[1:length(stas$name)])
      stas$col = stacols$col

        pcols = stacols$col[match(GH$STNS, stas$name)]

        GH$pcol = pcols

    onlyPY = twpx[twpx$phase=="P" | twpx$phase=="Y", ]
    
      s1 = RSEIS::secdifL(A1T$min,   onlyPY )
        
        ords1 = order(s1)
        
        osta =   onlyPY$name[ords1]


ohoh = list(dist=s1[ords1] , name=osta)
            
           jord =  RSEIS::seisorder(GH, ohoh, VNE="V")
            
  hret = RSEIS::swig(GH, sel=jord, APIX=twpx,  STDLAB=BUTLAB,  PADDLAB=buts)

    
return(hret)
    
  }
