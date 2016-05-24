`getvertsorder` <-
function(P, GU)
  {  ########  from pick file select firt arriving P-waves and return the order

    P = cleanpickfile(P)
    
    ustas = unique(P$STAS$name)
    verts  = vector()
    tees = as.numeric(P$STAS$sec)
    
    for(i in 1:length(ustas))
      {
        M = which(ustas[i] == P$STAS$name)
        w = which.min(P$STAS$sec[M])
        verts[i] = M[w]

      }

    #####  old way  verts = P$STAS$comp == "V"
    
    if(length(verts)<length(ustas))
      {
        sel=which(GU$COMPS=="V")
      }
    else
      {
        
        
        Vsta = P$STAS$name[verts]
        Vcomp = P$STAS$comp[verts]
        
        Varr = P$STAS$sec[verts]
        
        osta = Vsta[(order(Varr))]
        ocomp = Vcomp[(order(Varr))]
        
        
        msta = match(paste(sep=".", osta, ocomp), paste(sep=".", GU$STNS,GU$COMPS) )
        
        sel = msta
      }

    jd =getjul(P$LOC$yr, P$LOC$mo,P$LOC$dom )
    tadd =  secdif(GU$info$jd[msta], GU$info$hr[msta], GU$info$mi[msta], GU$info$sec[msta], jd, P$LOC$hr, P$LOC$mi, P$LOC$sec)
      
#### cat(sep="\n", paste(sep=' ',P$STAS$name, P$STAS$parr, tadd))
      win = c(max(min(tadd)-1, 0) , max(tadd)+1)
      


    return(list(sel=sel, win=win) )

  }

