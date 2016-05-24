COMPorder<-function(STNS, COMPS)
  {
    UNIsta = unique(STNS)
    
    oCOMPS = COMPS

    oCOMPS[grep("LD", COMPS)] = 1  
    oCOMPS[grep("I", COMPS)] = 1
    oCOMPS[grep("A", COMPS)] = 1
    
    oCOMPS[grep("V", COMPS)] = 2
    oCOMPS[grep("Z", COMPS)] = 2
    
    oCOMPS[grep("v", COMPS)] = 2
    oCOMPS[grep("z", COMPS)] = 2

    oCOMPS[grep("N", COMPS)] = 3
    oCOMPS[grep("n", COMPS)] = 3

    oCOMPS[grep("E", COMPS)] = 4
    oCOMPS[grep("e", COMPS)] = 4


    ords = match(STNS, sort(UNIsta))
    ordc = match(oCOMPS, sort(unique(oCOMPS) ))

    ordsel = order( ords+ordc/10)

    return(ordsel)

  }
