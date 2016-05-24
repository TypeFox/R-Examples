`PICK.DOC` <-
  function(w)
  {
    if(missing(w)) { w = NA }

    
    
   ALLLABS = c( 
"YPIX", "ROT.RT", "JustV", "JustE", "JustN", "WPIX", "NOPIX", 
"REPIX", "FILLPIX", "RIDPIX", "SEEPIX", "iNEXT", "PickWin", 
"pADDPIX", "Ppic", "Spic", "Apic", "POLSWITCH", "Pup", 
"Pnil", "Pdown", "NEXT", "PREV", "HALF", "CENTER", 
"MARK", "DOC", "REFRESH", "RESTORE", "ZOOM.out", "ZOOM.in", 
"LEFT", "RIGHT", "SCALE", "Xwin", "PSEL", "FLIP", 
"PTS", "FILT", "UNFILT", "fspread", "SPEC", "WWIN", 
"SGRAM", "WLET", "XTR", "Pinfo", "TSHIFT", "RMS", 
"LocStyle", "NA")
 

 

    
    N = length(ALLLABS)
    DOC = vector("list", N)
    
    names(DOC) = ALLLABS




   DOC[[1]] = 'generic (Y) pick' 
DOC[[2]] = 'Rotate seismogram to Radial Transverse' 
DOC[[3]] = 'Show only vertical'
DOC[[4]] = 'Show only East' 
DOC[[5]] = 'Show only North' 
DOC[[6]] = 'Window picks'
DOC[[7]] = 'Turn off Picks (set onoff to zero)' 
DOC[[8]] = 'Turn Picks back on' 
DOC[[9]] = 'Pick line spans vertical window'
DOC[[10]] = 'Remove picks' 
DOC[[11]] = 'print picks to screen'
DOC[[12]] = 'Do Nothing' 
DOC[[13]] = 'Spawn a 3-component picking window'
DOC[[14]] = 'add picks to list'
DOC[[15]] = 'P-wave pick' 
DOC[[16]] = 'S-wave pick' 
DOC[[17]] = 'Acoustic wave pick' 
DOC[[18]] = 'Switch Polarity'
DOC[[19]] = 'Up Polarity' 
DOC[[20]] = 'NUll Polarity'
DOC[[21]] = 'Down Polarity'
DOC[[22]] = 'Next BATCH of FILES'
DOC[[23]] = 'Previous BATCH of FILES'
DOC[[24]] = 'Shift Half a window'
DOC[[25]] = 'Center a window'
DOC[[26]] = 'Mark a trace' 
DOC[[27]] = 'Show documentation' 
DOC[[28]] = 'Refresh screen' 
DOC[[29]] = 'Restore from zoom' 
DOC[[30]] = 'Zoom out' 
DOC[[31]] = 'Zoom in' 
DOC[[32]] = 'Shift Left'
DOC[[33]] = 'Shift Right'
DOC[[34]] = 'Toggle Scale by trace/window'
DOC[[35]] = 'Delete all windows except main'
DOC[[36]] = 'Pick trace Sta/COMP to show' 
DOC[[37]] = 'Flip selected trace' 
DOC[[38]] = 'Show sample points' 
DOC[[39]] = 'Filter trace'
DOC[[40]] = 'Unfilter traces'
DOC[[41]] = 'do a filter spread on selection' 
DOC[[42]] = 'Display Spectrum' 
DOC[[43]] = 'Window' 
DOC[[44]] = 'Spectrogram' 
DOC[[45]] = 'Wavelet Transform'
DOC[[46]] = 'Extract single trace' 
DOC[[47]] = 'Pick information' 
DOC[[48]] = 'Shift traces to line up with first pick'
DOC[[49]] = 'Root Mean Square of selection'
DOC[[50]] = 'choose the locator style for picking in swig' 




J = 1:N
    
    if(is.character(w))
      {

        if(identical(w, "all") )
          {
            J = 1:N
          }
        else
          {
            J = match(w, ALLLABS)
          }


      }
    else
      {

        cp =  RPMG::chooser(ALLLABS, nsel=NA)
        
        ma = match(cp, ALLLABS)
        J = ma
        
      }
    
    if(length(J)<1) { print(DOC) }

    
    else
      {
        for(i in 1:length(J)){   cat(paste(ALLLABS[i], "=", DOC[[J[i]]])); cat("\n")  }
      }
#####    grep 'K\[Nclick\]\=\=match' PICK.R > ho.doc
#####  ho.doc = scan("/home/lees/Progs/R_stuff/ho.doc", what="", sep="\n")
#####   hot.doc = strsplit( ho.doc, split="\\\"")
#####   possiblebutts = unlist(lapply(hot.doc, getmem, 2))
    
#####      kli   = system("grep NOLAB /home/lees/R_PAX/RSEIS/R/swig.R | grep match", intern=TRUE)

#####     cp =  RPMG::chooser(blibs, nsel=NA)
#####     cp =  RPMG::chooser(blibs, nsel=0)

  }   


