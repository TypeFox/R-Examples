`local.file` <-
function(pref, suf)
  {
##X##     ###  get a file name that is in the local directory
##X##     ### and that is free, i.e. a new file name
##X##     ###  used to avoid writing over files
    ###  e.g.:  plfname = local.file("test","eps")
    i = 0
    indchar       = formatC(i, format="d", width=4,  flag="0")
    tfile = paste(sep='_', pref, indchar)    
    ofile = paste(sep='.', tfile, suf)
 
   ###  vof = system(paste(sep=" ", "ls", ofile), intern=TRUE, ignore.stderr=TRUE)

    vof = file.exists(ofile)
    
    i = 0
    while(vof==TRUE)
      {

        indchar       = formatC(i, format="d", width=4,  flag="0")
        tfile = paste(sep='_', pref, indchar)    
        ofile = paste(sep='.', tfile, suf)
 
       ###  ofile = paste(sep='_', pref, indchar, suf)
 
       ###  vof = system(paste(sep=" ", "ls", ofile), intern=TRUE, ignore.stderr=TRUE)
        vof = file.exists(ofile)
        
        i = i +1
      }
    return(ofile)

  }

