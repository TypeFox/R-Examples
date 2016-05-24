`GET.seis` <-
function(fnames, kind=1, Iendian=1, BIGLONG=FALSE, HEADONLY=FALSE , PLOT=-1, RAW=FALSE)
{
  ###  get a list of SEGY or SAC files from a directory and store in structure
  ####  kind 1=segy, 2=sac, kind=0->R format  ( 3=AH is no longer available)
  
  if(missing(PLOT)) { PLOT=-1 }
  if(missing(kind)) { kind=1 }
  if(missing(HEADONLY)) {HEADONLY=FALSE }
  if(missing(Iendian)) { Iendian=1 }
  if(missing(BIGLONG)) { BIGLONG=FALSE}
  if(missing(RAW)) { RAW=FALSE }

  
  tmpGIVE = as.list(1:length(fnames))

  ii = 1

  DATIM =  rep(0,length=4)
  n=1
  dt=0.025000
  sec = 0
  thesta="XXXXX"
  thecomp="XXXXX"


  for(i in 1:length(fnames))
    {

      fn = fnames[i]
      infile = fn
    ####   print(fn);
      ###  if this file does not exist, exit!
      if(file.exists(infile)==FALSE)
        {
         print(paste(sep=' ', "file does not exist", fn) ); 
          next;
        }
      else
        {
          ###  print(paste(sep=' ', "file exists", fn) );

        }

##################################  if the file is already in R format
      ###############  just load it it and skip to next file
      ##  this code assumes that the list is called DAT
      if(kind==0)
        {
          DAT = list()
          GED  = load(fn)
          assign("DAT", get(GED))
          if(is.null(names(DAT)) & length(DAT)>=1)
            {
              DAT = DAT[[1]]

            }
          
          if(HEADONLY) DAT$amp = NULL
          DAT$oldname = DAT$fn
          DAT$fn = fn
          tmpGIVE[[i]] = DAT

          
          next
        }

      if(kind==1)
        {
        #####   print(paste("RAW=", RAW))
          
          DAT  = segy2rseis(fn, Iendian=Iendian, HEADONLY=HEADONLY , BIGLONG=BIGLONG, PLOT=PLOT, RAW=RAW)
          tmpGIVE[[i]] = DAT[[1]]
          next
        }

      if(kind==2)
        {
          DAT  = sac2rseis(fn, Iendian=Iendian, HEADONLY=HEADONLY , BIGLONG=BIGLONG, PLOT=PLOT, RAW=RAW)
          tmpGIVE[[i]] = DAT[[1]]
          next
        }

      if(kind==3)
        {
          print("AH format currently not available")
          tmpGIVE[[i]] = NA
          next
        }

      #############  if no kind matches, better return NA
      tmpGIVE[[i]] = NA

      
    }
  invisible(tmpGIVE)
}

