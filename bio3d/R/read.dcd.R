"read.dcd" <-
function(trjfile, big=FALSE, verbose=TRUE, cell = FALSE){

  # Version 0.2 ... Tue Jan 18 14:20:12 PST 2011
  # Version 0.1 ... Thu Mar  9 21:18:54 PST 2005
  #
  # Description:
  #  Reads a CHARMM or X-PLOR/NAMD binary
  #  trajectory file with either big- or
  #  little-endian storage formats
  #
  # Details:
  #  Reading is accomplished with two different
  #  functions.
  #   1. 'dcd.header' which reads headder info
  #   2. 'dcd.frame' takes the header info and
  #         reads frame by frame producing a 
  #         nframes/natom*3 matrix of cartisean
  #         coordinates



#===DCD=FORMAT==============================================
#HDR  NSET  ISTRT NSAVC 5-ZEROS NATOM-NFREAT DELTA   9-ZEROS
#CORD files step1 step  zeroes  (zero)      timestep  zeroes
#C*4  INT   INT   INT    5INT    INT         DOUBLE   9INT
#  [CHARACTER*20]
#===========================================================
#NTITLE   TITLE
#INT      C*MAXTITL
#C*2      C*80
#===========================================================
#NATOM
#INT
#===========================================================
#CELL(I), I=1,6          (DOUBLE)
#===========================================================
#X(I), I=1,NATOM         (SINGLE)
#Y(I), I=1,NATOM         
#Z(I), I=1,NATOM         
#===========================================================


  dcd.header <- function(trj,...) {
    
    # Read DCD Header section
    
    end = .Platform$endian   # Check endianism
    check <- readBin(trj,"integer",1,endian=end)
    # first thing in file should be an '84' header
    if (check != 84) {
      # if not we have the wrong endianism
      if (end == "little") { end="big" } else { end="little" }
      check <- readBin(writeBin(check, raw()), "integer", 1, endian = end)

      if (check != 84) {
        close(trj)
        stop("PROBLEM with endian detection")
      }
    }

    hdr <- readChar(trj, nchars=4) # data => CORD or VELD

    # how big is the file 'end.pos'
    cur.pos <- seek(trj, where=1, origin = "end") # pos ?
    end.pos <- seek(trj, where=cur.pos, origin= "start")
  
    icntrl <- readBin(trj,"integer", 20, endian=end) # data => header info

    # header information:
    nframe = icntrl[1]  # number of frames
    first  = icntrl[2]  # number of previous steps 
    step   = icntrl[3]  # frequency of saving
    nstep  = icntrl[4]  # total number of steps
    nfile <- nstep/step  # number of files
    last  <- first + (step * nframe) # last step
    # 5 zeros
    ndegf  = icntrl[8]  # number of degrees of freedom
    nfixed = icntrl[9]  # number of fixed atoms 
    delta  = icntrl[10] # coded time step
    cryst  = icntrl[11] # crystallographic group
    block  = icntrl[12] # extra block?
    # 9 zeros
    vers     = icntrl[20]

    # flush to end of line
    a<-readBin(trj,"integer",1, endian=end) # should be '84' line tail
  ## cur.pos<-seek(trj, where=92, origin= "start") # position 92
    rm(icntrl)  # tidy up

  
    # Are we CHARMM or X-PLOR format
    charmm=FALSE; extrablock=FALSE; four.dims=FALSE
    if (vers != 0) {
      charmm=TRUE # charmm version number
      if (cryst == 1) { # check for
        extrablock = TRUE  # extra free 
      }                 # atom block &
      if (block == 1) { # extra four
        four.dims=TRUE     # dimensions 
      }
    } else {
      # re-read X-PLOR delta as a double
      cur.pos <- seek(trj, where=44, origin= "start")
      delta  = readBin(trj,"double", 1, endian=end)
      seek(trj, where=cur.pos, origin= "start")
    }
  

    #=======#
    # Title #
    a<-readBin(trj,"integer",1, endian=end)        # flush FORTRAN header 
    ntitle <- readBin(trj,"integer",1, endian=end) # data => Num title lines
    title<-NULL                                    # store title & date

    cur.pos <- seek(trj, where=NA)  ## 100
    for (i in 1:ntitle) {
### ==> !!!nasty hack due to invalid UTF-8 input (Jun 5th 07) !!! <=== ###    
      ll<-try(title<-c( title, readChar(trj,80) ),silent=TRUE)
    }
    # OR: title<- readChar(trj, (ntitle*80))

    if(class(ll)=="try-error") {
      warning("Check DCD header data is correct, particulary natom")
      ##cur.pos <- seek(trj, where=260, origin= "start") # pos 260
      cur.pos <- seek(trj, where=(80*ntitle+cur.pos), origin= "start")
    }
### == end hack
    a<-readBin(trj,"integer",1, endian=end)         # flush FORTRAN tail

    
    #=======#
    # Natom #
    a<-readBin(trj,"integer",1, endian=end)        # flush FORTRAN header 
    natom <- readBin(trj,"integer",1, endian=end)  # number of atoms
    a<-readBin(trj,"integer",1, endian=end)        # flush FORTRAN tail

  ##cur.pos <- seek(trj, where=276, origin= "start") # pos 276


    #=============#
    # Freeindexes #
    if (nfixed != 0) {
      # Free (movable) atom indexes if nfixed > 0
      a <- readBin(trj,"integer",1, endian=end)        # flush FORTRAN header 
      free.ind <- readBin(trj,"integer", (natom-nfixed), endian=end )
      a <- readBin(trj,"integer",1, endian=end)        # flush FORTRAN tail
      print("FIXED ATOMS IN SIMULATION => CAN'T READ YET")
    }
  

  
    if (verbose) {
## EDIT ## R version 2.11.0 does not like "\0", just remove for now - Apr 12 2010   
##      cat( sub(" +$","",gsub(pattern="\0", replacement="", x=title)),sep="\n" )
      cat(" NATOM =",natom,"\n")
      cat(" NFRAME=",nframe,"\n")
      cat(" ISTART=",first,"\n")
      cat(" last  =",last,"\n")
      cat(" nstep =",nstep,"\n")
      cat(" nfile =",nfile,"\n")
      cat(" NSAVE =",step,"\n")
      cat(" NDEGF =",ndegf,"\n")
      cat(" version",vers,"\n")
    }

    
    # Done with Header :-)
    header <- list(natom=natom,
                   nframe=nframe,
                   first=first,
                   last=last,
                   nstep=nstep,
                   nfile=nfile,
                   step=step,
                   ndegf=ndegf,
                   nfixed=nfixed,
                   charmm=charmm,
                   extrablock=extrablock,
                   four.dims=four.dims,
                   end.pos=end.pos,
                   end=end)

  }


  
  dcd.frame <- function(trj, head, cell) {
    
    # DCD step/frame data
    #  read one frame from the current conection 'trj'
    #  which should have been already through
    #  'dcd.header' so the "where" position is at
    #  the start of the cooedinate section

    #============#
    # Free atoms #
    # Uncomment the next two lines if reading cell
    # parameters only works with CHARMM DCD files
#     if(!head$charmm && cell)
#       stop("Cell parameters can only be read from CHARMM dcd files.")
    if ( head$charmm &&  head$extrablock) {
      # CHARMM files may contain lattice parameters
      a <- readBin(trj,"integer",1, endian=head$end)  # flush FORTRAN header 
      u <- readBin(trj, "numeric", size = 8, n = (a/8),endian = head$end)
      a <- readBin(trj,"integer",1, endian=head$end)  # flush FORTRAN tail
    }

    ##cur.pos <- seek(trj, where=332, origin= "start") # pos 332
  
    #========#
    # Coords #
    if (head$nfixed == 0) {
      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN header 
      x <- readBin(trj,"numeric",                    # read x coords
                   size=4, n=(a/4), endian=head$end)
      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN tail

      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN header 
      y <- readBin(trj,"numeric",       # read y coords
                   size=4, n=(a/4), endian=head$end)
      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN tail
    
      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN header 
      z <- readBin(trj,"numeric",       # read z coords
                   size=4, n=(a/4), endian=head$end)
      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN tail

    } else {
      # not implemented yet! => cant cope with fixed atoms
    }

    #===============#
    # 4th dimension #
    if (head$charmm && head$four.dims) {
      # CHARMM files may contain an extra block?
      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN header 
      seek(trj, where=a, origin= "current") # skip this block 
      a <- readBin(trj,"integer",1, endian=head$end) # flush FORTRAN tail    
    }

    # Done with coord frame :-)
    
    #coords <- list(x=x,
    #               y=y,
    #               z=z)

    if(cell) to.return <- c( u[c(1,3,6)], (180/pi)*acos(u[c(5,4,2)]))
    else to.return <- as.vector(rbind(x,y,z))
    class(to.return) = "xyz"
    return(to.return)
  }
  
  # Check if file exists
  if( !file.exists(trjfile) ) {
    stop(paste("No input DCD file found with name:", trjfile))
  }
  
  # Open file conection
  trj <- file(trjfile, "rb")
  #verbose=T
  head<-dcd.header(trj,verbose)

  nframes = head$nframe 
  natoms  = head$natom

  # blank xyz data structures
  # format: rows => nframes, cols => natoms
### ==> !!! Insert to read big dcd files (Sep 29th 08) !!! <=== ###    
  ###xyz <- matrix(NA, nrow=nframes,ncol=natoms*3)
  if(!big) {
    if(cell) to.return <- matrix(NA, nrow=nframes,ncol=6)
    else to.return <- matrix(NA, nrow=nframes,ncol=natoms*3)
  } else {
    ##-! Insert to read big dcd files (Sep 29th 08)
    oops <- requireNamespace("bigmemory", quietly = TRUE)  
    if(!oops) 
        stop("Please install the bigmemory package from CRAN")

    if(cell) to.return <- bigmemory::big.matrix(nrow=nframes,ncol=6, init = NA, type = "double")
    else to.return <- bigmemory::big.matrix(nrow=nframes,ncol=natoms*3, init = NA, type = "double")
  }
### ==> !!! end big.matrix insert
  
  if(verbose){ cat("Reading (x100)") }
  store <- NULL
  
  # fill xyz with frame coords
  if(verbose) pb <- txtProgressBar(1, nframes, style=3)
  for(i in 1:nframes) {
    curr.pos <- seek(trj, where=0, origin= "current")
    if (curr.pos <= head$end.pos) {
      to.return[i,]<-dcd.frame(trj,head,cell)
      if (verbose) {
        setTxtProgressBar(pb, i)
#        if(i %% 100==0) { cat(".") }
      }
#      print(paste("frame:",i,"pos:",curr.pos))
      store<-cbind(store,curr.pos)
    } else {
      print("Premature end of file")
      print(paste("  last frame:",i,
                  "nframe:",head$nframe ))
      break
    }
  }
#  if(verbose) { cat("done",sep="\n") }
  if(verbose)
    cat("\n")
  close(trj)

  ##class(to.return) = "xyz"
  return( as.xyz(to.return) ) 
}

