

########################################################################
#                        get.ncbi
#
# Try to connect to ncbi to get a list of complete bacterial genomes.
# Returns an n by 5 dataframe with 
#   species names
#   accesion number
#   size in bp
#   type (chromosome or plasmid)
#   Last update time
#
#########################################################################
get.ncbi <- function(repository = "ftp://ftp.ncbi.nih.gov/genomes/Bacteria/"  )
{
  #
  # First of all, check that this computer is not off the net:
  #
  if( ! capabilities("http/ftp") )
    stop("capabilities(\"http/ftp\") is not TRUE") 

#  
# BEGIN Proxy problem
#
# I have a problem here: the ftp connection apparently does
# not work when there is a proxy. I have fixed the bug
# this way, but this is a rather crude and unsatisfactory
# solution.
#
  ftp.proxy.bck <- Sys.getenv("ftp_proxy")

  if( ftp.proxy.bck != "" ) # there is a proxy
  {
    warning("I'am trying to neutralize proxies")
    Sys.setenv("no_proxy" = "") 
  }
#
# END Proxy problem
#  

  #
  # Try to get list of folder in ncbi repository. Note that R build-in ftp
  # does not allow to do this directly, so we rely on a system call that
  # will run only under Unix systems. Moreover, not all ftp client supports
  # this syntax. 
  #
  sysinfo <- Sys.info()[1]

  if( sysinfo == "Darwin" )
  {
    cmd <- sprintf("echo \"ls\" | ftp %s", repository)
    brut <- readLines(pipe(cmd))
  }
  else if( sysinfo == "SunOS" )
  {
    #
    # Build command file for ftp connection
    #
    tmpname <- tempfile(pattern="getncbi")
    tmpcmdfile <- file(tmpname, open="w")
    writeLines("user anonymous seqteam@biomserv.univ-lyon1.fr", tmpcmdfile)
    writeLines("cd genomes/Bacteria", tmpcmdfile)
    writeLines("dir", tmpcmdfile)
    writeLines("bye", tmpcmdfile)
    close(tmpcmdfile)
    #
    hostname <- unlist(strsplit(repository,split="/"))[3]
    cmd <- sprintf("ftp -v -n %s < %s", hostname, tmpname)
    brut <- readLines(pipe(cmd))
  }
  else if( sysinfo == "Linux" ){
    #
    # Build command file for ftp connection
    #
    tmpname <- tempfile(pattern="getncbi")
    tmpcmdfile <- file(tmpname, open="w")
    writeLines("user anonymous seqteam@biomserv.univ-lyon1.fr", tmpcmdfile)
    writeLines("cd genomes/Bacteria", tmpcmdfile)
    writeLines("dir", tmpcmdfile)
    writeLines("bye", tmpcmdfile)
    close(tmpcmdfile)
    #
    hostname <- unlist(strsplit(repository,split="/"))[3]
    cmd <- sprintf("ftp -v -n %s < %s", hostname, tmpname)
    brut <- readLines(pipe(cmd))

  }
  else
  {
    stop("Unimplemented platform")
  } 

  #
  # Keep only lines corresponding to folders:
  #
  brut <- brut[grep("dr-xr-xr-x", brut)]


  #
  # Now there should be a vector of chr in "brut", each line looking like:
  #
  # "dr-xr-xr-x   2 ftp      anonymous     4096 Jul 15 15:41 Aeropyrum_pernix"
  #

  brut <- sapply( brut, strsplit, split=" ")
  names(brut) <- NULL

  #
  # Now each element in "brut" should be splited as in:
  #
  # [1] "dr-xr-xr-x"       ""                 ""                 "2"               
  # [5] "ftp"              ""                 ""                 ""                
  # [9] ""                 ""                 "anonymous"        ""                
  # [13] ""                 ""                 ""                 "4096"            
  # [17] "Jul"              "15"               "15:41"            "Aeropyrum_pernix"
  #

  get.last <- function( vector )
  {
    return( vector[length(vector)] )
  }
  brut <- sapply( brut, get.last)
  if(length(grep("CLUSTERS",brut))!=0){
    brut<-brut[-grep("CLUSTERS",brut)]  # we remove the CLUSTERS folder, since it doesn't contain any annotated bacterial genomes
  }
  
  #
  # Now "brut" should contains folders names as in:
  # > brut[1:5]
  # [1] "Aeropyrum_pernix"                     
  # [2] "Agrobacterium_tumefaciens_C58_Cereon" 
  # [3] "Agrobacterium_tumefaciens_C58_UWash"  
  # [4] "Aquifex_aeolicus"                     
  # [5] "Archaeoglobus_fulgidus"               
  #

  #
  # Set vector types for results:
  #
  
  species <- character(0)
  accession <- character(0)
  size.bp <- integer(0)
  type <- character(0)
  lastupdate <- character(0)

  #
  # Main loop on folders to see what's inside
  #
  for( folder in brut )
  {
    if( sysinfo == "Darwin" )
    {
      where <- paste(repository, folder, "/", sep="", collapse="")
      cmd <- sprintf("echo \"ls\" | ftp %s", where)
      whatsin <- readLines(pipe(cmd))
      closeAllConnections()
    }
    else if( sysinfo == "SunOS" )
    {
      # Build command file for ftp connection:
      tmpname <- tempfile(pattern="getncbi")
      tmpcmdfile <- file(tmpname, open="w")
      writeLines("user anonymous seqteam@biomserv.univ-lyon1.fr", tmpcmdfile)
      writeLines(sprintf("cd genomes/Bacteria/%s", folder), tmpcmdfile)
      writeLines("dir", tmpcmdfile)
      writeLines("bye", tmpcmdfile)
      close(tmpcmdfile)
      #
      hostname <- unlist(strsplit(repository,split="/"))[3]
      cmd <- sprintf("ftp -v -n %s < %s", hostname, tmpname)
      whatsin <- readLines(pipe(cmd))
    }
    else if( sysinfo=="Linux"){
      # Build command file for ftp connection:
      tmpname <- tempfile(pattern="getncbi")
      tmpcmdfile <- file(tmpname, open="w")
      writeLines("user anonymous seqteam@biomserv.univ-lyon1.fr", tmpcmdfile)
      writeLines(sprintf("cd genomes/Bacteria/%s", folder), tmpcmdfile)
      writeLines("dir", tmpcmdfile)
      writeLines("bye", tmpcmdfile)
      close(tmpcmdfile)
      #
      hostname <- unlist(strsplit(repository,split="/"))[3]
      cmd <- sprintf("ftp -v -n %s < %s", hostname, tmpname)
      whatsin <- readLines(pipe(cmd))
    }
    else
    {
      stop("unimplemented platform")
    }
    whatsin <- whatsin[ grep("\\.gbk", whatsin)] # Keep only files with ".gbk" extension
    #
    # Remove backup files with % extension:
    #
    for( i in seq_len(length(whatsin)) )
      if( substr(whatsin[i], nchar(whatsin[i]), nchar(whatsin)[i]) == "%" )
        is.na(whatsin[i]) <- TRUE
    whatsin <- whatsin[!is.na(whatsin)]
    
    for( i in seq(from=1, to=length(whatsin), by=1 )) # Loop on sequences data
    {
      #
      # Try to get the accession number of this entry:
      #
      accname <- unlist(strsplit(whatsin[i], split=" "))
      accname <- accname[length(accname)]
      accname <- unlist(strsplit(accname, split="\\."))
      accname <- accname[1] # The accession number should be in this variable
      #
      # Try to get the size of this entry:
      #
      entry <- paste(repository, folder, "/", accname, ".gbk", sep="", collapse="")
      header <- readLines(entry, n=2)
      closeAllConnections()
      bp <- unlist(strsplit(header[1], split=" "))
      bp <- bp[nchar(bp) > 0]
      bp <- bp[3] # size in bp should be there
      #
      # Try to get the last update date of this entry
      #
      last <- unlist(strsplit(header[1], split=" "))
      last <- last[nchar(last) > 0]
      last <- last[length(last)] # last update time should be there

      #
      # Try to get the type (chromosome versus plasmid) of this entry:
      #
      if( length(grep("plasmid", tolower(header[2]))) != 0 )
        def <- "plasmid"
      else if(length(grep("chromosome", tolower(header[2]))) != 0)
        def <- "chromosome"
      else if(length(grep("genome", tolower(header[2]))) != 0)
        def <- "chromosome"
      else
        def <- NA
      #
      # Begin the human curated part:
      #
      if( accname == "NC_002528" ) def <- "chromosome"
      if( accname == "NC_003454" ) def <- "chromosome"
      if( accname == "NC_001732" ) def <- "plasmid"
      if( accname == "NC_001733" ) def <- "plasmid"
      if( accname == "NC_005042" ) def <- "chromosome"
      if( accname == "NC_005072" ) def <- "chromosome"
      if( accname == "NC_004631" ) def <- "chromosome"
      if( accname == "NC_004344" ) def <- "chromosome"
      if( accname == "NC_003902" ) def <- "chromosome"
      if( accname == "NC_005957" ) def <- "chromosome"
      if( accname == "NC_007984" ) def <- "chromosome"
      if( accname == "NC_002937" ) def <- "chromosome"
      if( accname == "NC_005863" ) def <- "plasmid"
      if( accname == "NC_008054" ) def <- "chromosome"
      if( accname == "NC_002942" ) def <- "chromosome"
      if( accname == "NC_005823" ) def <- "chromosome"
      if( accname == "NC_005824" ) def <- "chromosome"
      if( accname == "NC_000916" ) def <- "chromosome"
      if( accname == "NC_007633" ) def <- "chromosome"
      if( accname == "NC_006855" ) def <- "plasmid"
      if( accname == "NC_006856" ) def <- "plasmid"
      if( accname == "NC_006905" ) def <- "chromosome"
      if( accname == "NC_006511" ) def <- "chromosome"
      if( accname == "NC_003198" ) def <- "chromosome"
      if( accname == "NC_007350" ) def <- "chromosome"
      if( accname == "NC_007351" ) def <- "plasmid"
      if( accname == "NC_007352" ) def <- "plasmid"
      if( accname == "NC_003425" ) def <- "plasmid"
      if( accname == "NC_006833" ) def <- "chromosome"
      if( accname == "NC_006810" ) def <- "chromosome"
      if( accname == "NC_008529" ) def <- "chromosome"
      if( accname == "NC_008531" ) def <- "chromosome"
      if( accname == "NC_008346" ) def <- "chromosome"
      if( accname == "NC_008800" ) def <- "chromosome"
            
            

      
      #
      # Concatenate results:
      #
      species <- c(species, folder)
      accession <- c( accession, accname)
      size.bp <- c(size.bp, as.integer(bp))
      lastupdate <- c(lastupdate, last)
      type <- c(type, def)
      cat("\n",folder,accname,bp,def,last,"\n")
    }
  }

# shouldn't ftp_proxy be restored there ?
  return(data.frame(I(species), I(accession), size.bp, type, I(lastupdate)))
}




########################################################################
#             ncbi.fna.url 
#
#  Try to build urls to access complete genome sequences data
#  in fasta format from get.ncbi() output
#
########################################################################

ncbi.fna.url <- function( get.ncbi.out = get.ncbi() )
{
  build.url <- function( x )
  {
    ficname <- unlist(strsplit(x[2],"\\.")) # split prefix and suffix
    ficname <- ficname[1] # keep prefix
    ficname <- paste( ficname, ".fna", collapse="", sep="")
    urlname <- paste("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/", x[1],
                 "/",ficname, collapse="", sep="")
    return(urlname)
  }
  apply( get.ncbi.out, 1, build.url )  
}

########################################################################
#
#  Try to build urls to access complete genome sequences data
#  in genbank format from get.ncbi() output
#
########################################################################

ncbi.gbk.url <- function( get.ncbi.out = get.ncbi() )
{
  build.url <- function( x )
  {
    urlname <- paste("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/", x[1],
                 "/",x[2], collapse="", sep="")
    return(urlname)
  }
  apply( get.ncbi.out, 1, build.url )  
}
########################################################################
#  Try to build urls to access complete genome sequences data
#  file *.ptt from get.ncbi() output
#
########################################################################

ncbi.ptt.url <- function( get.ncbi.out = get.ncbi() )
{
  build.url <- function( x )
  {
    ficname <- unlist(strsplit(x[2],"\\.")) # split prefix and suffix
    ficname <- ficname[1] # keep prefix
    ficname <- paste( ficname, ".ptt", collapse="", sep="")
    urlname <- paste("ftp://ftp.ncbi.nih.gov/genomes/Bacteria/", x[1],
                 "/",ficname, collapse="", sep="")
    return(urlname)
  }
  apply( get.ncbi.out, 1, build.url )  
}

########################################################################
#
# Try to get the number of cds and genome size
#
########################################################################

ncbi.stats <- function( get.ncbi.out = get.ncbi() )
{
  gbkurls <- ncbi.gbk.url( get.ncbi.out )
  ptturls <- ncbi.ptt.url( get.ncbi.out )
  get.genome.size <- function( url )
  {
    header <- readLines( url, n = 1 )
    tmp <- unlist(strsplit(header, split=" "))
    tmp <- tmp[nchar(tmp)>0]
    tmp <- tmp[3]
    as.integer(tmp)
  }
  get.n.prot <- function( url )
  {
    lines <- readLines( url, n = 3 )
    lines <- lines[3]
    lines <- unlist(strsplit(lines, split=" "))
    print(lines[1])
    as.integer(lines[1])
  }

  gsize <- sapply( gbkurls, get.genome.size)
  nprot <- sapply( ptturls, get.n.prot )
  data <- data.frame( cbind(get.ncbi.out[,1], gsize, nprot) )
  return( data )
}
