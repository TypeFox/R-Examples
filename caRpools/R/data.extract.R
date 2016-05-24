data.extract = function(scriptpath=NULL, datapath=NULL, fastqfile=NULL, extract = FALSE, pattern = "default", machinepattern = "default", createindex = FALSE, referencefile = NULL, mapping = FALSE, reversecomplement = FALSE, threads = 1, bowtieparams = "", sensitivity = "very-sensitive-local", match = "perfect")
  
{
  # ON MAC, Rstudio has to be started different!
  # In the terminal: open -a RStudio
  
  # arguments expected
  #   scriptpath      -> path to where the CRISPR-extract.pl file is
  #   datapath        -> path to where data is expected to be
  #   fastqfile       -> name of fastqfile INCLUDING .fastq
  #   createindex     -> TRUE if bowtie needs to get a mapping index out of a fasta
  #   referencefile   -> filename of FASTA reference to be used if createindex == TRUE
  #                   -> filename of Bowtie2 FASTA index file if createindex == FALSE
  #   #sam            -> filename of SAM file -> is created by this script
  #   
  
  # Output will be filename of readcount files! can be passed on to variable and used for loading
  
  # needs installation of CRISPR-extract.pl for fastq data extraction
  # and                   CRISPR-mapping.pl for mapping
  
  # furthermore, bowtie2 needs to be installed (see bowtie2 manual)
  
  # for data extraction, the following information is necessary
  # $path to perl files
  # $path2 to fastq files from sequencer
  
  # usage of data extraction
  # ARGS:
  # CRISPR-extract.pl
  #   [pattern]
  #		[FastqFile] 
  #   [revcomp]
  
  #CRISPR-mapping.pl
  #		[FASTA of library] 
  #		[Bowtie2 SAM Library]
  
  
  # IDEA
  
  # user can set for extraction
  # -> files are generated and filename is returned at the end
  
  # if user wants data only to be extracted -> only filename is returned
  
  # if user wants to map data
  # -> if index is already there, just map and return filename
  # -> if index is not there, create index and directly passs it to the mapping -> filename is returned at the end
  
  
  if(is.null(scriptpath)) {stop("Path to script not set!")}
  if(is.null(datapath)) {stop("Path to data not set!")}
  if(is.null(fastqfile)) {stop("Name of .fastq file not set!")}
  toreturn = fastqfile
  
  
  if (identical(extract, TRUE))
  {
    
    # CALL
    # PATTERN, FILENAME, REVCOMP
    # First extract .fastq data and convert into .fasta
    extractstring = shQuote(file.path(scriptpath, "CRISPR-extract.pl"))
    filearg = shQuote(file.path(datapath, paste(fastqfile, ".fastq", sep="")))
    
    print(paste("Extraction command is sent to the system:", paste("perl", paste("'",extractstring, "'", sep=""),paste("'", pattern,"'", sep=""), filearg, reversecomplement, paste("'", machinepattern, "'", sep=""), sep=" "), sep=" "))
    system(paste("perl", paste("'",extractstring, "'", sep=""),paste("'", pattern,"'", sep=""), filearg, reversecomplement, paste("'", machinepattern, "'", sep=""), sep=" "))
    
    toreturn = paste(fastqfile, "_extracted.fastq", sep="")
    # change name for further use
    fastqfile = paste(fastqfile, "_extracted", sep="")
  }
  
  
  # Mapping
  if(identical(mapping, TRUE))
  {
    if(is.null(referencefile)) {stop("Reference file name not set!")}
    
    # is index file already there? if not, create index!
    if(identical(createindex, TRUE))
    {
      # get reference file
      filearg = shQuote(file.path(datapath, paste(referencefile, ".fasta", sep="")))
      bt2index = shQuote(file.path(datapath, referencefile))
      # call bowtie 2
      system(paste("bowtie2-build", "-f", filearg, paste(bt2index, bowtieparams, sep=" ")))
      # set name of bowtie2 indexfile
      
    }
    
    # do the mapping using the Bowtie2 Index File
    print(fastqfile)
    print(datapath)
    
    filearg.reference = shQuote(file.path(datapath, referencefile))
    filearg.fastq = shQuote(file.path(datapath, paste(fastqfile, ".fastq", sep="")))
    filearg.sam = shQuote(file.path(datapath, paste(fastqfile, ".sam", sep="")))
    
    bowtie2.params = paste("-p", threads, paste("--", sensitivity, sep=""), bowtieparams, sep=" ")
    
    # Call bowtie2 for mapping
    print(paste("Start Mapping using", paste("bowtie2", "-x", filearg.reference, "-U", filearg.fastq, "-S", filearg.sam,  bowtie2.params, sep=" "), sep=" "))
    system(paste("bowtie2", "-x", filearg.reference, "-U", filearg.fastq, "-S", filearg.sam,  bowtie2.params, sep=" "))
    # extract the information using the reference file and the Bowtie2 created SAM FILE and the perl script CRISPR-bowtie.pl
    
    # The script called is CRISPR-mapping.pl
    # Arguments:
    # fasta file of library
    # Bowtie2 SAM File
    filearg.reference = shQuote(file.path(datapath, paste(referencefile, ".fasta", sep = "")) )
    extractstring = shQuote(file.path(scriptpath, "CRISPR-mapping.pl"))
    print(paste("Start Extracting mapped data:", paste("perl", extractstring, filearg.reference, filearg.sam, match, sep=" "), sep=" "))
    system(paste("perl", extractstring, filearg.reference, filearg.sam, match, sep=" "))
    
    # return filename for designs file that is used later on, including the full path!
    toreturn = file.path(paste(fastqfile, "-designs.txt", sep=""))
    
  }
  return(toreturn)
}