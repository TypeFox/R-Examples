check.caRpools = function(packages=TRUE, files=TRUE, mageck=TRUE, bowtie2=TRUE, pandoc=TRUE, skip.updates=TRUE, template=NULL, scripts=TRUE, miaccs="MIACCS.xls")
{
  # check.caRpools will check for correct installation of caRpools to be used for data analysis and report generation
  # it is divided into checking packages, files, mageck and bowtie2 installation
  return.val=TRUE
  # check packages
  if(identical(packages,TRUE))
  {
    check.packages = load.packages(noupdate=skip.updates)
    
    if(length(check.packages) == 0)
    {
      stop("Package verification failed. Please make sure you are connected to the internet and all packages have been installed correctly.")
      return.val=FALSE
    }
    
  }
  # check MAGeCK
  if(identical(mageck,TRUE))
  {
    #open mageck command
    check = system("mageck -h",intern=FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE,
                   wait = TRUE, input = NULL, show.output.on.console = FALSE,
                   minimized = FALSE, invisible = TRUE)
    
    if(check != 0)
    {
      stop("MAGeCK does not work. Please make sure MAGeCK is installed correctly as explained on the MAGeCK homepage.")
      return.val=FALSE
    }
    
  }
  
  # check Pandoc
  if(identical(pandoc,TRUE))
  {
    #open mageck command
    check = system("pandoc -h",intern=FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE,
                   wait = TRUE, input = NULL, show.output.on.console = FALSE,
                   minimized = FALSE, invisible = TRUE)
    
    if(check != 0)
    {
      stop("Pandoc does not work. Please make sure Pandoc is installed correctly as explained on the Pandoc.org homepage.")
      return.val=FALSE
    }
    
  }
  
  # check Bowtie2 and bowtie2-build
  if(identical(bowtie2,TRUE))
  {
    #open Bowtie2 command
    check = system("bowtie2 -h",intern=FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE,
                   wait = TRUE, input = NULL, show.output.on.console = FALSE,
                   minimized = FALSE, invisible = TRUE)
    
    if(check != 0)
    {
      stop("Bowtie2 does not work. Please make sure Bowtie2 is installed correctly as explained in the Bowtie2 manual.")
      return.val=FALSE
    }
    
    #open Bowtie2-build command
    check = system("bowtie2-build -h",intern=FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE,
                   wait = TRUE, input = NULL, show.output.on.console = FALSE,
                   minimized = FALSE, invisible = TRUE)
    
    if(check != 0)
    {
      stop("Bowtie2-build does not work. Please make sure Bowtie2 is installed correctly as explained in the Bowtie2 manual.")
      return.val=FALSE
    }
    
  }
  
  # check for files being present
  if(identical(files,TRUE))
  {
    # Check for working directory and file in there:
    # MIACCS.xls
    # RMD file (standards + user defined ones)
    # perl files if necessary
    # return 0 is success or -1 for failure
    miaccs.check = file.access(paste(getwd(), "/", miaccs, sep=""))
    
    # check Readcount files or FASTQ
    miaccs.loaded = load.file(paste(getwd(), "/", miaccs, sep=""), type="xlsx")
    
    if(miaccs.check == 0)
    {
      miaccs.file = load.file(miaccs, type="xlsx")
      scriptpath=as.character(miaccs.file["carpools.scriptpath",3])
      datapath = as.character(miaccs.file["carpools.datapath",3])
      
    }
    else
    { scriptpath = paste(getwd(),"/","scripts", sep="")
      stop(paste("MIACCS.xls not present in the R working directory:", getwd()))
      return.val=FALSE
    }
    
    
    # check data files
    if(identical(as.logical(miaccs.loaded["carpools.extract.fastq",3]), TRUE))
    {
      file1 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.untreated1",3]), ".fastq" , sep=""))
      file2 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.untreated2",3]), ".fastq" , sep=""))
      file3 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.treated1",3]) , ".fastq", sep=""))
      file4 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.treated2",3]) , ".fastq", sep=""))
    }
    else
    {
      file1 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.untreated1",3]) , sep=""))
      file2 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.untreated2",3]) , sep=""))
      file3 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.treated1",3]) , sep=""))
      file4 = file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.treated2",3]) , sep=""))
    }
    
    # name of reference library .fasta file WITHOUT extension .fasta
    referencefile =  file.access(paste(datapath, "/", as.character(miaccs.loaded["carpools.lib.ref",3]), ".fasta" , sep=""))
    
    if(file1 == -1 || file2 == -1 || file3 == -1 || file4 == -1 || referencefile ==-1)
    {
      stop(paste("Any of the sample files or the reference library file could not be found or read within", datapath, ". Maybe you misstyped the name or added the .fasta to the reference file? Please also make sure there is no .fastq in the fastq file name. However, if you provide read count files, you must provide a file ending."))
      return.val=FALSE
    }

    
    
    
    if(identical(scripts,TRUE))
    {
      # check for presence of script files
      # open MIACCS and look for script folder
      
      
      mapping.file = file.access(paste(scriptpath, "/", "CRISPR-mapping.pl", sep=""))
      extraction.file = file.access(paste(scriptpath, "/", "CRISPR-extract.pl", sep=""))
      
      if(mapping.file == -1 || extraction.file == -1)
      {
        stop(paste("Either CRISPR-mapping.pl or CRISPR-extraction.pl not present in the script directory:", scriptpath, sep=" \n") )
        return.val=FALSE
      }
      
    }
    
    # Check for template files (standard one OR provided)
    if(is.null(template))
    {
      # check for standard files, at least one must be present in the R working dir
      if(file.access(paste(getwd(), "/", "CaRpools-extended-PDF.Rmd", sep="")) == 0 || file.access(paste(getwd(), "/", "CaRpools-extended-HTML", sep="")) ==0 || file.access(paste(getwd(), "/", "CaRpools-PDF.Rmd", sep="")) ==0 || file.access(paste(getwd(), "/", "CaRpools-HTML.Rmd", sep="")) ==0)
      {
        # any of them present 
      }
      else
      {
        # any of them not present
        retrun.val=FALSE
        stop(paste("Neither CaRpools-extended-PDF.Rmd, CaRpools-extended-HTML.Rmd, CaRpools-PDF.Rmd or CaRpools-HTML.Rmd have been found or have not access to it within the R working directory:"), getwd(), sep=" \n" )
        
      }
    }
    else
    {
      # check custom Rmd Template file
      if( file.access(paste(getwd(), "/", template, sep="")) != 0 )
      {
        return.val=FALSE
        stop(paste("Your template file",template, "does not provide read/write access or was not found within the R working directory:"), getwd(), sep=" \n" ) 
      }
    }
    
  }
  
  # Return TRUE or FALSE
  return(return.val)
  
}