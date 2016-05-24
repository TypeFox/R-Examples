ClustalW <- function(file.name, file.path="", type="DNA", aln.filetype="CLUSTALW", 
                     args=NULL, out.name=NULL, job.name=NULL, print.curl=FALSE,   
                     shared.username=NULL, suppress.Warnings=FALSE) {

  #  An approach for performing multiple alignments of large numbers of
  #   amino acid or nucleotide sequences is described. Input is a fasta file
  #   output is an alignment.
  #
  # Args:
  #   file.name (type = string): File name of input fasta file
  #   file.path (type = string): Path to where input file is located
  #   type (type = string): Either 'DNA' or 'PROTEIN'
  #   args (type = string): Optional for arguments (i.e. flags, like on command 
  #     line).  One string with commands separated by spaces. see help(ClustalW)
  #     for details.
  #   job.name (type = string): Job name adds a time stamp to make them unique
  #   a', ln.filetype (type = string): Seven options: 'CLUSTALW', 'PHYLIP_INT',
  #     'FASTA', 'NEXUS', 'GCG', 'GDE', 'PIR'.  This is output filetype
  #   print.curl (type = string): Prints the associated curl statement
  #   version (type = string): ClustalW version
  #   shared.username (type = string): Valid iPlant username with whom the files 
  #     are being shared.
  #   suppress.Warnings (type = boolean): Don't do any error checking
  #   email (type = boolean): By default an email will be sent to the user when 
  #     the job finishes.
  #
  # Returns:
  #   Returns the job id (number).  o/w an error

  type <- match.arg(type, c("DNA", "PROTEIN"))

  aln.filetype <- match.arg(aln.filetype, c("CLUSTALW", "PHYLIP_INT", "NEXUS", 
                                            "FASTA", "GCG", "GDE", "PIR"))
                                            
  if (rplant.env$api == "f") {
    privAPP=FALSE
    version="ClustalW2-2.1u1"
    
    # Get the appropriate flags for the function
    if (type == "DNA"){
      args <- append(args, "-TYPE=DNA")
    } else {
      args <- append(args, "-TYPE=PROTEIN")
    }
  
    if (aln.filetype == "FASTA"){
      args <- append(args, "-OUTPUT=FASTA")
    } else if (aln.filetype == "PHYLIP_INT"){
      args <- append(args, "-OUTPUT=PHYLIP")
    } else if (aln.filetype == "NEXUS"){
      args <- append(args, "-OUTPUT=NEXUS")
    } else if (aln.filetype == "GCG"){
      args <- append(args, "-OUTPUT=GCG")
    } else if (aln.filetype == "GDE"){
      args <- append(args, "-OUTPUT=GDE")
    } else if (aln.filetype == "PIR"){
      args <- append(args, "-OUTPUT=PIR")
    }    
      
    args <- paste(args, collapse=" ")  # make a single statement
  
    args <- list(c("arguments",args))
    out.name <- "clustalw2.aln"
  } else {
    version="Clustalw-2.1.0u2"  
    options=NULL
    # Get the appropriate flags for the function
    if (type == "DNA"){
      options <- append(options, list(c("T", "DNA")))
    } else {
      options <- append(options, list(c("T", "PROTEIN")))
    }
  
    if (aln.filetype == "CLUSTALW"){
      options <- append(options, list(c("format", "CLUSTAL")))
      if (is.null(out.name)) {
        out.name <- "clustalw2.aln"
      }
    } else if (aln.filetype == "FASTA"){
      options <- append(options, list(c("format", "FASTA"))) 
      if (is.null(out.name)) {
        out.name <- "fasta.aln"
      }
    } else if (aln.filetype == "PHYLIP_INT"){
      options <- append(options, list(c("format", "PHYLIP")))
      if (is.null(out.name)) {
        out.name <- "phylip_interleaved.aln"
      }
    } else if (aln.filetype == "NEXUS"){
      options <- append(options, list(c("format", "NEXUS")))
      if (is.null(out.name)) {
        out.name <- "nexus.aln"
      }
    } else if (aln.filetype == "GCG"){
      options <- append(options, list(c("format", "GCG")))
      if (is.null(out.name)) {
        out.name <- "gcg.aln"
      }
    } else if (aln.filetype == "GDE"){
      options <- append(options, list(c("format", "GDE")))
      if (is.null(out.name)) {
        out.name <- "gde.aln"
      }
    } else if (aln.filetype == "PIR"){
      options <- append(options, list(c("format", "PIR")))
      if (is.null(out.name)) {
        out.name <- "pir.aln"
      }
    }  
    options <- append(options, list(c("outname", out.name))) 

    if (!is.null(args)) {
      args <- paste(args, collapse=" ")  # make a single statement
      args <- append(options,list(c("arguments",args)))
    } else {
      args <- options
    }
  }

  # Get proper input.list
  nprocs=1
  App <- GetAppInfo(version)[[3]]
  input.list <- vector("list",1)
  input.list[[1]] <- App[,2][1]
  
  if (is.null(job.name)) {
#   job.name <- paste(rplant.env$user,"_",version,"_viaR", sep="")
    job.name <- version
  }
    
  myJob<-SubmitJob(application=version, job.name=job.name, nprocs=nprocs, 
                   file.list=list(file.name), file.path=file.path, 
                   input.list=input.list, suppress.Warnings=suppress.Warnings,
                   print.curl=print.curl, shared.username=shared.username,
                   args.list=args)

  cat(paste("Result file: ", out.name, "\n", sep=""))
  return(myJob)
}
