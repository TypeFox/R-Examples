Fasttree <- function(file.name, file.path="", job.name=NULL, out.name=NULL, 
                     args=NULL, type="DNA", model=NULL, gamma=FALSE, stat=FALSE,
                     print.curl=FALSE, shared.username=NULL, suppress.Warnings=FALSE) {
                     

  # FastTree infers approximately-maximum-likelihood phylogenetic trees 
  #   from alignments of nucleotide or protein sequences.
  #
  # Args:
  #   file.name (type = string): File name of input fasta file
  #   file.path (type = string): Path to where input file is located
  #   job.name (type = string): Job name adds a time stamp to make them unique
  #   version (type = string): ClustalW version
  #   type (type = string): Either 'DNA' or 'PROTEIN'
  #   args (type = string): Optional for arguments (i.e. flags, like on command 
  #     line).  One string with commands separated by spaces. see help(ClustalW)
  #     for details.
  #   model (type = string): Name of substitution model.  For DNA the choices 
  #     are GTRCAT, and JCCAT.  For protein the choices are JTTCAT, and WAGCAT.
  #   gamma (type = boolean): Use this option (about 5\% slower) if you want 
  #     to rescale the branch lengths and compute a Gamma20-based likelihood.
  #   stat (type = boolean): This allows for statistical comparisons (when 
  #     gamma=TRUE) of the likelihood of different topologies.
  #   print.curl (type = string): Prints the associated curl statement
  #   version (type = string): ClustalW version
  #   shared.username (type = string): Valid iPlant username with whom the files 
  #     are being shared.
  #   suppress.Warnings (type = boolean): Don't do any error checking
  #   email (type = boolean): By default an email will be sent to the user when 
  #     the job finishes.

  # Returns:
  #   Returns the job id (number).  o/w an error

  type <- match.arg(type, c("DNA", "PROTEIN"))

  if (type == "DNA"){

    if (is.null(model)){
      model="GTRCAT"
    }

    model <- match.arg(model, c("GTRCAT", "JCCAT"))

#   if (is.null(job.name)){
#     job.name <- paste(rplant.env$user, "_Fasttreedna_", model, "_viaR",
#                       sep="")
#   }

  } else {

    if (is.null(model)){
      model="JTTCAT"
    }

    model <- match.arg(model, c("JTTCAT", "WAGCAT"))

#   if (is.null(job.name)){
#     job.name <- paste(rplant.env$user, "_Fastreeprotein_", model, "_viaR", 
#                       sep="")
#   }
  }

  if (!is.null(args)) {
    args <- c(args)
  }

  if (gamma == TRUE) {
    args <- append(args, "-gamma")
  }

  if (stat == TRUE) {
    args <- append(args, "-log logfile")
  }
    
  if (rplant.env$api == "f") {

    privAPP=FALSE
    version="fasttreeDispatcher-1.0.0u1"

    if (model == "GTRCAT"){
      args <- append(args, "-gtr -nt")
    } else if (model == "WAGCAT"){
      args <- append(args, "-wag")
    }
  
    args <- paste(args, collapse=" ")  # make a single statement
  
    args <- list(c("arguments",args))

  } else {

    version="FasttreeDispatcher-2.1.4u2" 

    options=NULL
  
    if (model == "GTRCAT"){
      options <- append(options, list(c("model", "-gtr -nt")))
    } else if (model == "WAGCAT"){
      options <- append(options, list(c("model", "-wag")))
    }

    if (!is.null(out.name)) {
      options <- append(options, list(c("outname", out.name)))
    }
  
    if (!is.null(args)){
      args <- c(args)
      args <- append(options, list(c("arguments",args)))
    } else {
      args <- options
    }
    

  }
  
  nprocs <- 1
  App <- GetAppInfo(version)[[3]]
  input.list <- vector("list",1)
  input.list[[1]] <- App[,2][1]

  if (is.null(job.name)){
#   job.name <- paste(rplant.env$user,"_",version,"_viaR", sep="")
    job.name <- version
  }

  myJob<-SubmitJob(application=version, job.name=job.name, nprocs=nprocs,
                   file.list=list(file.name), file.path=file.path, 
                   input.list=input.list, suppress.Warnings=suppress.Warnings,
                   print.curl=print.curl, shared.username=shared.username,
                   args.list=args)

  return(myJob)
}
