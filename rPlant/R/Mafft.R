Mafft <- function(file.name, file.path="", type="DNA", aln.filetype="FASTA", args=NULL, 
                  out.name=NULL, print.curl=FALSE, job.name=NULL, shared.username=NULL, 
                  suppress.Warnings=FALSE) {

  type <- match.arg(type, c("DNA", "PROTEIN"))

  aln.filetype <- match.arg(aln.filetype, c("CLUSTALW", "FASTA"))

  if (rplant.env$api == "f") {

    privAPP=FALSE
    version="mafftDispatcher-1.0.13100u1"
  
    if (!is.null(args)) {
      args <- c(args)
    }
  
    if (type == "DNA"){
      args <- append(args, "--nuc")
    } else {
      args <- append(args, "--amino")
    }
  
    if (aln.filetype == "CLUSTALW"){
      args <- append(args, "--clustalout")
    }
  
    args <- paste(args, collapse=" ")  # make a single statement
  
    args <- list(c("arguments",args))
    
  } else {
    privAPP=TRUE
    version="mafft-beta-7.221"  
    
    options=NULL
    if (type == "DNA"){
      options <- append(options, list(c("T", "--nuc")))
    } else {
      options <- append(options, list(c("T", "--amino")))
    }
  
    if (aln.filetype == "CLUSTALW"){
      options <- append(options, list(c("format", "--clustalout")))
    }

    if (!is.null(out.name)){
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
#   job.name <- paste(rplant.env$user,"_",version,"_viaR", sep="")\
    job.name <- version
  }

  myJob<-SubmitJob(application=version, job.name=job.name, nprocs=nprocs,
                   file.list=list(file.name), file.path=file.path,
                   input.list=input.list, suppress.Warnings=suppress.Warnings,
                   print.curl=print.curl, shared.username=shared.username,
                   args.list=args, private.APP=privAPP)

  if ((!is.null(out.name)) && (rplant.env$api == "a")){
    cat(paste("Result file:", out.name, "\n"))
  } else if (rplant.env$api == "f"){
    cat("Result file: mafft.fa\n")
  } else {
    cat("Result file: mafft.aln\n")
  }

  return(myJob)
}
