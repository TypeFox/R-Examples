Muscle <- function(file.name, file.path="", job.name=NULL, args=NULL,
                   aln.filetype="PHYLIP_INT", shared.username=NULL, out.name=NULL,
                   suppress.Warnings=FALSE, print.curl=FALSE) {

  fm=NULL
  aln.filetype <- match.arg(aln.filetype, c("PHYLIP_INT", "PHYLIP_SEQ", "HTML", "FASTA", "CLUSTALW", "MSF"))

  if (rplant.env$api == "f") {
    privAPP=FALSE
    version="Muscle-3.8.32u4"
    
    if (aln.filetype == "PHYLIP_INT"){
      args <- append(args, "-phyiout")
      aln.name <- "phylip_interleaved.aln"
    } else if (aln.filetype == "PHYLIP_SEQ"){
      args <- append(args, "-physout")
      aln.name <- "phylip_sequential.aln"
    } else if (aln.filetype == "FASTA"){
      args <- append(args, "-fastaout")
      aln.name <- "fasta.aln"
    } else if (aln.filetype == "CLUSTALW"){
      args <- append(args, "-clwout")
      aln.name <- "clustalw.aln"
    } else if (aln.filetype == "MSF"){
      args <- append(args, "-msfout")
      aln.name <- "msf.aln"
    } else if (aln.filetype == "HTML"){
      args <- append(args, "-htmlout")
      aln.name <- "html.aln"
    }
  
    args <- paste(args, collapse=" ")  # make a single statement
  
    App <- GetAppInfo(version)[[3]]
    input.list <- vector("list",1)
    input.list[[1]] <- App[,2][1]
  
    if (is.null(job.name))
      job.name <- paste(rplant.env$user,"_",version,"_viaR", sep="")
  
    args <- list(c("arguments",args))

  } else {
    version="Muscle-hpc-3.8.31u1"  

    if (aln.filetype == "PHYLIP_INT"){
      fm <- append(fm, "-phyiout")
    } else if (aln.filetype == "PHYLIP_SEQ"){
      fm <- append(fm, "-physout")
    } else if (aln.filetype == "FASTA"){
      fm <- append(fm, "-fastaout")
    } else if (aln.filetype == "CLUSTALW"){
      fm <- append(fm, "-clwout")
    } else if (aln.filetype == "MSF"){
      fm <- append(fm, "-msfout")
    } else if (aln.filetype == "HTML"){
      fm <- append(fm, "-htmlout")
    }

    if (is.null(out.name)) {
      if (aln.filetype == "PHYLIP_INT"){
        out.name <- "phylip_interleaved.aln"
      } else if (aln.filetype == "PHYLIP_SEQ"){
        out.name <- "phylip_sequential.aln"
      } else if (aln.filetype == "FASTA"){
        out.name <- "fasta.aln"
      } else if (aln.filetype == "CLUSTALW"){
        out.name <- "clustalw2.aln"
      } else if (aln.filetype == "MSF"){
        out.name <- "msf.aln"
      } else if (aln.filetype == "HTML"){
        out.name <- "html.aln"
      }
    } 

    App <- GetAppInfo(version)[[3]]
    input.list <- vector("list",1)
    input.list[[1]] <- App[,2][1]
   #Sinput.list[[1]] <- "stdin"

    if (is.null(job.name)) {
#     job.name <- paste(rplant.env$user,"_",version,"_viaR", sep="")
      job.name <- version
    }
    
    if (!is.null(args)){
      args <- paste(args, collapse=" ")  # make a single statement
      args <- list(c("format", fm), c("outname", out.name), c("arguments",args))
    } else {
      args <- list(c("format", fm), c("outname", out.name))  
    }

  }

  myJob<-SubmitJob(application=version, job.name=job.name, args.list=args,
                   file.list=list(file.name), file.path=file.path,
                   input.list=input.list, suppress.Warnings=suppress.Warnings,
                   print.curl=print.curl, shared.username=shared.username)

  cat(paste("Result file: ", out.name, "\n", sep=""))
  return(myJob)
}
