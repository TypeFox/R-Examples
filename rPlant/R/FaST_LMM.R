FaST_LMM <- function(input.file.list="", ALL.file.path="", print.curl=FALSE,
                     sim.file.list=NULL, pheno.file.name=NULL, mpheno=1,
                     args=NULL, covar.file.name=NULL, job.name=NULL, 
                     shared.username=NULL, suppress.Warnings=FALSE,
                     out.basename=NULL) {

  if (rplant.env$api == "a") {
#     privAPP=TRUE
#     version="FaST-LMM-beta-2.07"

    privAPP=FALSE
    version="FaST-LMM-hpc-2.07u1"
  } else {
    privAPP=FALSE
    version="FaST-LMM-1.09u1"
  }

  input.len <- length(input.file.list)
  input.list <- list()
  if ((input.len) == 3){
    ext1 <- unlist(strsplit(input.file.list[[1]], "\\."))[2]
    ext2 <- unlist(strsplit(input.file.list[[2]], "\\."))[2]
    ext3 <- unlist(strsplit(input.file.list[[3]], "\\."))[2]
    input.type="B"
    input.list[[1]] <- find.input(ext1)
    input.list[[2]] <- find.input(ext2)
    input.list[[3]] <- find.input(ext3)
  } else {
    ext1 <- unlist(strsplit(input.file.list[[1]], "\\."))[2]
    ext2 <- unlist(strsplit(input.file.list[[2]], "\\."))[2]
    input.list[[1]] <- find.input(ext1)
    input.list[[2]] <- find.input(ext2)
    if ((ext1 == "tfam" ) || (ext1 == "tped")) {
      input.type="T"
    } else {
      input.type="R"
    }
  }

  args <- c(args)

  if (input.type=="T"){
    options <- list(c("T",TRUE))
    if (!is.null(sim.file.list)){
      input.file.list <- append(input.file.list,sim.file.list)
      input.list <- append(input.list,c("Sim1","Sim2"))
      options <- append(options,list(c("S",TRUE)))
    }
  } else if (input.type=="B") {
    options <- list(c("B",TRUE))
    if (!is.null(sim.file.list)){
      input.file.list <- append(input.file.list,sim.file.list)
      input.list <- append(input.list,c("Sim1","Sim2","Sim3"))
      options <- append(options,list(c("S",TRUE)))
    }
  } else {
    options <- NULL
    if (!is.null(sim.file.list)){
      input.file.list <- append(input.file.list,sim.file.list)
      input.list <- append(input.list,c("Sim1","Sim2"))
      options <- append(options,list(c("S",TRUE)))
    }
  }

  if (is.null(out.basename)){
    out.basename <- unlist(strsplit(input.file.list[[1]], "\\."))[1]
  }

  args <- append(args, c("-out",out.basename))

  if (is.null(job.name)){
    job.name <- out.basename
  }

  if (!is.null(pheno.file.name)){
    input.file.list <- append(input.file.list, pheno.file.name)
    input.list <- append(input.list,"inputPHENO")
    options <- append(options,list(c("P",TRUE),c("mpheno",mpheno)))
  }

  if (!is.null(covar.file.name)){
    input.file.list <- append(input.file.list, covar.file.name)
    input.list <- append(input.list,"inputCOVAR")
    options <- append(options,list(c("C",TRUE)))
  }

  # make a single statement
  args <- paste(args, collapse=" ") 
  options <- append(options, list(c("arguments",args)))
  # Submit
  myJob<-SubmitJob(application=version, args.list=options, job.name=job.name,
                   file.list=input.file.list, file.path=ALL.file.path, private.APP=privAPP,
                   input.list=input.list, shared.username=shared.username,
                   print.curl=print.curl, suppress.Warnings=suppress.Warnings)

  return(myJob)

}
