PLINK <- function(file.list="", file.path="", job.name=NULL, out.basename=NULL,
                  assoc.method="Standard", no.sex=TRUE, args=NULL,
                  print.curl=FALSE, multi.adjust=TRUE, shared.username=NULL, 
                  suppress.Warnings=FALSE) {
                  
  assoc.method <- match.arg(assoc.method, c("Standard", "QTM", "GxE", "Fisher", "AlternateModel", "Breslow-Day", "Homog", "Linear", "Logistic", "NONE"))
   
  if (rplant.env$api == "a") {
#     privAPP=TRUE
#     version="plink-beta-1.07"

    privAPP=FALSE
    version="PLINK-hpc-1.07u1"
  } else {
    privAPP=FALSE
    version="plink-1.07u1"
  }
  
  input.len <- length(file.list)
  input.list <- list()
  if ((input.len) == 3){
    ext1 <- unlist(strsplit(file.list[[1]], "\\."))[2]
    ext2 <- unlist(strsplit(file.list[[2]], "\\."))[2]
    ext3 <- unlist(strsplit(file.list[[3]], "\\."))[2]
    input.type="B"
    input.list[[1]] <- find.input(ext1)
    input.list[[2]] <- find.input(ext2)
    input.list[[3]] <- find.input(ext3)
  } else {
    ext1 <- unlist(strsplit(file.list[[1]], "\\."))[2]
    ext2 <- unlist(strsplit(file.list[[2]], "\\."))[2]
    input.list[[1]] <- find.input(ext1)
    input.list[[2]] <- find.input(ext2)
    if ((ext1 == "tfam" ) || (ext1 == "tped")) {
      input.type="T"
    } else {
      input.type="R"
    }
  }  

  if (assoc.method == "Standard") {
    method="--assoc"
  } else if (assoc.method == "QTM") {
    method="--assoc --qt-means"
  } else if (assoc.method == "GxE") {
    method="--gxe"
  } else if (assoc.method == "Fisher") {
    method="--fisher"
  } else if (assoc.method == "AlternateModel") {
    method="--model"
  } else if (assoc.method == "Breslow-Day") {
    method="--bd"
  } else if (assoc.method == "Homog") {
    method="--homog"
  } else if (assoc.method == "Linear") {
    method="--linear"
  } else if (assoc.method == "Logistic") {
    method="--logistic"
  } else {
    method=""
  }
  
  if (multi.adjust){args <- append(args, c("--adjust"))}

  if (no.sex){args <- append(args, c("--allow-no-sex"))}

  if (input.type=="T"){
    options <- list(c("T",TRUE))
    if (is.null(out.basename)){
      BASE1 <- substr(file.list[[1]],1,nchar(file.list[[1]])-5)
      EXT1 <- substr(file.list[[1]],nchar(file.list[[1]])-3,nchar(file.list[[1]]))
      BASE2 <- substr(file.list[[2]],1,nchar(file.list[[2]])-5)
      if (EXT1 == "tfam"){
        out.basename <- paste(BASE1,"_",BASE2,"_", assoc.method, sep="")
      } else {
        out.basename <- paste(BASE2,"_",BASE1,"_", assoc.method, sep="")
      }
      out.basename <- paste(unlist(strsplit(out.basename, " ")), collapse="")
    }
  } else if (input.type=="B") {
    options <- list(c("B",TRUE))
    if (is.null(out.basename)){
      BASE1 <- substr(file.list[[1]],1,nchar(file.list[[1]])-4)
      EXT1 <- substr(file.list[[1]],nchar(file.list[[1]])-2,nchar(file.list[[1]]))
      BASE2 <- substr(file.list[[2]],1,nchar(file.list[[2]])-4)
      EXT2 <- substr(file.list[[2]],nchar(file.list[[2]])-2,nchar(file.list[[1]]))
      BASE3 <- substr(file.list[[3]],1,nchar(file.list[[3]])-4)
      EXT3 <- substr(file.list[[3]],nchar(file.list[[3]])-2,nchar(file.list[[1]]))
      if (EXT1 == "bed" && EXT3 == "bam"){
        out.basename <- paste(BASE1, "_", BASE3, "_", BASE2, "_", assoc.method, sep="")
      } else if (EXT1 == "bed" && EXT2 == "bam") {
        out.basename <- paste(BASE1, "_", BASE2, "_", BASE3, "_", assoc.method, sep="")
      } else if (EXT2 == "bed" && EXT3 == "bam"){
        out.basename <- paste(BASE2, "_", BASE3, "_", BASE1, "_", assoc.method, sep="")
      } else if (EXT3 == "bed" && EXT1 == "bam"){
        out.basename <- paste(BASE3, "_", BASE1, "_", BASE2, "_", assoc.method, sep="")
      } else if (EXT2 == "bed" && EXT1 == "bam") {
        out.basename <- paste(BASE2, "_", BASE1, "_", BASE3, "_", assoc.method, sep="")
      } else {
        out.basename <- paste(BASE3, "_", BASE2, "_", BASE1, "_", assoc.method, sep="")
      }
      out.basename <- paste(unlist(strsplit(out.basename, " ")), collapse="")
    }
  } else {
    options <- NULL
    if (is.null(out.basename)){
      BASE1 <- substr(file.list[[1]],1,nchar(file.list[[1]])-4)
      out.basename <- paste(BASE1,"_", assoc.method, sep="")
    }
    out.basename <- paste(unlist(strsplit(out.basename, " ")), collapse="")
  }

  if (is.null(options)) {
    options <- list(c("method", method))
  } else {
    options <- append(options, list(c("method", method)))
  }
  options <- append(options, list(c("inputOUT", out.basename)))
  
  if (is.null(job.name)){
    job.name <- out.basename
  }

  # make a single statement
  args <- paste(args, collapse=" ") 
  options <- append(options, list(c("arguments",args)))

  # Submit
  myJob<-SubmitJob(application=version, args.list=options, job.name=job.name,
                   file.list=file.list, file.path=file.path, 
                   input.list=input.list, shared.username=shared.username,
                   print.curl=print.curl, suppress.Warnings=suppress.Warnings,
                   private.APP=privAPP)

  return(myJob)

}
