### An wrapper of RUE13 functions.

gen.reu13.df <- function(seq.string, phi.df = NULL,
    aa.names = .CF.GV$amino.acid,
    split.S = TRUE, drop.X = TRUE, drop.MW = TRUE, drop.1st.codon = TRUE){
  if(is.null(phi.df)){
    phi.df <- cbind(names(seq.string), rep(1, length(seq.string)))
    names(phi.df) <- c("ORF", "phi.value")
  }

  if(split.S){
    if("S" %in% aa.names){ 
      if(! "Z" %in% aa.names){
        aa.names <- c(aa.names, "Z")
      }
    } else{
      split.S <- FALSE
    }
  } else{
    if(all(c("S", "Z") %in% aa.names)){
        split.S <- TRUE
    }
  }

  if(drop.X){
    aa.names <- aa.names[aa.names != "X"]
  }

  if(drop.MW){
    aa.names <- aa.names[!(aa.names %in% c("M", "W"))]
  }

  if(drop.1st.codon){
    seq.string <- lapply(seq.string, function(x){ x[-1] })
  }

  aa.names <- sort(aa.names)
  ret <- build.reu13.df(seq.string, phi.df, aa.names, split.S = split.S)

  ret <- rearrange.reu13.df(ret)
  ret
} # End of gen.reu13.df().


build.reu13.df <- function(seq.string, phi.df, aa.names, split.S = split.S){
  ### Subset genes which are both in common.
  names.seq <- as.character(names(seq.string))
  names.ORF <- as.character(phi.df[, 1])
  phi.df <- phi.df[names.ORF %in% names.seq,]
  n.seq <- length(seq.string)

  ### Reorder genes and phi.df.
  seq.string <- seq.string[order(names.seq)]
  # names.seq <- data.frame(ORF = names.seq, stringsAsFactors = FALSE)
  # tmp.phi <- merge(names.seq, phi.df, all.x = TRUE)

  ### Find AA position for all genes.
  ret <- lapply(aa.names, build.aa.df,
                seq.string = seq.string, phi.df = phi.df, split.S = split.S)
  # ret <- list()
  # for(i.aa in 1:length(aa.names)){
  #   ret[[i.aa]] <- build.aa.df(aa.names[i.aa], seq.string, phi.df, split.S)
  # }
  names(ret) <- aa.names

  ret
} # End of build.reu13.df()


build.aa.df <- function(aa, seq.string, phi.df, split.S = TRUE){
  if(split.S){
    synonymous.codon <- .CF.GV$synonymous.codon.split[[aa]]
  } else{
    synonymous.codon <- .CF.GV$synonymous.codon[[aa]]
  }

  ### Build all AA data.frame from aa.names.
  # func.get.seq.df <- function(i.gene){
  #   id <- which(seq.string[[i.gene]] %in% synonymous.codon)
  #   n.match <- length(id)
  #   if(n.match > 0){
  #     tmp <- data.frame(
  #              ORF = rep(as.character(phi.df[i.gene, 1]), n.match),
  #              phi = rep(as.double(phi.df[i.gene, 2]), n.match),
  #              # Amino.Acid = rep(as.character(aa), n.match),
  #              Pos = as.double(id),
  #              Codon = as.character(seq.string[[i.gene]][id]),
  #              stringsAsFactors = FALSE)
  #   } else{
  #     tmp <- data.frame(ORF = NULL, phi = NULL,
  #                       # Amino.Acid = NULL,
  #                       Pos = NULL, Codon = NULL,
  #                       stringsAsFactors = FALSE)
  #   }
  #   tmp
  # } # End of func.get.seq.df().
  # if(.CF.CT$parallel[1] == "lapply"){
  #   ret <- lapply(1:length(seq.string), func.get.seq.df)
  # } else if(.CF.CT$parallel[1] == "mclapply"){
  #   ret <- parallel::mclapply(1:length(seq.string), func.get.seq.df)
  # } else if(.CF.CT$parallel[1] == "task.pull"){
  #   ret <- pbdMPI::task.pull(1:length(seq.string), func.get.seq.df)
  # } else if(.CF.CT$parallel[1] == "pbdLapply"){
  #   ret <- pbdMPI::pbdLapply(1:length(seq.string), func.get.seq.df)
  # } else{
  #   stop("parallel method is not found.")
  # }
  ret <- lapply(1:length(seq.string),
           function(i.gene){
             id <- which(seq.string[[i.gene]] %in% synonymous.codon)
             n.match <- length(id)
             if(n.match > 0){
               tmp <- data.frame(
                        ORF = rep(as.character(phi.df[i.gene, 1]), n.match),
                        phi = rep(as.double(phi.df[i.gene, 2]), n.match),
                        # Amino.Acid = rep(as.character(aa), n.match),
                        Pos = as.double(id),
                        Codon = as.character(seq.string[[i.gene]][id]),
                        stringsAsFactors = FALSE)
             } else{
               tmp <- data.frame(ORF = NULL, phi = NULL,
                                 # Amino.Acid = NULL,
                                 Pos = NULL, Codon = NULL,
                                 stringsAsFactors = FALSE)
             }
             tmp
           })

  ret <- do.call("rbind", ret)

  ret
} # End of find.aa().

