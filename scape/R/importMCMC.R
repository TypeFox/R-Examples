importMCMC <- function(dir, coda=FALSE, quiet=TRUE, pretty.labels=FALSE, l.choose=NULL, p.choose=NULL)
{
  ## 1  Define functions
  prettyL <- function(l)
  {
    l[l=="CPUE"           ] <- "CPUE"
    l[l=="CALL"           ] <- "CAc"
    l[l=="CLLL"           ] <- "CLc"
    l[l=="Survey"         ] <- "Survey"
    l[l=="surveyCALL"     ] <- "CAs"
    l[l=="surveyCLLL"     ] <- "CLs"
    l[l=="surveynosexCLLL"] <- "CLn"
    l[l=="negpen"         ] <- "Negative"
    l[l=="Priorlike"      ] <- "Prior"
    l[l=="f"              ] <- "Total"
    return(l)
  }

  prettyP <- function(p)
  {
    p <- gsub("plusscale",           "Plus",     p)
    p <- gsub("Sfullest",            "cSfull",   p)
    p <- gsub("Sfulldelta",          "cSdelta",  p)
    p <- gsub("log_varLest",         "cSleft",   p)
    p <- gsub("log_varRest",         "cSright",  p)
    p <- gsub("log_qCPUE",           "logqCPUE", p)
    p <- gsub("log_qsurvey",         "logq",     p)
    p <- gsub("surveySfull",         "sSfull",   p)
    p <- gsub("surveySfulldeltaest", "sSdelta",  p)
    p <- gsub("log_surveyvarL",      "sSleft",   p)
    p <- gsub("log_surveyvarR",      "sSright",  p)
    p <- if(!any(grep("_2",p))) gsub("_1","",p) else p  # remove _1 if no _2 exist
    return(p)
  }

  getL <- function()
  {
    if(!quiet) cat("Likelihoods  ")
    L <- read.table(paste(dir,"mcmclike.out",sep="/"), header=TRUE)
    if(!quiet) cat("file...")
    if(pretty.labels)
    {
      names(L) <- prettyL(names(L))
      if(!quiet) cat("labels...")
    }
    if(!is.null(l.choose))
    {
      L <- L[,match(l.choose,names(L))]
      if(!quiet) cat("choose...")
    }
    if(!quiet) cat("OK\n")
    return(L)
  }

  getP <- function()
  {
    if(!quiet) cat("Parameters   ")
    P <- read.table(paste(dir,"params.pst",sep="/"), header=TRUE)
    if(!quiet) cat("file...")
    if(pretty.labels)
    {
      names(P) <- prettyP(names(P))
      if(!quiet) cat("labels...")
    }
    if(!is.null(p.choose))
    {
      if(!identical(p.choose, ""))
        P <- P[,match(p.choose, names(P))]
      if(!quiet) cat("choose...")
    }
    if(!quiet) cat("OK\n")
    return(P)
  }

  getB <- function()
  {
    if(!quiet) cat("Biomass      ")
    B <- read.table(paste(dir,"spawbiom.pst",sep="/"), header=TRUE)
    if(!quiet) cat("file...")
    names(B) <- substring(names(B), 8)
    if(!quiet) cat("labels...")
    if(!quiet) cat("OK\n")
    return(B)
  }

  getR <- function()
  {
    if(!quiet) cat("Recruitment  ")
    R <- read.table(paste(dir,"recruits.pst",sep="/"), header=TRUE, fill=TRUE)
    if(!quiet) cat("file...")
    R <- R[,-ncol(R)]
    names(R) <- substring(names(R), 8)
    names(R) <- as.integer(names(R)) - 1
    if(!quiet) cat("labels...")
    if(!quiet) cat("OK\n")
    return(R)
  }

  ## 2  Verify that files exist
  files <- paste(dir, c("mcmclike.out","params.pst","spawbiom.pst","recruits.pst"), sep="/")
  sapply(files, function(f)
         if(!file.exists(f)) stop("File ",f," does not exist. Please check the 'dir' argument."))

  ## 3  Parse files
  if(!quiet) cat("\nParsing files in directory ", dir, ":\n\n", sep="")
  L <- getL()
  P <- getP()
  B <- getB()
  R <- getR()
  if(!quiet) cat("\n")

  ## 4  Collect and convert L, P, B, R
  output <- list(L=L, P=P, B=B, R=R)
  if(coda)
    output <- lapply(output, mcmc)

  ## 5  Create attributes
  attr(output, "call") <- match.call()

  return(output)
}
