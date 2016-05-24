PerFit.PFS <- function(matrix, method=NULL, simplified=TRUE, 
                   NA.method="Pairwise", Save.MatImp=FALSE, 
                   IP=NULL, IRT.PModel=NULL, Ability=NULL, Ability.PModel=NULL, mu=0, sigma=1)
{
  matrix   <- as.matrix(matrix)
  N        <- dim(matrix)[1]; I <- dim(matrix)[2];
  dico.PFS <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "D.KB", "r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar")
  poly.PFS <- c("Gpoly", "Gnormed.poly", "U3poly", "lzpoly")
  
  # Sanity check - Are data dichotomous or polytomous?
  data.type <- NA
  if (class(try(Sanity.dma(matrix, N, I), silent=TRUE)) != "try-error")
  {
    data.type <- "dico"
    Ncat      <- 2
  } else
  {
    Ncat <- max(matrix, na.rm = TRUE) + 1
    M    <- Ncat - 1
    if (class(try(Sanity.dma.poly(matrix, N, I, M), silent=TRUE)) != "try-error")
    {
      data.type <- "poly"
    } else
    {
      stop('The data matrix is not dichotomous (0/1 scores) nor 
           polytomous (scores in {0, 1, ..., Ncat-1}, including 0 and (Ncat - 1)). Aborted.')
    }
  }
  
  # Sanity check - Were any PFS methods added?
  if (length(method) == 0)
  {
    stop('Please add your PFSs of choice to vector "method" before proceeding. Aborted.')
  }
  
  # Sanity check - Are the methods in accordance to the type of data?
  if ((data.type == "dico") & !all(method %in% dico.PFS))
  {
    stop('One or more PFSs declared in parameter "method" are not suitable to dichotomous data. Aborted.')
  } else
  {
    if (is.null(IRT.PModel)) {IRT.PModel <- "2PL"}
    if (is.null(Ability.PModel)) {Ability.PModel <- "ML"}
  }
  # 
  if ((data.type == "poly") & !all(method %in% poly.PFS))
  {
    stop('One or more PFSs declared in parameter "method" are not suitable to polyotomous data. Aborted.')
  } else
  {
    if (is.null(IRT.PModel)) {IRT.PModel <- "GRM"}
    if (is.null(Ability.PModel)) {Ability.PModel <- "EAP"}
  }
  
  # Compute PFSs:
  res <- vector("list", length(method))
  if (data.type == "dico")
  {
    for (i in 1:length(method))
    {
      res[[i]] <- eval(parse(text = method[i]))(matrix, 
                                                NA.method, Save.MatImp, 
                                                IP, IRT.PModel, Ability, Ability.PModel, mu, sigma)
    }
  }
  # 
  if (data.type == "poly")
  {
    for (i in 1:length(method))
    {
      res[[i]] <- eval(parse(text = method[i]))(matrix, Ncat, 
                                                NA.method, Save.MatImp, 
                                                IP, IRT.PModel, Ability, Ability.PModel)
    }
  }
  
  if (simplified == TRUE)
  {
    rownames.bckp <- rownames(res[[1]]$PFscores)
    res           <- data.frame(matrix(unlist(lapply(res, function(lst) {lst[[1]]})), nrow = N))
    colnames(res) <- method
    rownames(res) <- rownames.bckp
  }
  
  return(res)
}




