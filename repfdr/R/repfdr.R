repfdr <- function(pdf.binned.z, binned.z.mat, non.null = c('replication','meta-analysis','user.defined'),
                   non.null.rows = NULL, Pi.previous.result=NULL, control = em.control())
{
  if (is.null(Pi.previous.result)) {
    Pi <- piem(pdf.binned.z, binned.z.mat, control)$last.iteration
  } else {
    Pi <- Pi.previous.result
  }
  
  H  <- Pi[,-dim(Pi)[2],drop=FALSE]
  
  h0 <- switch(match.arg(non.null), replication = which(apply(H,1,function(y){ sum(y==1)<=1 & sum(y==-1)<=1 })),
               "meta-analysis" = which(rowSums(abs(H))==0),
               user.defined = (1:dim(H)[1])[-non.null.rows])
  
  if (non.null=='user.defined' & is.null(non.null.rows))
    stop("'user.defined' is selected but rows are not specified.")
  
  if (non.null!='user.defined' & !is.null(non.null.rows))
    warning(sprintf("%s is selected, supplied rows are ignored.",non.null))
  
  if (length(non.null.rows) >= dim(H)[1])
    stop("Number of selected configurations is larger than possible.")
  
  chunksize <- 20000
  if (dim(binned.z.mat)[1] <= chunksize)
  {
    chunkbegin <- 1
    chunkend   <- dim(binned.z.mat)[1]
  }
  else
  {
    chunkbegin <- c(1,seq(from = chunksize,to = dim(binned.z.mat)[1],by = chunksize))
    chunkend   <- c(chunkbegin[-1]-1,dim(binned.z.mat)[1])
  }
  
  fdr <- NULL
  for (b in 1:length(chunkbegin)) {
    fh <- ldr(pdf.binned.z, binned.z.mat[chunkbegin[b]:chunkend[b],,drop=FALSE],Pi = Pi, h.vecs = h0)
    if (dim(fh)[1]==1)
      fdr <- c(fdr,fh[,-(1:dim(H)[2])])  
    else
      fdr <- c(fdr,colSums(fh[,-(1:dim(H)[2]),drop=FALSE]))  
  }
  
  o <- order(fdr)
  ro <- order(o)
  Fdr <- (cumsum(fdr[o])/(1:length(fdr)))[ro]
    
  return (list(mat = cbind(fdr = fdr, Fdr = Fdr) ,Pi = Pi))
}

piem <- function(pdf.binned.z, binned.z.mat, control = em.control())
{
  inputchk(pdf.binned.z, binned.z.mat)
  
  n.studies <- dim(pdf.binned.z)[1]
  n.bins <- dim(pdf.binned.z)[2]
  n.association.status=dim(pdf.binned.z)[3]
  
  H <- hconfigs(n.studies,n.association.status)
  
  if (n.association.status == 2)
  {
    pbz <- array(c(array(0,dim=c(n.studies,n.bins,1)),pdf.binned.z),
                 dim=c(n.studies,n.bins,3))
  }
  
  if (n.association.status == 3)
    pbz <- pdf.binned.z
    
  if (is.null(control$pi.initial)) {
    control$pi.initial <- rep(0.1 / (dim(H)[1] - 1), dim(H)[1])
    zer.idx <- apply((H == 0), 1, all)
    control$pi.initial[zer.idx] <- 0.9
  }
  
  h     <- apply(H + 1, 2, as.integer)
  bzm <- apply(binned.z.mat - 1, 2, as.integer)
  
  out <- .Call('REM', h, pi.initial = control$pi.initial, pbz, bzm, max.iter = control$max.iter,
               tol = control$tol, nr.threads = control$nr.threads, verbose = control$verbose)
  EMi   <- out[1, control$max.iter + 1]
  out <- out[, 1:EMi,drop=FALSE]
  if (EMi == control$max.iter) {
    warning('EM did not converge within tolerance after maximum number of iterations specified.', call. = FALSE)
  } else {
    cat("Converged after",EMi,"EM iterations.")
  }
  return(list(all.iterations=cbind(H,out),
              last.iteration=cbind(H,Pi=out[,dim(out)[2]])))
}  

ldr <- function(pdf.binned.z, binned.z.mat ,Pi, h.vecs = NULL)
{
  if (is.vector(binned.z.mat)) #when only one test is selected
    binned.z.mat <- matrix(binned.z.mat,nrow=1,ncol=length(binned.z.mat))
  
  inputchk(pdf.binned.z, binned.z.mat)
  
  n.studies <- dim(pdf.binned.z)[1]
  n.bins <- dim(pdf.binned.z)[2]
  n.association.status <- dim(pdf.binned.z)[3]
  
  if (n.association.status == 2)
  {
    pbz <- array(c(array(0,dim=c(n.studies,n.bins,1)),pdf.binned.z),
                 dim=c(n.studies,n.bins,3))
  }
  
  if (n.association.status == 3)
    pbz <- pdf.binned.z
    
  H  <- Pi[,-dim(Pi)[2],drop=FALSE]
  nrowH <- dim(H)[1]
  ncolH <- dim(H)[2]
  
  if(ncolH!=n.studies)
    stop('First dimension of pdf.binned.z do not match the number of studies indicated by Pi.')
  
  nrowbz <- dim(binned.z.mat)[1]
    
  if (is.null(h.vecs))
    h.vecs = 1:dim(H)[1] else
      if (max(h.vecs) > nrowH)
        stop("Number of configurations to report is larger than possible.")
  
  fzh <- matrix(nrow = nrowH, ncol = nrowbz)
  for(i in 1:nrowH)
  {
    term <- vector(length = nrowbz)
    term[] <- 1
    for (j in 1:ncolH)
      term <- term*pbz[j,binned.z.mat[,j],H[i,j]+2]
    fzh[i,] <- term
  }
  
  Pih0     <- Pi[,dim(Pi)[2]]
  Pih0[-h.vecs] <- 0
  
  fh   <- matrix(t(t(fzh * Pih0) / as.numeric(Pi[,dim(Pi)[2]] %*% fzh ))[h.vecs,],
                 nrow=length(h.vecs),ncol=dim(binned.z.mat)[1])
  colnames(fh) <- rownames(binned.z.mat)
  return(cbind(H[h.vecs,,drop=FALSE],fh))
}

hconfigs <- function(n.studies,n.association.status = 3,studies.names=NULL)
{
  if (n.association.status == 3) Hstates <- c(0,-1,1)
  if (n.association.status == 2) Hstates <- c(0,1) 
  
  H <- as.matrix(expand.grid(replicate(n.studies,Hstates,FALSE)))
  
  if (is.null(studies.names))
    colnames(H) <- sprintf("Study %i",1:n.studies) else
    colnames(H) <- studies.names
  return(H)
}

inputchk <- function(pdf.binned.z, binned.z.mat){
  
  # pdf.binned.z dimension restriction
  if (!(dim(pdf.binned.z)[3] %in% c(2,3)))
    stop('dim(pdf.binned.z)[3] should be equal to 2 or 3')
  
  # binned.z.mat value restriction (not NA or NaN):
  if (any(is.na(binned.z.mat)) | any(!is.finite(binned.z.mat)))
      stop('binned.z.mat should contain positive integers only.')
  
  # binned.z.mat value restriction (positive integer):
  if (any(!is.int(binned.z.mat)) | !any(as.logical(binned.z.mat)))
    stop('binned.z.mat should contain positive integers only.')
  
  # binned.z.mat value restriction (maximum bin number):
  if (max(binned.z.mat) > dim(pdf.binned.z)[2])
    stop('bin number(s) in \'binned.z.mat\' is out of bin\'s range specified in \'pdf.binned.z.\'')
  
  # binned.z.mat & pdf.binned.z march in number of studies:
  if (dim(binned.z.mat)[2] != dim(pdf.binned.z)[1])
    stop('number of studies in pdf.binned.z is different than in binned.z.mat')
  
  # pdf.binned.z value restriction (0<=p<=1):
  if (any(pdf.binned.z<0 & pdf.binned.z>1))
    stop('pdf.binned.z should contain values between 0 and 1.')
  
  return(NULL)
}

is.int <- function(x, tol = .Machine$double.eps^0.5)  (x - round(x)) < tol

em.control <- function(pi.initial = NULL, max.iter = 1e4, tol = 1e-12, nr.threads = 0, verbose = TRUE)
{
  if (!is.null(pi.initial)) if (any(pi.initial<0 | pi.initial>1))
    stop("pi.initial values should be between 0 and 1.")

  if (!is.int(max.iter))
    stop('max.iter should be integer.')
  
  if (!is.int(nr.threads))
    stop('nr.threads should be integer.')
  
  if (tol <= 0)
    stop('tol should be positive number.')
  
  max.iter <- as.integer(max.iter)
  nr.threads = as.integer(nr.threads)
  verbose = as.integer(verbose)
  return (list(pi.initial = pi.initial, max.iter = max.iter, tol = tol, nr.threads = nr.threads, verbose = verbose))
}