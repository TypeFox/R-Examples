coerce_mat <- function(mat,cfun=as.simple_triplet_sym_matrix) if (is.null(dim(mat))) mat else cfun(mat)
coerce_blkmat <- function(Aij,cfun=as.simple_triplet_sym_matrix) lapply(Aij,coerce_mat,cfun=cfun)
coerce_const <- function(Ai,cfun=as.simple_triplet_sym_matrix) lapply(Ai,coerce_blkmat,cfun=cfun)

getdf_mat <- function(x)
  {
    if(!is.null(dim(x))) {
      x <- as.simple_triplet_sym_matrix(x)
      if (length(x$v) == 0)
        return(NULL)
      else
        return(cbind(x$j,x$i,x$v))
    }
    else {
      ind <- which(x != 0)
      if (length(ind) == 0)
        return(NULL)
      else
        return(cbind(ind,ind,x[ind]))
    }
  }

getdf_blkmat <- function(x)
  {
    do.call(rbind,lapply(1:length(x),function(j) {df <- getdf_mat(x[[j]]); if (is.null(df)) return(NULL) else return(cbind(j,df))}))
  }

getdf_const <- function(x)
  {
    do.call(rbind,lapply(1:length(x), function(i) {df <- getdf_blkmat(x[[i]]); if (is.null(df)) return(NULL) else return(cbind(i,df))}))
  }

readsdpa <- function(file="",verbose=FALSE)
  {
    if (file=="")
      stop("'file' argument must be a non-empty string")
    
    ret <- .Call("readsdpa",
                 as.character(file),
                 as.integer(verbose),
                 PACKAGE="Rcsdp")

    names(ret) <- c("C","A","b","K");
    names(ret$K) <- c("type","size");
    ret$K$type <- c("s","l")[vector_csdp2R(ret$K$type)];
    ret$K$size <- vector_csdp2R(ret$K$size);

    m <- length(ret$b)-1;
    prob.info <- get.prob.info(ret$K,m)
    ret$C <- blkmatrix_csdp2R(ret$C,prob.info);
    ret$A <- constraints_csdp2R(ret$A,prob.info);
    ret$b <- vector_csdp2R(ret$b);
    ret
  }

readsdpa.sol <- function(K,C,m,file="")
  {
    if (file=="")
      stop("'file' argument must be a non-empty string")

    prob.info <- get.prob.info(K,m);
    ret <- .Call("readsdpa_sol",
                 as.character(file),
                 as.integer(sum(prob.info$block.sizes)),
                 as.integer(m),
                 blkmatrix_R2csdp(C,prob.info),
                 PACKAGE="Rcsdp");
    names(ret) <- c("X","y","Z");
    ret$X <- blkmatrix_csdp2R(ret$X,prob.info);
    ret$y <- vector_csdp2R(ret$y);
    ret$Z <- blkmatrix_csdp2R(ret$Z,prob.info);
    ret
  }

writesdpa <- function(C,A,b,K,file="")
  {
    if (file=="")
      stop("'file' argument must be a non-empty string")

    prob.info <- get.prob.info(K,length(b));
    validate.data(C,A,b,prob.info)
    prob.data <- prepare.data(C,A,b,prob.info)
    .Call("writesdpa",
          as.character(file),
          as.integer(sum(prob.info$block.sizes)),
          as.integer(prob.info$nconstraints),
          as.integer(prob.info$nblocks),
          as.integer(c(0,prob.info$block.types)),
          as.integer(c(0,prob.info$block.sizes)),
          prob.data$C,
          prob.data$A,
          prob.data$b,
          PACKAGE="Rcsdp")
  }

writesdpa.sol <- function(X,Z,y,K,file="")
  {
    if (file=="")
      stop("'file' argument must be a non-empty string")

    prob.info <- get.prob.info(K,length(y));

    .Call("writesdpa_sol",
          as.character(file),
          as.integer(sum(prob.info$block.sizes)),
          as.integer(prob.info$nconstraints),
          blkmatrix_R2csdp(X,prob.info),
          vector_R2csdp(y),
          blkmatrix_R2csdp(Z,prob.info),
          PACKAGE="Rcsdp")
  }


