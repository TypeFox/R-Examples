mosek_lptoprob =
function(f=NA, A=NA, b=NA, Aeq=NA, beq=NA, lb=NA, ub=NA) {

  #
  # Input validation
  #

  stopifnot(all(!is.na(f)));

  if (all(!is.na(A)) || all(!is.na(b))) {
    stopifnot(all(!is.na(A)) && all(!is.na(b)));

    # Convert to a supported matrix format (COO or CSC)
    if (!is(A, 'TsparseMatrix')) {
      A <- as(A, 'CsparseMatrix');
    }

    stopifnot(nrow(A) == length(b));
    stopifnot(ncol(A) == length(f));

  } else {
    A <- Matrix(0, nrow=0, ncol=length(f), sparse=TRUE);
    b <- numeric(0);
  }

  if (all(!is.na(Aeq)) || all(!is.na(beq))) {
    stopifnot(all(!is.na(Aeq)) && all(!is.na(beq)));

    # Convert to a supported matrix format (COO or CSC)
    if (!is(Aeq, 'TsparseMatrix')) {
      Aeq <- as(Aeq, 'CsparseMatrix');
    }

    stopifnot(nrow(Aeq) == length(beq));
    stopifnot(ncol(Aeq) == length(f));
  } else {
    Aeq <- Matrix(0, nrow=0, ncol=length(f), sparse=TRUE);
    beq <- numeric(0);
  }

  stopifnot(nrow(A) + nrow(Aeq) >= 1);

  stopifnot(all(!is.na(lb)));
  stopifnot(length(lb) == length(f));

  stopifnot(all(!is.na(ub)));
  stopifnot(length(ub) == length(f));


  #
  # Transformation
  #

  prob <- list(sense = "min");
  prob$c <- f;

  nrA <- nrow(A);
  prob$A <- rBind(A, Aeq);

  prob$bc <- rbind(blc = c(rep(-Inf,nrA),  beq),
                   buc = c(            b,  beq));

  prob$bx <- rbind(lb, ub);

  return(prob);
}



mosek_qptoprob =
function(F=NA, f=NA, A=NA, b=NA, Aeq=NA, beq=NA, lb=NA, ub=NA) {

  #
  # Input validation
  #

  stopifnot(all(!is.na(F)));

  stopifnot(all(!is.na(f)));
  stopifnot(length(f) == ncol(F));

  if (all(!is.na(A)) || all(!is.na(b))) {
    stopifnot(all(!is.na(A)) && all(!is.na(b)));

    # Convert to a supported matrix format (COO or CSC)
    if (!is(A, 'TsparseMatrix')) {
      A <- as(A, 'CsparseMatrix');
    }

    stopifnot(nrow(A) == length(b));
    stopifnot(ncol(A) == length(f));
  } else {
    A <- Matrix(0, nrow=0, ncol=length(f), sparse=TRUE);
    b <- numeric(0);
  }

  if (all(!is.na(Aeq)) || all(!is.na(beq))) {
    stopifnot(all(!is.na(Aeq)) && all(!is.na(beq)));

    # Convert to a supported matrix format (COO or CSC)
    if (!is(Aeq, 'TsparseMatrix')) {
      Aeq <- as(Aeq, 'CsparseMatrix');
    }

    stopifnot(nrow(Aeq) == length(beq));
    stopifnot(ncol(Aeq) == length(f));
  } else {
    Aeq <- Matrix(0, nrow=0, ncol=length(f), sparse=TRUE);
    beq <- numeric(0);
  }

  stopifnot(all(!is.na(lb)));
  stopifnot(length(lb) == length(f));

  stopifnot(all(!is.na(ub)));
  stopifnot(length(ub) == length(f));


  #
  # Transformation
  #

  prob <- list(sense = "min");

  nt <- nrow(F);
  nx <- ncol(F);

  nrA <- nrow(A);
  nrEQ <- nrow(Aeq);

  #                    x                 v                 w               t  
  # ----------------------------------------------------------------------------
  prob$c <- c(         f,                1,                0,         rep(0,nt) );
  prob$A <- rBind(
              cBind(   A,  Matrix(0,nrA,1),  Matrix(0,nrA,1),  Matrix(0,nrA,nt) ),
              cBind( Aeq, Matrix(0,nrEQ,1), Matrix(0,nrEQ,1), Matrix(0,nrEQ,nt) ),
              cBind(   F,  Matrix(0, nt,1),  Matrix(0, nt,1),   -1*Diagonal(nt) )
            );


  #                                  con1  con2      con3
  # -------------------------------------------
  prob$bc <- rbind(blc = c(rep(-Inf,nrA),   beq, rep(0,nt)),
                   buc = c(            b,   beq, rep(0,nt)));


  #                          x    v  w            t
  # ---------------------------------------
  prob$bx <- rbind(blx = c( lb,   0, 1, rep(-Inf,nt) ),
                   bux = c( ub, Inf, 1, rep( Inf,nt) ));


  prob$cones <- matrix(nrow=2, dimnames=list(c("type","sub"),c()),
                  list( "RQUAD", nx+(1:(2+nt)) ),
                );

  return(prob);
}

