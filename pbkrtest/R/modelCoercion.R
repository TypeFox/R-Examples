
.restrictedModelMatrix<-function(B,L) {
    ##cat("B:\n"); print(B); cat("L:\n"); print(L)
    ## find A such that <A>={Bb| b in Lb=0}

    ## if (!is.matrix(L))
    ##     L <- matrix(L, nrow=1)

    if ( !inherits(L, c("matrix", "Matrix")) )
        L <- matrix(L, nrow=1)
    L <- as(L, "matrix")
    if ( ncol(B) != ncol(L) ) {
        print(c( ncol(B), ncol(L) ))
        stop('Number of columns of B and L not equal \n')
    }
    A <- B %*% .orthComplement(t(L))
    A
}

.restrictionMatrixBA<-function(B,A) {
  ## <A> in <B>
  ## determine L such that  <A>={Bb| b in Lb=0}
  d <- rankMatrix(cbind(A,B)) - rankMatrix(B)
  if (d > 0) {
    stop('Error:  <A> not subspace of <B> \n')
  }
  Q  <- qr.Q(qr(cbind(A,B)))
  Q2 <- Q[,(rankMatrix(A)+1):rankMatrix(B)]
  L  <- t(Q2)  %*% B
  ##make rows of L2 orthogonal
  L <-t(qr.Q(qr(t(L))))
  L
}

.model2restrictionMatrix <- function (largeModel, smallModel) {
  L <- if(is.matrix(smallModel)) {
    ## ensures  that L is of full row rank:
    LL <- smallModel
    q  <- rankMatrix(LL)
    if (q < nrow(LL) ){
      t(qr.Q(qr(t(LL)))[,1:qr(LL)$rank])
    } else {
      smallModel
    }
  } else  { #smallModel is mer model
    .restrictionMatrixBA(getME(largeModel,'X'),  getME(smallModel,'X'))
  }
  L<-.makeSparse(L)
  L
}


model2restrictionMatrix <- function (largeModel, smallModel) {
    UseMethod("model2restrictionMatrix")
}


model2restrictionMatrix.merMod <-
    model2restrictionMatrix.mer <-
        function (largeModel, smallModel) {
        L <- if(is.matrix(smallModel)) {
            ## ensures  that L is of full row rank:
            LL <- smallModel
            q  <- rankMatrix(LL)
            if (q < nrow(LL) ){
                t(qr.Q(qr(t(LL)))[,1:qr(LL)$rank])
            } else {
                smallModel
            }
        } else  { #smallModel is mer model
            .restrictionMatrixBA(getME(largeModel,'X'),  getME(smallModel,'X'))
        }
        L<-.makeSparse(L)
        L
    }

model2restrictionMatrix.lm <- function (largeModel, smallModel) {
  L <- if(is.matrix(smallModel)) {
    ## ensures  that L is of full row rank:
    LL <- smallModel
    q  <- rankMatrix(LL)
    if (q < nrow(LL) ){
      t(qr.Q(qr(t(LL)))[,1:qr(LL)$rank])
    } else {
      smallModel
    }
  } else  { #smallModel is mer model
    .restrictionMatrixBA(model.matrix( largeModel ),  model.matrix( smallModel ))
  }
  L<-.makeSparse(L)
  L
}








.formula2list <- function(form){
  lhs <- form[[2]]
  tt  <- terms(form)
  tl  <- attr(tt, "term.labels")
  r.idx <- grep("\\|", tl)

  if (length(r.idx)){
    rane  <- paste("(", tl[r.idx], ")")
    f.idx <- (1:length(tl))[-r.idx]
    if (length(f.idx))
      fixe  <- tl[f.idx]
    else
      fixe  <- NULL
  } else {
    rane <- NULL
    fixe <- tl
  }

  ans <- list(lhs=deparse(lhs),
              rhs.fix=fixe,
              rhs.ran=rane)
  ans
}

restrictionMatrix2model <- function(largeModel, LL){
  UseMethod("restrictionMatrix2model")
}

restrictionMatrix2model.merMod <-
    restrictionMatrix2model.mer <-
  function(largeModel, LL){

  XX.lg 	 <- getME(largeModel, "X")

  form <- as.formula(formula(largeModel))
  attributes(XX.lg)[-1] <- NULL
  XX.sm <- .restrictedModelMatrix(XX.lg, LL)

  ncX.sm  <- ncol(XX.sm)
  colnames(XX.sm) <- paste(".X", 1:ncX.sm, sep='')

  rhs.fix2 <- paste(".X", 1:ncX.sm, sep='', collapse="+")

  fff  <- .formula2list(form)
  new.formula <- as.formula(paste(fff$lhs, "~ -1+", rhs.fix2, "+", fff$rhs.ran))
  new.data    <- cbind(XX.sm, eval(largeModel@call$data))

##   ans <- lmer(eval(new.formula), data=new.data, REML=getME(largeModel, "is_REML"))
  ans <- update(largeModel, eval(new.formula), data=new.data)
  ans
}

restrictionMatrix2model.lm <- function(largeModel, LL){

  form <- as.formula(formula(largeModel))
  XX.lg 	 <- model.matrix(largeModel)
  attributes(XX.lg)[-1] <- NULL
  XX.sm <- zapsmall( .restrictedModelMatrix(XX.lg, LL) )

  ncX.sm  <- ncol(XX.sm)
  colnames(XX.sm) <- paste(".X", 1:ncX.sm, sep='')

  rhs.fix2 <- paste(".X", 1:ncX.sm, sep='', collapse="+")
  fff  <- .formula2list(form)
  new.formula <- as.formula(paste(fff$lhs, "~ -1+", rhs.fix2))
  new.data    <- as.data.frame(cbind(XX.sm, eval(largeModel$model)))
  #print(new.data)
  ans <- update(largeModel, eval(new.formula), data=new.data)
  ans
}


