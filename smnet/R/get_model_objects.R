
get_model_objects<-function(formula, data, adjacency, 
                                 response.locs, formulaout,
                                 control, ord = NULL, rid_data = NULL, 
                                 netID = NULL, weight = NULL){
  
  # INTEPRET THE FORMULA, FIND THE NAMES AND NUMBER OF DIFFERENT MODEL COMPONENTS
  # ------------------------------------------------------------------------------
  adj         <- adjacency$adjacency
  n           <- nrow(data)
  net         <- !is.null(adj)
  n.segments  <- ifelse(net, nrow(adj), "NULL")
  n.terms     <- length(attr(terms(formulaout$formula), "variables")) - 2
  term.names  <- attr(terms(formulaout$formula), "term.labels")
  st3         <- substr(term.names, 1,3)
  st2         <- substr(term.names, 1,2)
  sm.names    <- term.names[st2 == "m("]
  net.names   <- term.names[st3 == "net"]
#   lin.names   <- term.names[!((st2 == "m(") | (st3 == "net"))]
  n.smooth    <- length(sm.names)
  n.net       <- length(net.names)
#   n.linear    <- length(lin.names)
  
  # need to obtain the number of columns in the model matrix that are simple linear ones
  # do this by working out the number associated with smooth/ network terms and subtract
  # this from number of columns in the data mat in formulaout$mf
  smooth.variables  <- sum(unlist(lapply(formulaout$gp$smooth.spec, function(L) length(L$term))))
  network.variables <- sum(unlist(lapply(formulaout$gp$smooth.spec, function(L) !is.null(L$weight))))
  n.linear          <- ncol(as.matrix(formulaout$mf)) - smooth.variables - 1
  if(n.linear > 0) lin.names         <- colnames(as.matrix(formulaout$mf))[2:(n.linear + 1)] else lin.names <- NULL
  n.terms           <- n.linear + n.smooth + n.net

  variables   <- as.matrix(formulaout$mf[,-1])
  response    <- formulaout$mf[,1]
  var.names   <- colnames(as.matrix(formulaout$mf[1,]))[-1]
  colnames(variables) <- var.names
  sm.terms.names <- vector("list")

  
  # POPULATE A LIST CONTAINING THE DIFFERENT PARTS OF THE DESIGN MATRIX
  # ------------------------------------------------------------------------------
  X.list  <-  vector("list", length = (n.terms+1))
  # intercept first - currently included by default
  X.list[[1]]  <- spam(1, ncol = 1, nrow = n)
  # subsequent n.linear columns are the linear covariates
  # remember to remove the means of each linear variable
  if(n.linear > 0){
    lin.means <- vector("numeric", length = n.linear)
    for(i in 1:n.linear){
      lin.means[i]  <- 0#mean(variables[,i])
      variables[,i] <- variables[,i]# - lin.means[i]
      X.list[[i+1]] <- as.spam(variables[,i])
    }
  } else {lin.means = NULL}

  # the next components of X.list are the univariate and bivariate smooth terms
  if(n.smooth>0){
    sm.basis <- varbs.len <- sm.cyclic <- vector("numeric", length=n.smooth)
    varbs.loc<- vector("list", length=n.smooth)
    for(i in 1:n.smooth){
      if(net){
        smooth.only <- formulaout$gp$smooth.spec[-which(unlist(lapply(formulaout$gp$smooth.spec, class)) == "network.spec")]
      } else smooth.only <- formulaout$gp$smooth.spec
      # interpret the user specified basis, if supplied
      user_basis  <- as.numeric(smooth.only[[i]]$bs.dim)
      sm.cyclic[i]<- smooth.only[[i]]$cyclic
      sm.basis[i] <- ifelse(user_basis == -1, 10, user_basis)
      # this is the name of the i^th smooth term
      varbs.char     <- smooth.only[[i]]$term
      varbs.loc[[i]] <- match(varbs.char, var.names)
      varbs.len[i]   <- length(varbs.loc[[i]])
      varbs          <- variables[,varbs.loc[[i]]]
      sm.terms.names[[i]] <- varbs.char
      # if univariate construct a 1-D B-spline basis
      if(varbs.len[i] == 1) X.list[[i+1+n.linear]]<-b_spline_basis(x=varbs, nseg = (sm.basis[i] - 3))#, 
#                                                                    cyclic = sm.cyclic[i], deg=3)
      # if bivariate construct a tensor B-spline basis - current higher order not supported
      if(varbs.len[i] == 2){
        a1<-b_spline_basis(x=varbs[,1], nseg = (sm.basis[i]-3), deg=3)#, cyclic = (sm.cyclic[i] == 1)|(sm.cyclic[i] == T))
        a2<-b_spline_basis(x=varbs[,2], nseg = (sm.basis[i]-3), deg=3)#, cyclic = (sm.cyclic[i] == 2)|(sm.cyclic[i] == T))
        X.list[[i + 1 + n.linear]] <- b1<-make_spam(not_sparse_box_product(a1, a2))
        constraint<-crossprodspam(b1)
      }
    }
  }
  
  # if specified, the network matrix is added as the final set of columns in the model matrix
  if(net){
    # sometimes the rid starts at 0?  in which case add one to each identifier
    add.one<-ifelse(min(as.numeric(response.locs)) == 0, 1, 0)
    # spatial component; currently hardcoded to the last of X.list
    X.list[[n.linear+n.smooth+2]]<-spam(
      list(i = 1:n, j = (as.numeric(response.locs)+add.one),  rep(1, n)), 
      nrow = n, ncol = n.segments)
    # this identity matrix is going to be a ridge penalty for the network component
    networkStabiliser<-diag.spam(1, n.segments)
  }
  
  # summaries of X.list, convert to a sparse spam object, concatenate and crossprod
  X.dim          <- lapply(X.list, ncol)
  X.list         <- lapply(X.list, make_spam)
  X.spam         <- Reduce("cbind", X.list)
  XTX.spam       <- info <- t(X.spam)%*%X.spam
  
  # POPULATE A LIST CONTAINING THE DIFFERENT PENALTY MATRICES
  # IF I^TH COMPONENT HAS MORE THAN ONE PENALTY, THESE ARE CONTAINED IN A
  # SUBLIST WHICH IS THE I^TH ELEMENT OF P.LIST
  # ------------------------------------------------------------------------------
  # each penalty matrix needs to be n*n where n is the number of colums in X.list
  # easist to do this by first creating a list of zero matrices 
  if(n.smooth>0 | net) blockzero<-lapply(X.dim, make_sparse)
  # output lists.  Pwee stores only the individual penalties submatrices, P.list stores
  # full n*n matrices, *containing* the submatrix in Pwee
  P.list <- Pwee.list <- vector("list")
  # first start with univariate smooth terms
  if(n.smooth>0){
    for(i in 1:n.smooth){
      if(varbs.len[i] == 1){
        P         <- get_penalty(sm.basis[i], cyclic = F)#sm.cyclic[i])
        Pwee.list <- c(Pwee.list, list(P=P))
        P.list    <- c(P.list, get_block_penalty(P, blockzero, (i+1+n.linear)))
      }
      if(varbs.len[i] == 2){
        cyclic_1 <- F# if((sm.cyclic == 1)|(sm.cyclic == T)) cyclic_1 <- T else cyclic_1 <- F
        cyclic_2 <- F#if((sm.cyclic == 2)|(sm.cyclic == T)) cyclic_2 <- T else cyclic_2 <- F
        Pmat1     <- get_penalty(sm.basis[i], cyclic = cyclic_1)
        Pmat2     <- get_penalty(sm.basis[i], cyclic = cyclic_2)
        Id        <- diag.spam(1, sm.basis[i]) 
        Pwee1     <- kronecker.spam(Pmat1,  Id)
        Pwee2     <- kronecker.spam(Id, Pmat2)
        P1        <- get_block_penalty(Pwee1, blockzero, (i+1+n.linear))
        P2        <- get_block_penalty(Pwee2, blockzero, (i+1+n.linear))
        Pwee.list <- c(Pwee.list, list(Pw1=Pwee1 , Pw2=Pwee2))
        P.list[[length(P.list)+1]]<-list(P1 , P2)
      }
    }
    P.flat<-make_flat(P.list)
  }
  

  # CONSTRUCT THE SPATIAL NETWORK PENALTY, THE MOST COMPLEX PART HERE
  # THIS INVOLVES CAREFULLY POPULATING THE MATRIX USING THE ADJACENCY MATRIX 
  # AND THE VECTOR OF WEIGHTS SPECIFIED BY THE USER THAT IS CONTAINED IN THE DATA
  # ------------------------------------------------------------------------------
  if(net){
    adj.spam <- make_spam(adj)
    # create 2 penlaties, one across confluences
    pseudo.inds  <- which(colSums.spam(adj.spam) == 1)
    ij.nzero.adj <- triplet(adj.spam)$indices
    in.pseudo    <- ij.nzero.adj[,2] %in% pseudo.inds
    ij.confl     <- ij.nzero.adj#[!in.pseudo,]
    n.nzero      <- nrow(ij.confl)
    p.row.ind    <- rep(1:n.nzero, each = 2)
    p.col.ind    <- c(t(ij.confl))
    p.val        <- weight[rep(ij.confl[,1], each = 2)]*rep(c(-1, 1), n.nzero)
    D2           <- spam(list(i=p.row.ind, j=p.col.ind, p.val), nrow = n.nzero, ncol = n.segments)
    D2           <- crossprodspam(D2)
    P2           <- get_block_penalty(D2, blockzero, i = n.terms + 1)
    # put these matrices into a single object by growing two lists
    P.seg        <- Pwee <- list()
    Pwee         <- c(Pwee, D2 = D2, ID = networkStabiliser)
    Pwee.list    <- c(Pwee.list, list(Pwee))
    P.seg        <- c(P.seg, P2 = P2, get_block_penalty(networkStabiliser, blockzero, i = n.terms + 1))
    P.list       <- c(P.list, list(P.seg))
  }
  
  output<-list(X.list = X.list, X.spam = X.spam, XTX.spam = XTX.spam, 
               variables = variables, adjacency = adjacency, P.list = P.list,
               Pwee.list = Pwee.list, response = response, net = net, n.linear = n.linear, 
               n.smooth = n.smooth, lin.names = lin.names, 
               sm.names = sm.names, control = control, 
               ord = ord, rid_data = rid_data, netID = netID, 
               weight = weight, lin.means = lin.means)
  if((n.smooth)>0){
    output<-c(output, list(P.list = P.list, sm.basis = sm.basis, 
                           varbs.len = varbs.len, varbs.loc = varbs.loc, 
                           sm.terms.names=sm.terms.names, sm.cyclic = sm.cyclic))
  }
  output
}

