DiscML <- function (x, phy, CI = FALSE, model = "ER", reversible =FALSE,  zerocorrection=FALSE, 
                    ip = 0.1, rootprobability=FALSE, irootprobability=FALSE,  alpha=FALSE, ialpha= 0.5, ngamma = 8,kappa = 1,
                    characters= TRUE,
                    simplify = FALSE ,individualrates = FALSE, plotmu = FALSE,  plotloglik = FALSE ) 
{ 
  start <- proc.time()
  #The codes inside midpoint function were taken from phangorn package
  midpoint<- function(tree)
  {
    allAncestors <- function (x) 
    {
      x = reorder(x, "postorder")
      parents <- x$edge[, 1]
      child <- x$edge[, 2]
      l = length(parents)
      res <- vector("list", max(x$edge))
      for (i in l:1) {
        pa = parents[i]
        res[[child[i]]] = c(pa, res[[pa]])
      }
      res
    }
    reroot <- function (tree, node) 
    {
      anc = Ancestors(tree, node, "all")
      l = length(anc)
      if (is.na(match(node, tree$edge[, 1]))) 
        stop("node not in tree")
      if (l == 0) 
        return(tree)
      ind = match(c(node, anc[-l]), tree$edge[, 2])
      tree$edge[ind, c(1, 2)] = tree$edge[ind, c(2, 1)]
      root = anc[l]
      tree$edge[tree$edge == root] = 0L
      tree$edge[tree$edge == node] = root
      tree$edge[tree$edge == 0L] = node
      tree <- collapse.singles(tree)
      attr(tree, "order") <- NULL
      reorder(tree, "postorder")
    }
    
    Ancestors<-  function (x, node, type = c("all", "parent")) 
    {
      parents <- x$edge[, 1]
      child <- x$edge[, 2]
      pvector <- numeric(max(x$edge))
      pvector[child] <- parents
      type <- match.arg(type)
      if (type == "parent") 
        return(pvector[node])
      anc <- function(pvector, node) {
        res <- numeric(0)
        repeat {
          anc <- pvector[node]
          if (anc == 0) 
            break
          res <- c(res, anc)
          node <- anc
        }
        res
      }
      if (length(node) == 1) 
        return(anc(pvector, node))
      else allAncestors(x)[node]
    }
    node2root <- function(x) {
      x = reorder(x, "postorder")
      el = numeric(max(x$edge))
      parents <- x$edge[, 1]
      child <- x$edge[, 2]
      el[child] = x$edge.length
      l = length(parents)
      res <- numeric(max(x$edge))
      for (i in l:1) {
        res[child[i]] = el[child[i]] + res[parents[i]]
      }
      res
    }
    #The codes below were taken from 'ace' function from 'ape' package,
    #and were modified for the purpose of DiscML.
    tree = unroot(tree)
    nTips = length(tree$tip)
    maxD1 = node2root(tree)[1:nTips]
    ind = which.max(maxD1)
    tmproot = Ancestors(tree, ind, "parent")
    tree = reroot(tree, tmproot)
    el = numeric(max(tree$edge))
    el[tree$edge[, 2]] = tree$edge.length
    maxdm = el[ind]
    tree$edge.length[tree$edge[, 2] == ind] = 0
    maxD1 = node2root(tree)[1:nTips]
    tree$edge.length[tree$edge[, 2] == ind] = maxdm
    ind = c(ind, which.max(maxD1))
    maxdm = maxdm + maxD1[ind[2]]
    rn = max(tree$edge) + 1
    edge = tree$edge
    el = tree$edge.length
    children = tree$edge[, 2]
    left = match(ind[1], children)
    tmp = Ancestors(tree, ind[2], "all")
    tmp = c(ind[2], tmp[-length(tmp)])
    right = match(tmp, children)
    if (el[left] >= (maxdm/2)) {
      edge = rbind(edge, c(rn, ind[1]))
      edge[left, 2] = rn
      el[left] = el[left] - (maxdm/2)
      el = c(el, maxdm/2)
    }
    else {
      sel = cumsum(el[right])
      i = which(sel > (maxdm/2))[1]
      edge = rbind(edge, c(rn, tmp[i]))
      edge[right[i], 2] = rn
      eltmp = sel[i] - (maxdm/2)
      el = c(el, el[right[i]] - eltmp)
      el[right[i]] = eltmp
    }
    tree$edge.length = el
    tree$edge = edge
    tree$Nnode = tree$Nnode + 1
    attr(tree, "order") <- NULL
    reorder(reroot(tree, rn), "postorder")
  }
  rateparameters <- p <- TRUE
  mu <- TRUE
  imu <- 1
  if(is.character(phy))
    phy <- read.tree2(text=phy)
  if(!is.rooted (phy) )
  {
    warning("Warning: The tree in the input is an unrooted tree. The tree will be midpoint rooted.", call. = TRUE)
    phy<- midpoint(phy)
    if(!is.null(phy$dollarData))
      phy$dollarData<-midpoint(phy$dollarData)
  }
  if( !is.null(phy$dollarData))
  {
    muphy <- reorder(phy$dollarData, "postorder") 
    tempchars<- phy$tempchars
    tempvalues <- phy$tempvalues
    nmuphy<- length(levels(factor(phy$dollarData$edge.length)))
  }
  else
  {
    tempchars <- "mu0"
    tempvalues <- NA
    muphy<- list()
    muphy$edge.length = 1
    nmuphy <- 1
  }
  if(!identical(plotmu, FALSE) && identical(individualrates, FALSE))
    stop("You can set either 'plotmu=TRUE' or 'plotmu= \"a string\"' only when 'individualrates = TRUE'")
  if(!identical(plotmu, FALSE) && identical(mu, FALSE))
    stop("You have to set 'mu = TRUE' (or a string) in order to set 'plotmu = TRUE.'")
  if(!identical(plotloglik, FALSE) && identical(individualrates, FALSE))
    stop("You can set 'plotloglik= TRUE' only when 'individualrates = TRUE'") 
  count <- 0 
  talphacount <-0
  startalpha <-0
  multiplier <-1
  if(identical(plotmu, TRUE))
    mu <-TRUE
  cnames <- NULL
  if(is.vector(x))
    cnames <- names(x)
  else
    cnames <- colnames(x)
  if(is.data.frame(x))
  {
    tempnames <- cnames
    x <- suppressWarnings(as.matrix(x))
    if(is.character(x))
    {   
      tempx <- x
      tx<- suppressWarnings(matrix(as.integer(tempx), nrow(tempx), ncol(tempx)))
      if(is.na(sum(tx[1,]) ))
        tempnames <- x[1,]          
      x<- tempx[!is.na(rowSums(tx)),] 
    } 
    cnames<- tempnames
    TIPS <-1:ncol(x)
    NROW <-nrow(x)
    NCOL <-ncol(x)
  }
  else if(is.matrix(x))
  {
    tempnames <- cnames
    if(is.character(x))
    {   
      tempx <- x
      tx<- suppressWarnings(matrix(as.integer(tempx), nrow(tempx), ncol(tempx)))
      if(is.na(sum(tx[1,]) ))
        tempnames <- x[1,]
      x<- tempx[!is.na(rowSums(tx)),] 
    } 
    cnames<- tempnames
    TIPS <-1:ncol(x)
    NROW <-nrow(x)
    NCOL <-ncol(x)
  }
  if(identical(simplify, TRUE))
  {
    if(is.matrix(x))
      x<- matrix( as.numeric(!(x==0)), nrow(x), ncol(x) )
    if(is.vector(x))
      x<-as.numeric(!(x==0))
  }
  if(is.matrix(x))
  { 
    colnames(x) <- cnames
    tempx<- x
    for(i in 1:length(cnames))
      tempx[,which( toupper(phy$tip.label) == toupper(colnames(x)[i]))] <- x[,i]
    x<-matrix(as.integer(tempx), nrow(x), ncol(x))
    colnames(x) <- phy$tip.label
  }
  if(is.vector(x))
  {
    names(x) <- cnames
    tempx<- x
    for(i in 1:length(cnames))
    {
      tempx[which(   toupper(phy$tip.label) == toupper((names(x)[i]))) ] <- x[i]
    }
    x<- as.numeric(tempx)
    TIPS <-1:length(x)
    NROW <-1
    NCOL <-length(x)
    names(x) <- phy$tip.label
  }
  
  if(is.matrix(x))
  {
    if (is.numeric(characters))
    {
      characters <- sort( characters)
      if (!is.factor(characters)) 
        characters2 <- factor(characters)
      lvls <- levels(characters2)
      if(!all(x %in% as.integer(lvls))) 
        stop("The argument 'characters' should contain all possible discrete integer characters.")
      nl <- nlevels(characters2) 
    }
    if(identical(characters,TRUE))
    {
      if(!is.factor(x))
        kk <- factor(x)
      else
        kk <- x      
      lvls <- levels(kk)
      nl <- nlevels(kk)
    }
  }
  if(is.vector(x))
  {
    if (is.numeric(characters))
    {
      characters <- sort( characters)
      if (!is.factor(characters)) 
        kk<- factor(characters)
      else
        kk<- characters
      lvls <- levels(kk)
      if(!all(x %in% as.integer(lvls))) 
        stop("The argument 'characters' should contain all possible discrete integer characters.")
      nl <- nlevels(kk)
    }
    if(identical(characters,TRUE))
    {
      if (!is.factor(x)) 
        tx <- factor(x)
      nl <- nlevels(tx)
      lvls <- levels(tx)
    } 
  } 
  if(length(lvls) ==1)
    stop("You need at least 2 different character states.")
  if(is.vector(x)&& identical(individualrates,TRUE))
    individualrates <- FALSE
  if(identical(individualrates, TRUE) && !identical(alpha ,FALSE))
    stop("You cannot set both individualrates = TRUE and alpha = TRUE/(a numeric value) at the same time.")    
  if(is.matrix(model))
  { 
    diag(model) <- 0
    np <- max(model[col(model)!= row(model)] )
    model[model==0] <- np+1    
  }
  nlmStart <- FALSE
  if(!identical(T, TRUE))
    stop("'T' should be set to be 'TRUE' for DiscML to work. type T <- TRUE in your interpereter to set 'T' to 'TRUE'.")
  if(!identical(F, FALSE))
    stop("'F' should be set to be 'FALSE' for DiscML to work. type F <- FALSE in your interpereter to set 'F' to 'FALSE'.")
  subcall <- match.call()
  if(identical( model , "GTR"))
    reversible <- TRUE
  if(identical(reversible, TRUE))
    rootprobability<-TRUE  
  command  <- FALSE
  if(!identical(T, TRUE))
    stop("The object T is not set to a logical value of TRUE. Before using this function, please change the value of T to TRUE by typing T <- TRUE")
  if(!identical(F, FALSE))
    stop("The object F is not set to a logical value of FALSE. Before using this function, please change the value of F to FALSE by typing F <- FALSE")
  templevels = NULL
  type ="discrete"
  method = "ML"
  .getSEs <- function(out)
  {
    h <- out$hessian
    if (any(diag(h) == 0)) {
      warning("The likelihood gradient seems flat in at least one dimension (gradient null):\ncannot compute the standard-errors of the transition rates.\n")
      se <- rep(NaN, nrow(h))
    }
    else {
      se <- suppressWarnings(sqrt(diag(solve(h))))
    }
    se
  } 
  if(is.matrix(model)|| toupper(model) =="GTR"|| toupper(model) == "ER" || toupper(model) == "BDBI"||toupper(model) =="ARD" || toupper(model) == "SYM"|| toupper(model) == "BDER" || toupper(model) == "BDSYM"||  toupper(model) =="BDARD"||  toupper(model) =="BDIER" || toupper(model) == "BDISYM"|| toupper(model) =="BDIARD")
  {}
  else{
    stop("You must input right 'model'. Please check your spelling")
  }  
  if(identical(alpha,FALSE) )
    talpha <- 1
  if (!inherits(phy, "phylo")) 
    stop("object \"phy\" is not of class \"phylo\"")
  if (is.null(phy$edge.length)) 
    stop("tree has no branch lengths")
  type <- match.arg(type, c("continuous", "discrete"))
  nb.tip <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (nb.node != nb.tip - 1) 
    stop("\"phy\" is not rooted AND fully dichotomous.")
  if (kappa != 1) 
    phy$edge.length <- phy$edge.length^kappa 
  if(nl==1)
    newnl<- 1
  if(nl > 1)
    newnl <- nl-1 
  rprobf <- function(...)
  {
    if(nl ==1)
      r1 <- paste("cos(", "t",")", 1, sep="")
    
    if(nl == 2)
      r1 <- paste("sin(","t",1,")", sep="")
    
    if(nl == 3)
    {
      r1 <- paste("sin(","t",1,")", sep="")
      for(v in 2:(nl-1))
      {
        r1 <- paste(r1, "*sin(","t", v, ")",sep ="")
      }
      r2 <- paste("sin(t1)*cos(t2)")
    }
    if(nl > 3 )
    { 
      r1 <- paste("sin(","t",1,")", sep="")
      
      for(v in 2:(nl-1))
      {
        r1 <- paste(r1, "*sin(","t", v, ")",sep ="")
      }
      for(vv in 2:(nl-2))
      {
        temp<- paste("sin(","t",1,")", sep="")
        for(v in 2:(nl-(vv)))
        {
          temp <- paste(temp, "*sin(","t", v, ")",sep ="")
        }
        assign( paste( "r", vv , sep="") ,  paste(temp, "*cos(","t", nl-(vv-1),")" ,sep=""))
      }
      assign(paste("r", nl-1, sep=""), paste("sin(t1)*cos(t2)"))
    }
    assign( paste("r", nl, sep="") , paste("cos(","t", 1, ")", sep=""))
    g<- list(...)
    if( length(g[[1]])>1)
    {
      number <- length(g[[1]])
      for(nnn in 1:number)
      {
        assign( paste("t", nnn, sep="") , g[[1]][nnn] )
      }
    }
    else
    {
      number <- length(g)
      for(nnn in 1:number)
      {
        assign( paste("t", nnn, sep="") , g[[nnn]] )
      }
    }
    for(i in 1:nl)
      assign( paste("x", i, sep=""),   eval(parse(text = get(paste("r", i, sep = "") ) )))
    temp <- paste( "c(", get("x1") )
    if(nl>1)
      for(i in 2:nl)
        temp <- paste( temp, ",", eval(parse(text = paste("x", i, sep="")  )  )  )
    temp <- paste(  temp, ")") 
    temp2<-vector()
    for(ii in 1:nl)
      temp2[ii] <-  paste( "(",eval(parse(text = (paste("r", ii, sep = "") ) ))  ,")^2", sep="" )
    if(length(g)>1 && is.character(g[[2]]))
    {
      return(  temp2)
    }
    else
      return( (eval(parse(text = temp)))^2)
  }
  
  if( identical(zerocorrection,TRUE) && !any(as.integer(lvls)==0))
    stop("Remember: If zerrocorection is set to TRUE, unobservable character state should be denoted by '0'.
         Please make sure that character state 0 is included in all possible character states.")
  if (method != "ML") 
    stop("only ML estimation is possible for discrete characters.")
  if (any(phy$edge.length <= 0)) 
    stop("some branches have length zero or negative")      
  if (!is.matrix(model)) {
    rate <- matrix(NA, nl, nl)
    switch(model, ER = np <- rate[] <- 1, ARD = {
      np <- nl * (nl - 1)
      rate[col(rate) != row(rate)] <- 1:np
    }, SYM = {
      np <- nl * (nl - 1)/2
      sel <- col(rate) < row(rate)
      rate[sel] <- 1:np
      rate <- t(rate)
      rate[sel] <- 1:np
    }, BDER = {
      np <- 1
      rate[row(rate)!=col(rate)] = np+1
      selg=row(rate)==(col(rate)-1)
      sell=row(rate)==(col(rate)+1)
      rate[selg]=1
      rate[sell]=1
      rate[is.na(rate)]=2
    }, BDARD = {          
      np <- 2*(nl-1)
      rate[row(rate)!=col(rate)] = np+1          
      for(i in 1:(nl-1))
      {
        rate[i, i+1]=i
        rate[i+1, i]= nl-1+i
      }
      rate[is.na(rate)]=np+1          
    }, BDSYM = {
      np <- (nl-1)
      rate[row(rate)!=col(rate)] = np+1
      for(i in 1:(nl-1))
      {
        rate[i, i+1]=i
        rate[i+1, i]=i
      }
      rate[is.na(rate)]=np+1          
    }, BDBI = {
      np <- 2
      rate[row(rate)!=col(rate)] = np+1
      for(i in 1:(nl-1))
      {
        rate[i, i+1]=1
        rate[i+1, i]=2
      }          
      rate[is.na(rate)]=np+1
    }, GTR = {
      np <- nl * (nl - 1)/2
      sel <- col(rate) < row(rate)
      rate[sel] <- 1:np
      rate <- t(rate)
      rate[sel] <- 1:np        
    },  BDISYM = {
      np <- 2
      rate[row(rate)!=col(rate)] = np+1
      rate[1, 2] =1
      rate[2, 1] =1
      for(i in 2:(nl-1))
      {
        rate[i, i+1]=2
        rate[i+1, i]=2
      }
      rate[is.na(rate)]=np+1          
    }, BDIER = {
      np <- 1
      rate[row(rate)!=col(rate)] = np+1
      selg=row(rate)==(col(rate)-1)
      sell=row(rate)==(col(rate)+1)
      rate[selg]=1
      rate[sell]=1
      rate[is.na(rate)]=2          
    }, BDIARD = {
      np <- 4
      rate = matrix(np+1,nl,nl)
      rate[1, 2] =1
      rate[2, 1] =2
      for(i in 2:(nl-1))
      {
        rate[i, i+1]=3
        rate[i+1, i]=4
      }
      rate[is.na(rate)]=np+1    
    })
  }
  else {
    rate = matrix(NA,nl,nl)
    if (ncol(model) != nrow(model)) 
      stop("the matrix given as 'model' is not square")
    if (ncol(model) != nl) 
      stop("the matrix 'model' must have as many rows as the number of categories in 'x'")
    rate[] <- suppressWarnings(as.integer(model))
    rate[is.na(rate)] <- 0       
    temp <- suppressWarnings(as.integer(model[col(model) != row(model)]))
    np <- max(temp[!is.na(temp)])
    diag(rate) <- np+1
    rate[rate==0] <- np+1 
  }
  if(!identical(rateparameters,TRUE)&&!(length(rateparameters)==np) && !reversible )
    stop("The number of rate parameters does not fit into the rate matrix, which is," ,np)
  if(identical(rateparameters,FALSE))
    stop("The argument 'p' cannot be set to FALSE") 
  if( identical(reversible,TRUE) ||  identical(model ,"GTR"))
  {
    if(is.matrix(model))
    {
      test <- identical(model, t(model))
      if(!test)
        stop("If 'reversible' is set to be TRUE, you must use symmetric metrix as a model. It is because when 'reversible' is set to be TRUE,
             DiscML converts Symmetric matrix into reversible matrix. (Please see description documents).")  
    }
    else if ( model == "ARD"|| model =="BDARD" || model =="BDIARD" )
      stop("If 'reversible' is set to be TRUE, you must use symmetric metrix as a model. It is because when 'reversible' is set to be TRUE,
           DiscML converts Symmetric matrix into reversible matrix. (Please see description documents).")
    }
  index.matrix <- rate
  tmp <- cbind(1:nl, 1:nl)
  index.matrix[tmp] <- NA
  rate[tmp] <-  np + 1
  TIPS <- 1:nb.tip
  if(identical(rootprobability, 1))
  {
    command <-TRUE
    rootprobability <- rep(1,nl)
  }
  if(is.numeric(mu))
  {
    if(length(mu) != nmuphy)
      stop("The number of elements in mu does not match the number of different mu's defined in the tree")
  }   
  phy <- reorder(phy, "postorder")
  
  Q <- matrix(0, nl, nl)
  e1 <- phy$edge[, 1]
  e2 <- phy$edge[, 2]
  EL <- phy$edge.length
  zliks <- matrix(0, NROW, nl)
  if(is.numeric(rootprobability) )
  {
    if(length(rootprobability)!= nl )
    {      
      stop("The number of elements of the root probability you have imputed does not match the number of possible character states, which is ", nl)      
    }
  }
  setIndividualratesTOFalse <-FALSE
  if(identical(individualrates,TRUE))
    setIndividualratesTOFalse <-TRUE
  if(is.numeric(irootprobability)&&identical(rootprobability,TRUE)&&!is.matrix(irootprobability))
  {
    if(length(irootprobability)!=nl)
      stop("The number of elements in 'irootprobability' does not match the number of possible character states, whiche is ", nl, "\n")
  }
  DiscML2 <- function (x, phy, CI = FALSE, model = "ER",rateparameters =TRUE, 
                       reversible =FALSE,  kappa = 1,
                       ip = 0.1, alpha=FALSE, ialpha= 0.5, ngamma = 8, zerocorrection=FALSE,
                       rootprobability=FALSE, irootprobability=FALSE, characters= TRUE,
                       plotmu= FALSE, simplify = FALSE ,individualrates = FALSE, plotloglik = FALSE) 
  {
    if(identical(setIndividualratesTOFalse,TRUE))
      NROW <-1
    count <- count+ 1
    likslist<- list()
    zlikslist<- list()
    likelivector <- numeric(0)
    obj <- list()
    xx<-x
    for( jj in 1:NROW )
    {
      if(is.matrix(xx))
      {
        x <-  xx[jj,]
      }
      if (length(x) != nb.tip) 
        stop("length of phenotypic and of phylogenetic data do not match.")
      if (!is.null(names(x))) {
        if (all(names(x) %in% phy$tip.label)) 
          x <- x[phy$tip.label]
        else if(!is.null(names(x)))
          stop("Some of names defined in 'x' does not match the tip labels.")
        else if(count == 1)
          warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
      }  
      if(is.matrix(xx))
      {
        if (is.numeric(characters))
        {
          xtemp<-factor(c(x,as.integer(lvls)))
          x <- as.integer(xtemp[1:length(x)]) 
        }
        if(identical(characters,TRUE))
        {
          if(!is.factor(x))
            kx <- factor(x, as.integer(lvls))
          else
            kx <- x
          x <- as.integer(kx)[1:length(x)]
        }
      }
      if(is.vector(xx))
      {
        if (is.numeric(characters))
        {
          if(!is.factor(xx))
            kx <- factor(xx, as.integer(lvls))
          
          x <- as.integer(kx)[1:length(x)]
        }
        if(identical(characters,TRUE))
        {
          if (!is.factor(x)) 
            x <- factor(x)
          x <- as.integer(x)
        } 
      }    
      liks <- matrix(0, nb.tip + nb.node, nl)
      liks[cbind(TIPS, x)] <- 1
      likslist[[jj]] <- liks
    }
    dev <- function(r, output.liks = FALSE, output.zliks= FALSE, output.mu = FALSE, output.Q=FALSE , output.div=FALSE, 
                    manualExceptAlpha=FALSE,  manualExceptMu = FALSE, manualExceptRprob =FALSE, manualExceptRates = FALSE,
                    estimateMu=rep(1, 2*nb.tip-2) , setAlphaToFalse =FALSE ) {
      if(setAlphaToFalse)
        alpha <- FALSE
      if(any(r<0))
        return(1e50)
      #print("#######################################################################START")
      #print("r")
      #print(r)
      talpha <- 1
      if (any(is.nan(r)) || any(is.infinite(r))) 
        return(1e+50)
      if(!identical(manualExceptMu, FALSE) )
      {      
        #print("#########################MU")
        length <- 0
        if(identical(rateparameters,TRUE))
        {
          p <- manualExceptMu[(length+1):(length+np)]
          length <- length + np
        }
        if(is.numeric(rateparameters))
        {
          p <- rateparameters
        }
        if(identical(rootprobability, TRUE)) 
        {
          rprob <- rprobf(manualExceptMu[(length+1):(length+(newnl))])
          length <- length + newnl
        }
        else if(identical(rootprobability, FALSE))
        {
          rprob <- rep(1/nl,nl)
        }
        else
        {
          rprob <- rootprobability/sum(rootprobability)
        }
        if(identical(alpha, TRUE) )        
          talpha <- manualExceptMu[length +1]
        if(is.numeric(alpha))
          talpha <- alpha
        for( i in 1:nmuphy)
        {
          estimateMu[muphy$edge.length==i] <- r[i] 
          if(!is.na(tempvalues[i]))
            estimateMu[muphy$edge.length ==i] <-as.numeric(tempvalues[i])
        }
      }
      else if(!identical(manualExceptAlpha,FALSE))
      {
        #print("#########################alpha")
        length <- 0
        if(identical(rateparameters,TRUE))
        {
          p <- manualExceptAlpha[(length+1):(length+np)]
          length <- length + np
        }
        if(is.numeric(rateparameters))
        {
          p <- rateparameters
        }
        if(identical(rootprobability, TRUE)) 
        {
          rprob <- rprobf(manualExceptAlpha[(length+1):(length+newnl)])
          length <- length + newnl
        }
        else if(identical(rootprobability, FALSE))
        {
          rprob <- rep(1/nl,nl)
        }
        else
        {
          rprob <- rootprobability/sum(rootprobability)
        }
        if(identical(mu,TRUE))
        {
          for( i in 1:nmuphy)
          {
            estimateMu[muphy$edge.length==i] <- manualExceptAlpha[length + i] 
            if(!is.na(tempvalues[i]))
              estimateMu[muphy$edge.length ==i] <-as.numeric(tempvalues[i])
          }
          
          length<- length + nmuphy
        }
        talpha <- r
        if(is.numeric(alpha))
          talpha <- alpha
      }
      else if(!identical(manualExceptRprob,FALSE))
      {
        #print("#########################rprob")
        
        length <- 0
        if(identical(rateparameters,TRUE))
        {
          p <- manualExceptRprob[(length+1):(length+np)]
          length <- length + np
        }
        if(is.numeric(rateparameters))
        {
          p <- rateparameters
        }
        if(identical(mu,TRUE))
        {
          for( i in 1:nmuphy)
          {
            estimateMu[muphy$edge.length==i] <- manualExceptRprob[length + i] 
            if(!is.na(tempvalues[i]))
              estimateMu[muphy$edge.length ==i] <- as.numeric(tempvalues[i])
          }
          length<- length + nmuphy
        }
        if(identical(alpha, TRUE) )
        {        
          talpha <- manualExceptRprob[length +1]  
          length <- length +1
        }
        else if(is.numeric(alpha))
        {
          talpha <- alpha
        }
        rprob <- rprobf(r)
      }
      else if(!identical(manualExceptRates,FALSE))
      {
        #print("#########################rate")
        length <- 0
        if(identical(rootprobability, TRUE)) 
        {
          rprob <- rprobf(manualExceptRates[(length+1):(length+newnl)])
          length <- length + newnl
        }
        else if(identical(rootprobability, FALSE))
        {
          rprob <- rep(1/nl,nl)
        }
        else
        {
          rprob <- rootprobability/sum(rootprobability)  
        }
        if(identical(mu,TRUE))
        {
          for( i in 1:nmuphy)
          {
            estimateMu[muphy$edge.length==i] <- manualExceptRates[length + i] 
            if(!is.na(tempvalues[i]))
              estimateMu[muphy$edge.length ==i] <- as.numeric(tempvalues[i])
          }
          length<- length + nmuphy
        } 
        else if(is.numeric(mu))
        {
          for( i in 1:nmuphy)
            estimateMu[muphy$edge.length ==i] <- mu[i]
        }
        if(identical(alpha, TRUE) )
        {        
          talpha <- manualExceptRates[length +1]  
          length <- length +1
        }
        else if(is.numeric(alpha))
        {
          talpha <- alpha
        }
        p <- r
      }
      else
      {
        #print("#########################free")
        length <- 0
        if(identical(rateparameters,TRUE))
        {
          p <- r[(length+1):(length+np)]
          length <- length + np
        }
        if(is.numeric(rateparameters))
        {
          p <- rateparameters
        }
        if(identical(rootprobability, TRUE))
        {
          rprob <- rprobf(r[(length+1):(length+newnl)])
          length <- length + newnl
        } 
        else if(identical(rootprobability, FALSE))
        {
          rprob <- rep(1/nl,nl)
        }
        else
        {
          rprob <- rootprobability/sum(rootprobability)
        }      
        if(identical(mu,TRUE))
        {
          for( i in 1:nmuphy)
          {
            estimateMu[muphy$edge.length==i] <- r[length + i] 
            if(!is.na(tempvalues[i]))
              estimateMu[muphy$edge.length ==i] <-as.numeric(tempvalues[i])
          }
          length<-length + nmuphy
        }
        if(identical(alpha, TRUE))
        {
          talpha <- r[length+1]
          length <- length +1
        }
        else if(is.numeric(alpha))
          talpha <- alpha   
        
        if(talpha > 50 && ialpha <50 && identical(alpha, TRUE) )
        {
          if(talpha > startalpha)
          {
            if(talphacount < 3)
              talphacount <<- talphacount + 1
            if(talphacount ==1)
              startalpha <<- talpha
          }
        } 
        if(talpha > startalpha)
        {
          multiplier <<- 1.1^(talpha - startalpha)
          talpha <- talpha*multiplier
        }
      }
      Q[] <- c(p, 0)[rate]
      diag(Q) <- rep(0, nl)   
      if(  reversible || identical(model , "GTR") )
        Q <- t(t(Q)*rprob)
      diag(Q) <- -rowSums(Q)
      div <- -sum(diag(Q)*rprob)        
      if(div==0)
        return(1e50)
      if(identical(mu, TRUE))
      {
        tQ <- Q/div      
      }
      if(identical(mu , FALSE))
        tQ <- Q
      if(!identical(alpha,FALSE))
      {
        #gamma rate 
        rateht<-matrix(NA, 1,ngamma)
        for(ii in 1:ngamma)
        {
          num1=(ii-0.5)/ngamma
          rateht[ii] <- qgamma(num1, shape=talpha, rate=talpha)
        } 
        rateht <- rateht/sum(rateht)*ngamma     
        if(any(is.na(rateht)))
          return(1e+50)      
        #gamma rate end
      }
      is(sum(rprob)!=1)
      rprob <- rprob/sum(rprob) 
      
      if(command)
        rprob<- rep(1,nl)
      decompo <- eigen(Q)
      lambda <- decompo$values
      GAMMA <- decompo$vectors
      tempInvGAMMA <- try(solve(GAMMA))
      if( class(tempInvGAMMA)=="try-error"   )
        return(1e50)
      else
        invGAMMA<- tempInvGAMMA  
      for(jj in 1:NROW)
      {
        liks <- likslist[[jj]]   
        if(identical(alpha,FALSE))
        {
          comp <- numeric(nb.tip + nb.node+1)
          for (i in seq(from = 1, by = 2, length.out = nb.node)) 
          {
            j <- i + 1L
            anc <- e1[i]
            des1 <- e2[i]
            des2 <- e2[j]
            v.l <- GAMMA %*% diag(exp(lambda * EL[i] * estimateMu[i]  )) %*% 
              invGAMMA %*% liks[des1, ]
            v.r <- GAMMA %*% diag(exp(lambda * EL[j] * estimateMu[j]  )) %*% 
              invGAMMA %*% liks[des2, ]
            v <- v.l * v.r
            comp[anc] <- sum(v)       
            liks[anc, ] <- v/comp[anc]      
          }
          
          likslist[[jj]] <- liks
          comp[nb.tip+nb.node +1] <- sum(liks[nb.tip+1, ]  * rprob)        
          if(!zerocorrection)
          {
            likeli =  sum(log(comp[-TIPS]))
          }       
          if(zerocorrection)
          {
            zliks = cbind(matrix(1, nrow(liks),1),matrix(0, nrow(liks), ncol(liks)-1))
            for (i in seq(from = 1, by = 2, length.out = nb.node)) 
            {
              j <- i + 1L
              anc <- e1[i]
              des1 <- e2[i]
              des2 <- e2[j]    
              v.l <- GAMMA %*% diag(exp(lambda * EL[i] * estimateMu[i]  )) %*% 
                invGAMMA %*% zliks[des1, ]
              v.r <- GAMMA %*% diag(exp(lambda * EL[j] * estimateMu[j] )) %*% 
                invGAMMA %*% zliks[des2, ]
              v <- v.l * v.r
              zliks[anc, ] <- v            
            }          
            zlikslist[[jj]] <- zliks
            zlikeli<- sum(zliks[nb.tip+1, ]*rprob)       
            if(zlikeli==1 || is.na(zlikeli))
            {
              return(1e+50)
            }
            likeli <- sum(log(comp[-TIPS]))-log(1-zlikeli)
          }
          likelivector[jj]<- likeli
        }
        
        if(!identical(alpha,FALSE))
        {
          likeli <-0
          for(w in 1:ngamma)
          {       
            for (i in seq(from = 1, by = 2, length.out = nb.node))           {
              j <- i + 1L
              anc <- e1[i]
              des1 <- e2[i]
              des2 <- e2[j]
              v.l <- GAMMA %*% diag(exp(lambda * estimateMu[i]* rateht[w] * EL[i])) %*% 
                invGAMMA %*% liks[des1, ]
              v.r <- GAMMA %*% diag(exp(lambda * estimateMu[j]* rateht[w] * EL[j])) %*% 
                invGAMMA %*% liks[des2, ]
              v <- v.l * v.r
              liks[anc, ] <- v
            }          
            likeli <- likeli + sum(liks[nb.tip+1, ]*rprob)/ngamma      
          }        
          likslist[[jj]] <- liks        
          if(!zerocorrection)
            likeli <- log(likeli)        
          if(zerocorrection)
          {        
            zlikeli <-0
            zliks = cbind(matrix(1, nrow(liks),1),matrix(0, nrow(liks), ncol(liks)-1))
            for(w in 1:ngamma)
            {
              for (i in seq(from = 1, by = 2, length.out = nb.node)) 
              {
                j <- i + 1L
                anc <- e1[i]
                des1 <- e2[i]
                des2 <- e2[j]    
                v.l <- GAMMA %*% diag(exp(lambda * EL[i] *rateht[w]*  estimateMu[i]  )) %*% 
                  invGAMMA %*% zliks[des1, ]
                v.r <- GAMMA %*% diag(exp(lambda * EL[j] *rateht[w]*  estimateMu[j] )) %*% 
                  invGAMMA %*% zliks[des2, ]
                v <- v.l * v.r
                zliks[anc, ] <- v            
              }          
              zlikeli<- zlikeli +sum(zliks[nb.tip+1, ]*rprob)/ngamma
            }
            zlikslist[[jj]] <- zliks
            likeli <- ( log(likeli) -log(1-zlikeli)  )          
          }        
        }      
        likelivector[jj] <-  likeli      
      } 
      #print("liks")
      #print(liks)
      #print("p")
      #print(p)
      #print("Q")
      #print(Q)
      if(!identical(alpha,FALSE))
      {
        #print("alpha")
        #print(talpha)
      }
      #print("mu")
      #print(estimateMu)
      #print("tempchars")
      #print(tempchars)
      #print("tempvalues")
      #print(tempvalues)
      #print("rprob")
      #print(rprob)
      #print("rate")
      #print(rate)
      #print("edge")
      #print(phy$edge)
      #print("muphy edge length")
      #print(muphy$edge.length)
      #print("estimateMu")
      #print(estimateMu)
      #print("likelivector")
      #print(likelivector)
      output = -2*sum(likelivector)    
      if(output.liks)
        return(likslist)
      if(output.zliks)
        return(zlikslist)   
      if(output.Q)
        return(tQ)  
      if(output.div)
        return(div)  
      if(output.mu)
        return(estimateMu)  
      if (is.na(output)||is.infinite(output)) 
        10e50
      else 
      {
        #print("out")
        #print(-1/2*output)
        Re(output)
      }
    } 
    length <- 0
    if(identical(rateparameters,TRUE))
    {
      rateStart <- length+1
      rateEnd <- length+np
      rateIndicator <- 1
      length <- length +np
    }
    if(!identical(rateparameters,TRUE))
    {
      rateIndicator <- 0
    }
    
    if(identical(rootprobability, TRUE))
    {
      rprobStart <- length+1
      rprobEnd <- length+newnl  
      rprobIndicator <-1
      length <- length + newnl
    }
    if(!identical(rootprobability, TRUE))
    {
      rprobIndicator <-0
    }
    if(identical(mu, TRUE))
    {
      muStart <- length+1
      muEnd <- length + nmuphy
      muIndicator <- 1
      length <- length + nmuphy
    }
    if(!identical(mu, TRUE))
    {
      muIndicator <- 0
    }
    if(identical(alpha, TRUE))
    {
      alphaStart <- length+1
      alphaEnd <- length+1
      alphaIndicator <- 1
      length <- length + 1
    }
    if(!identical(alpha, TRUE))
    {
      alphaIndicator <- 0
    }
    if(is.numeric(irootprobability))
    {
      irprob <- irootprobability/sum(irootprobability)
    }
    else
      irprob <- rep(1/nl, nl)
    
    if(nl == 1)
      irprob<- nlminb( 0.1 , function(r) {sum(abs(rprobf(r)-irprob)) }, lower= -2*pi, upper = 4*pi)$par
    if(nl > 1)
      irprob<- nlminb( rep(0.1, nl-1) , function(r) {sum(abs(rprobf(r)-irprob)) } , lower= rep(-2*pi, nl-1), upper = c(rep(2*pi,nl-2), 4*pi))$par
    if(nl==1)
    {
      lowerRprob <- -2*pi
      upperRprob <- 4*pi
    }
    if(nl >1)
    {
      lowerRprob <- rep(-2*pi, nl-1)
      upperRprob <- c(rep(2*pi,nl-2), 4*pi)
    }
    obj$rootprobability <- rootprobability
    out <- list()
    temp <- list()
    pars  <- numeric()
    
    lowerPoint <- c( rep(rep(0,np), rateIndicator), rep(lowerRprob, rprobIndicator),  rep(rep(0.0001,nmuphy), muIndicator) , rep(0.0001, alphaIndicator) )
    upperPoint <- c( rep(rep(1e50,np), rateIndicator), rep(upperRprob, rprobIndicator),  rep(rep(1e50, nmuphy), muIndicator) , rep(1e10, alphaIndicator) )
    #print("#########################Optimize#######################################")
    initialPoint <- c( rep(rep(ip,np), rateIndicator), rep(irprob, rprobIndicator), rep(rep(imu,nmuphy), muIndicator) , rep(ialpha, alphaIndicator) ) 
    if(!identical(initialPoint, numeric(0)))
      out <- nlminb(initialPoint, function(r)dev(r), lower= lowerPoint , upper = upperPoint  )
    else
      out <- nlminb(1, function(r)dev(r), lower= 1 , upper = 1  )
    if(identical(alpha,TRUE) && out$par[length(out$par)] > 3000)
    {
      lowerPoint <- c( rep(rep(0,np), rateIndicator), rep(lowerRprob, rprobIndicator),  rep(rep(0.0001,nmuphy), muIndicator) , rep(0.0001, 0) )
      upperPoint <- c( rep(rep(1e50,np), rateIndicator), rep(upperRprob, rprobIndicator),  rep(rep(1e50, nmuphy), muIndicator) , rep(1e10, 0) )
      #print("#########################Optimize#######################################")
      initialPoint <- c( rep(rep(ip,np), rateIndicator), rep(irprob, rprobIndicator), rep(rep(imu,nmuphy), muIndicator) , rep(ialpha, 0) )
      if(!identical(initialPoint, numeric(0)))
        out2 <- nlminb(initialPoint, function(r)dev(r, setAlphaToFalse =TRUE), lower= lowerPoint , upper = upperPoint  )
      else
        out2 <- nlminb(1, function(r)dev(r, setAlphaToFalse =TRUE), lower= 1 , upper = 1  )
      if(out2$objective < out$objective)
      {
        out <<- out
        out2 <<- out2
        out <- out2
        out$par[length(out$par)+1] <- 1e50
      }
    }
    obj$lvls <- lvls
    if(identical(reversible,TRUE))
    {
      rootprobability <- TRUE
    }
    obj$srprob <- matrix( rep(0, 2*(newnl)), ncol =2)
    obj$isrprob<- rootprobability
    obj$isrates<- rateparameters
    obj$ismu <- mu
    obj$isalpha <- alpha
    obj$exrprob <- vector()
    checkQ <- dev(out$par , output.Q=T)
    checkDiv <- dev(out$par ,  output.div=T)   
    
    
    if(is.numeric(rootprobability))
    {
      obj$rprob<- rootprobability/sum(rootprobability)
    }
    else if(identical(rootprobability, FALSE))
    {
      obj$rprob<- rep(1/nl, nl)
    }
    else if(identical(rootprobability, TRUE) )
    {
      obj$rprob <-  rprobf(out$par[rprobStart:rprobEnd] )
      obj$srprob[,1] <- out$par[rprobStart:rprobEnd]
      obj$exrprob <- rprobf(out$par[rprobStart:rprobEnd], "char" )
    }
    
    if(identical( rateparameters,TRUE))
    {
      obj$rates  <- matrix(c(out$par[rateStart:rateEnd], rep(0, np)), ncol = 2)  
    }
    else if(is.numeric(rateparameters) && is.vector(rateparameters))
    {
      obj$rates <- matrix(c(rateparameters, rep(0, np)), ncol = 2)
    }
    if(identical(mu,TRUE))
    {
      out$par[muStart:muEnd] <- out$par[muStart:muEnd] *checkDiv
      obj$rates[,1] <- obj$rates[,1]/checkDiv
      if(identical( rateparameters,TRUE))
      {
        out$par[rateStart:rateEnd]  <-       out$par[rateStart:rateEnd]/checkDiv
      }
      for( i in 1:nmuphy)
      {
        if(!is.na(tempvalues[i]))
          out$par[muStart:muEnd][i] <- tempvalues[i]
      }
      obj$mu <- out$par[muStart:muEnd]
    }
    else if(identical(mu,FALSE))
    {
      obj$mu <- rep(1, nmuphy )
    }
    alphaInf <- FALSE
    if(identical(alpha,TRUE))
    {
      if(out$par[alphaStart:alphaEnd] == 1e50)
        alphaInf <- TRUE
      if(out$par[alphaStart:alphaEnd] > 50 && out$par[alphaStart:alphaEnd] > startalpha)
        out$par[alphaStart:alphaEnd] <-  out$par[alphaStart:alphaEnd] * multiplier
      obj$alpha <- matrix( c( out$par[alphaStart:alphaEnd] , rep(0, 1) ), ncol = 2 )
    }
    else if(is.numeric(alpha) )
    {
      obj$alpha <- matrix( c(alpha, rep(0, 1)) , ncol = 2)
    } 
    obj$loglik <- -out$objective/2
    
    oldwarn <- options("warn")
    options(warn = -1)
    out.nlm <- list()
    nlmStart <- TRUE
    if(identical(rateparameters,TRUE))
    {
      #print("####rate nlm###################################")
      out.nlm <- try(nlm(function(r) dev(r, manualExceptRates= out$par[ -(rateStart:rateEnd) ] ), p = obj$rates[,1], 
                         iterlim = 1, stepmax = 0, hessian = TRUE))
      options(oldwarn)
      obj$rates[,2] <- if (class(out.nlm) == "try-error") {
        warning("model fit suspicious: gradients apparently non-finite")
        c(rep(NaN, np))
      }
      else if( class(try(suppressWarnings(.getSEs(out.nlm)))) =="try-error")
        c(rep(NaN, np))
      else
        suppressWarnings(.getSEs(out.nlm ))  
    }
    else
      obj$rates[,2] <- rep(0, np)
    
    
    
    if(identical(mu, TRUE))
      obj$rates[,2] <- rep(0,np)
    if(  is.matrix(model) &&!reversible )
    {
      QQ <- dev(out$par , output.Q=TRUE)
      rate2<- model
      diag(rate2)<- rep(NA, nl)
      
      rownamesRates <- c(1:np)
      index.matrix <- rate2
      temp <- as.numeric(QQ[row(QQ)!=col(QQ)])
    } 
    obj$rprob <- obj$rprob / sum(obj$rprob)
    if(identical(rootprobability, TRUE))
    {
      #print("####rprob nlm###################################")
      rprob.nlm <- try(nlm(function(r) dev(r, manualExceptRprob= out$par[ -(rprobStart:rprobEnd) ]), p = out$par[ rprobStart:rprobEnd ], 
                           iterlim = 1, stepmax = 0, hessian = TRUE))
      options(oldwarn)
      obj$srprob[,2]<- if (class(rprob.nlm ) == "try-error") {
        warning("model fit suspicious: gradients apparently non-finite")
        c(rep(NaN, newnl))
      }
      else if( class(try(suppressWarnings(.getSEs(rprob.nlm )))) =="try-error")
        c(rep(NaN, newnl))
      else
        suppressWarnings(.getSEs(rprob.nlm ))  
    }
    else
      obj$srprob[,2] <- rep(0, newnl)
    
    if(identical(rootprobability, TRUE))
    {      
      obj$srprob <- cbind( paste("t", 1:(newnl), sep="") , as.data.frame(obj$srprob))
      colnames(obj$srprob) <- c("angle", "estimate(radian)", "std-err")      
      obj$rprob <- data.frame( as.integer(lvls), obj$exrprob, obj$rprob)    
      mm <- function(x){ 
        for(ww in 1:length(obj$srprob[,1]))
        {
          assign(  levels(unlist(obj$srprob[,1]))[ww] , x[ww] )
        }
        return(eval(parse(text= levels(obj$rprob[,2])[zz] )))
      } 
      LL <- numeric()
      UU <- numeric()
      for(zz in 1:length(obj$rprob[,1]))
      {
        UU[zz] <-max(c( mm(    obj$srprob[,2] + obj$srprob[,3]  ) ,mm(obj$srprob[,2] - obj$srprob[,3] )))
        LL[zz] <-min(c( mm(      obj$srprob[,2] + obj$srprob[,3] ) ,mm(obj$srprob[,2] - obj$srprob[,3] )))
      }
      diffLL <- abs(LL-mm(obj$srprob[,2]))
      diffUU <- abs(UU-mm(obj$srprob[,2]))
      
      ffff <- numeric()
      fino <-matrix(c(diffLL, diffUU), ncol = 2)
      for(hh in 1:nrow(fino))
        ffff[hh]<- max(fino[hh,]) 
      obj$rprob <- cbind(obj$rprob, ffff)
      colnames(obj$rprob) <- c("characters" ,"estimate", "std-err")
      obj$rprob <- obj$rprob[ ,-2] 
    }
    else
    {
      obj$rprob <- cbind( as.integer(lvls), obj$rprob)
      colnames(obj$rprob) <- c("characters", "value")
      rownames(obj$rprob) <- rep("", nl)
    }
    if(identical(alpha, TRUE))
    {
      #print("####alpha nlm###################################")
      alpha.nlm <- try(nlm(function(r) dev(r, manualExceptAlpha= out$par[ -(alphaStart:alphaEnd) ]), p = obj$alpha[,1], 
                           iterlim = 1, stepmax = 0, hessian = TRUE))
      options(oldwarn)
      obj$alpha[,2] <- if (class(alpha.nlm) == "try-error") {
        warning("model fit suspicious: gradients apparently non-finite")
        c(rep(NaN, 1))
      }
      else if( class(try(suppressWarnings(.getSEs(alpha.nlm)))) =="try-error")
        c(rep(NaN, 1))
      else
        suppressWarnings(.getSEs(alpha.nlm))  
    }
    else if(is.numeric(alpha))
      obj$alpha[,2] <- rep(0, 1)
    
    if(alphaInf)
    {
      obj$alpha[,1] <-Inf
      obj$alpha[,2] <- NaN
    }
    if(identical(mu, TRUE))
    {
      
      #print("####mu nlm###################################")
      mu.nlm <- try(nlm(function(r) dev(r, manualExceptMu= out$par[ -(muStart:muEnd) ]), p = obj$mu, 
                        iterlim = 1, stepmax = 0, hessian = TRUE))
      options(oldwarn)
      obj$smu <- if (class(mu.nlm) == "try-error") {
        warning("model fit suspicious: gradients apparently non-finite")
        c(rep(NaN, nmuphy))
      }
      else if( class(try(suppressWarnings(.getSEs(mu.nlm)))) =="try-error")
        c(rep(NaN, nmuphy))
      else
        suppressWarnings(.getSEs(mu.nlm))  
    }
    else
      obj$smu <- rep(0, nmuphy)    
    if(reversible)
    {
      if( !is.matrix(model)&& toupper(substring(model, 1,2)) =="BD")
      {
        QQ <- dev(out$par, output.Q =TRUE)
        np2 <- 2*(nl-1)
        rate2 <- rate
        rate2[row(rate2)!=col(rate2)] <- np2+1
        for(i in 1:(nl-1))
        {
          rate2[i, i+1]=i
          rate2[i+1, i]= nl-1+i
        }
        rate2[is.na(rate2)] <- np2+1
        diag(rate2) <- rep(np2+1, nl)
        index.matrix <- rate2
        temp<- c( QQ[ col(QQ) == row(QQ)+1 ] , QQ[ col(QQ) == row(QQ)-1 ]    )
        stemp<-numeric()
        for(mm in 1:np2)
        {
          tt<-c(rate[ col(rate) == row(rate)+1 ] , rate[ col(rate) == row(rate)-1 ]  )[mm]
          if( tt!=0)
            stemp[mm] <-  obj$rates[,2][tt]
          else
            stemp[mm] <- 0
        }
        obj$rates <-matrix( c(temp,stemp) , ncol=2)
      }
      else
      {
        rownamesRates <- 1:(nl*(nl-1))
        QQ <- dev(out$par, output.Q =TRUE)
        np2 <- nl * (nl - 1)
        rate2= matrix(np2+1,nl,nl)
        rate2[col(rate2) != row(rate2)] <- 1:np2
        index.matrix <- rate2
        temp<- as.numeric(QQ[row(QQ)!=col(QQ)])
        stemp<-numeric()
        for(mm in 1:np2)
        {
          tt<-suppressWarnings(as.integer(rate[row(rate)!= col(rate)][mm]))
          if( tt!=0)
            stemp[mm] <-  obj$rates[,2][tt]
          else
            stemp[mm] <- 0
        }
        obj$rates <- matrix(c(temp, rep(0, length(temp))), ncol = 2)
        obj$rates[,2] <- stemp
      }
    }
    obj$call <- subcall
    class(obj) <- "DiscML"
    if(reversible)
    {
      index.matrix[ is.na(index.matrix)] <- np2+1
      index.matrix[index.matrix == np2+1] <- rep(NA, length(index.matrix[index.matrix == np2+1]))
    }
    else
    {
      index.matrix[ is.na(index.matrix)] <- np+1
      index.matrix[index.matrix == np+1] <- rep(NA, length(index.matrix[index.matrix == np+1]))
      if(is.matrix(model))
      {
        index.matrix[ is.na(index.matrix)] <- np+1
        index.matrix[row(index.matrix) == col(index.matrix)] <- rep(NA, nl)
      }
    }
    obj$index.matrix <- index.matrix
    tempQQ <- matrix(NA, nl ,nl)
    
    if(!is.matrix(model))
    {
      rate3<- rate
      rate3[rate3==0] <- rep(np+1 , length(rate3[rate3==0]) )
      
      tempQQ[] <-  c(obj$rates[ ,2], 0)[rate3]
      if( !(toupper(substring(model, 1,2)) == "BD") && (if(!is.matrix(model)){ toupper(model) == "ARD" }  || reversible ) )
        obj$rates[,2] <- c(tempQQ[col(tempQQ) != row(tempQQ)])
      if( toupper(substring(model, 1, 2) ) =="BD" && reversible || toupper(model) =="BDARD")
      {
        obj$rates[,2]<- c(tempQQ[col(tempQQ) ==1 + row(tempQQ)], tempQQ[col(tempQQ) == -1 + row(tempQQ)]  )
      }
    }
    if (CI)
    {
      obj$tobj = matrix(0,NROW, nl)
      obj$lik.anc <- dev(out$par,  output.liks=TRUE)
      for(i in 1:NROW)
      {
        colnames(obj$lik.anc[[i]]) <- lvls
        obj$tobj[i,] <- obj$lik.anc[[i]][nb.tip+1,]
        obj$lik.anc[[i]] <- obj$lik.anc[[i]]
        names(obj$lik.anc)[i] <- paste("Data ", i, " :")
      }
      cc<- numeric()
      colnames(obj$tobj) <- lvls
      for( i in 1:NROW)
        cc[i] <- paste("Data", i, ": ") 
      rownames(obj$tobj) <- cc
      obj$tobj <- obj$tobj/rowSums(obj$tobj)
    }
    npar <- length(obj$rates[,1])
    if(!is.matrix(model))
      obj$rates <- data.frame(1:npar, obj$rates[,1], obj$rates[,2][1:npar])
    if(is.matrix(model))
      obj$rates <- data.frame(rownamesRates, obj$rates[,1], obj$rates[,2][1:npar])
    colnames(obj$rates) <- c("rate index", "estimate", "std-err")
    if(!identical(alpha,FALSE) )
    {
      colnames(obj$alpha) <- c("estimate", "std-err")
      rownames(obj$alpha) <- "alpha: "
    } 
    obj$mu <- data.frame(tempchars, obj$mu, obj$smu )
    colnames(obj$mu) <- c("mu(ID)", "estimate", "std-err")
    Order <- order(obj$mu[,1])
    trates<-list()
    tmu <- obj$mu
    trprob <- list()
    
    for(gg in 1:nmuphy)
    {
      tmu[gg,] <- obj$mu[Order[gg], ]
      trprob[[gg]] <- obj$rprob[[Order[gg]]]
    }
    obj$mu <- tmu
    obj$trprob <- trprob
    obj$Order <- Order 
    obj$nmuphy <- nmuphy
    obj$nl <- nl
    obj$NROW <- NROW
    obj$individualrates <- individualrates
    obj$lik.anc <- obj$tobj
    obj$tobj <- NULL
    obj$tempchars <- tempchars
    obj$irootprobability <- irootprobability    
    end <- proc.time() - start
    obj$time <- end[3]
    obj$temprates <- obj$rates
    if(identical(mu, TRUE))
    {
      obj$rates <- obj$rates[, colnames(obj$rates)!= "std-err"]
      colnames(obj$rates)<- c("rate.index", "value")
    }
    if(identical(rootprobability, TRUE))
      colnames(obj$rprob) <- c("characters", "estimate","std-err")
    return(obj)
  }
  if(identical(individualrates , FALSE))
  {
    #sink("DiscML Debug.txt")
    return( DiscML2(x =x, phy =phy, CI = CI, model = model, rateparameters =rateparameters,
                    reversible =reversible,  kappa = kappa,
                    ip = ip, alpha=alpha, ialpha= ialpha, ngamma = ngamma, zerocorrection=zerocorrection,
                    rootprobability=rootprobability, irootprobability=irootprobability, characters= characters,
                    plotmu =plotmu, simplify = simplify ,individualrates = individualrates, plotloglik = plotloglik) )
  }
  if(identical(individualrates ,TRUE))
  {
    #sink("DiscML Debug.txt")
    column <- numeric(2*NROW)
    rcolumn <- numeric(NROW)
    for(i in 1:NROW)
    {
      column[  2*i-1] <- paste( "estimate", "#", i, sep="" ) 
      column[  2*i] <- paste( "std-err", "#", i, sep="" )   
      rcolumn[i] <- paste("estimate", "#", i, sep = "")
    }
    final <- list()
    ultimate <- list()
    ultimate$loglik <- matrix(0, 1,NROW)
    if(identical(rootprobability,TRUE))
      ultimate$rprob <- matrix(0, NROW, nl*2)
    else
      ultimate$rprob <- matrix(0, nl, NROW) 
    ultimate$srprob <- matrix(0, nl-1, 2*NROW)
    ultimate$individualrates <- individualrates
    ultimate$mu <- matrix(0, nmuphy, NROW*2)
    colnames(ultimate$mu) <- column
    if(!identical(alpha,FALSE))
    {
      ultimate$alpha <- matrix(0, 1, NROW*2)
    } 
    xx<-x
    final[[1]] <- DiscML2(x =xx[1,], phy = phy, CI = CI, model = model, rateparameters = rateparameters,                          
                          reversible =reversible,  kappa = kappa,
                          ip = ip, alpha=alpha, ialpha= ialpha, ngamma = ngamma, zerocorrection=zerocorrection,
                          rootprobability=rootprobability, irootprobability=irootprobability, characters= if(identical(characters,TRUE)){as.integer(lvls)}else{ characters},
                          plotmu = plotmu, simplify = simplify ,individualrates = FALSE,plotloglik = plotloglik )
    npar <- length(final[[1]]$temprates[,1])
    ultimate$lvls <- final[[1]]$lvls
    ultimate$index.matrix <- final[[1]]$index.matrix
    ultimate$lik.anc <- final[[1]]$lik.anc
    ultimate$ismu <- final[[1]]$ismu
    ultimate$isalpha <- final[[1]]$isalpha
    if(identical(mu, TRUE))
      ultimate$tempchars <- final[[1]]$tempchars
    
    if(identical(mu,TRUE))
    {
      ultimate$rates <- matrix(0, npar, NROW)
      ultimate$rates[,1] <- as.matrix(final[[1]]$temprates[, 2])
    }
    else
    {
      ultimate$rates <- matrix(0, npar, 2*NROW)
      ultimate$rates[,1:2] <- as.matrix(final[[1]]$temprates[, 2:3])
    }
    ultimate$loglik[1,1] <- final[[1]]$loglik
    if(identical(rootprobability, TRUE))
    {
      for(rr in 1:nl)
      {
        ultimate$rprob[ 1,2*rr-1 ]<- (final[[1]]$rprob[ rr,2])
        ultimate$rprob[ 1,2*rr ]<- (final[[1]]$rprob[ rr,3])      
      }
      ultimate$srprob[, 1:2]<- as.matrix(final[[1]]$srprob[ ,2:3])
      ultimate$exrprob <- final[[1]]$exrprob
    }
    else
      ultimate$rprob[, 1]<- as.numeric(final[[1]]$rprob[ ,2])
    ultimate$mu[, (2*1-1):(2*i)]<- as.matrix(final[[1]]$mu[, 2:3])
    if(!identical(alpha,FALSE))
      ultimate$alpha[, (2*1-1):(2*1)] <- final[[1]]$alpha  
    for(i in 2:NROW)
    {
      final[[i]] <- DiscML2(x =xx[i,], phy = phy, CI = CI, model = model, rateparameters = rateparameters,                            
                            reversible =reversible,  kappa = kappa,
                            ip = ip, alpha=alpha, ialpha= ialpha, ngamma = ngamma, zerocorrection=zerocorrection,
                            rootprobability=rootprobability, irootprobability=irootprobability , characters= if(identical(characters,TRUE)){as.integer(lvls)}else{ characters},
                            plotmu = plotmu ,simplify = simplify ,individualrates = FALSE,plotloglik = plotloglik )
      ultimate$loglik[1,i] <- final[[i]]$loglik
      if(identical(mu,TRUE))
        ultimate$rates[  , i] <- as.matrix(final[[i]]$temprates[ ,2])    
      else
        ultimate$rates[  , (2*i-1):(2*i)] <- as.matrix(final[[i]]$temprates[ ,2:3])
      if(identical(rootprobability, TRUE))
      {
        for(rr in 1:nl)
        {
          ultimate$rprob[ i,2*rr-1 ]<- (final[[i]]$rprob[ rr,2])
          ultimate$rprob[ i, 2*rr ]<- (final[[i]]$rprob[ rr,3])      
        }
        ultimate$srprob[, (2*i-1):(2*i)]<- as.matrix(final[[i]]$srprob[ ,2:3])
      }
      else
        ultimate$rprob[, i]<- as.numeric(final[[i]]$rprob[ ,2])
      ultimate$mu[, (2*i-1):(2*i)]<- as.matrix(final[[i]]$mu[, 2:3])
      if(!identical(alpha,FALSE))
        ultimate$alpha[, (2*i-1):(2*i)] <- final[[i]]$alpha  
    }
    npar <- length(ultimate$rates[,1])
    ultimate$rates<- (cbind(1:npar,   ultimate$rates ))
    if(identical(mu,TRUE))
      colnames(ultimate$rates) <- c("rate index", rcolumn) 
    else
      colnames(ultimate$rates) <- c("rate index", column)
    rownames(ultimate$rates) <- rep("", npar) 
    if(identical(rootprobability, TRUE))
    {
      ultimate$srprob <- cbind( paste("t", 1:(nl-1), sep="") , as.data.frame(ultimate$srprob))
      colnames(ultimate$srprob) <- c("angles(radian)", column)
    }
    else
    {
      ultimate$rprob <- cbind( as.integer(lvls) , ultimate$rprob)
      colnames(ultimate$rprob) <- c("character", rcolumn)
      rownames(ultimate$rprob) <- rep("", nl)
    }
    #print(2)
    logliknames<-numeric()
    for (i in 1:NROW)
      logliknames[i] <- paste( "loglik", "#", i, sep="" ) 
    colnames(ultimate$loglik) <- logliknames
    rownames(ultimate$loglik) <- ""
    ultimate$mu <- data.frame( tempchars, as.data.frame(ultimate$mu)  ) 
    colnames(ultimate$mu) <- c("mu(ID)", column)
    if(!identical(alpha,FALSE) )
    {
      colnames(ultimate$alpha) <- column
      rownames(ultimate$alpha) <- "alpha: "
    }  
    class(ultimate) <- "DiscML"
    ultimate$call <- subcall
    tempn <- 1:NROW
    
    ultimate$nmuphy <- nmuphy
    ultimate$tempchars <- tempchars
    Order <- order(ultimate$mu[,1])
    trates<-list()
    tmu <- ultimate$mu
    trprob <- list()
    for(hh in 1:nmuphy)
    {
      tmu[hh,] <- ultimate$mu[Order[hh], ]
      trprob[[hh]] <- ultimate$rprob[[Order[hh]]]
    }
    ultimate$NROW <- NROW
    ultimate$mu <- tmu
    ultimate$trprob <- trprob
    ultimate$Order <- Order 
    ultimate$isrprob<- rootprobability
    ultimate$irootprobability <- irootprobability
    
    tplotmu <- 0
    if(!identical(plotmu,FALSE))
    {
      if(identical(plotmu, TRUE) && identical(mu,TRUE))
      {
        tplotmu <- sort(tempchars)
      }
      else if(is.character(plotmu))
      {
        tplotmu<- gsub("\n","", plotmu)
        tplotmu<- gsub(" " , "" ,tplotmu)
        tplotmu<- sort(unlist(strsplit(tplotmu, ',')))
      }
    }
    if( identical(plotloglik , TRUE) )
    {
      if( !identical(plotmu , FALSE) )
        par(mfcol=c(length(tplotmu) + 1 ,1))
      
      if( identical(plotmu , FALSE) )
        par(mfcol=c( 1 ,1))
    }
    if( identical(plotloglik , FALSE) )
    {
      if( !identical(plotmu , FALSE) )
        par(mfcol=c(length(tplotmu)  ,1))
      if( identical(plotmu , FALSE) )
      {}
    }
    if(!identical(plotmu,FALSE))
    {
      temp <- list()
      tv<- numeric()
      for( i in 1:length(tempchars))
        if(!(tempchars[i] %in% tplotmu))
          tv[i] <- i
      tv <- tv[!is.na(tv)]
      if(identical(tv, numeric(0)))
        mat<- t(as.matrix(ultimate$mu[,-1]))
      else
        mat<-  t(as.matrix(ultimate$mu[,-1]))[,-(tv)]  
      mat <- mat[ (row(mat)%%2!=0)[,1] ,  ]
      if(!is.vector(mat))
        mat <- matrix(as.numeric(mat), ncol= ncol(mat), nrow = nrow(mat))
      else
        mat <- as.numeric(mat)
      Data <- data.frame(tempn ,mat )
      
      colnames(Data) <- c("Sites", sort(tplotmu))    
      scinot <- function(x, digits = 1) {
        if (length(x) > 1) {
          return(append(scinot(x[1]), scinot(x[-1])))
        }
        if (!x) return(0)
        exponent <- floor(log10(x))
        base <- round(x / 10^exponent, digits)
        as.expression(substitute(base %*% 10^exponent, 
                                 list(base = base, exponent = exponent)))
      }    
      for(aa in 1:length(tplotmu))
      {
        par(las=1)
        if(all(Data[,1+aa] >999)||all(Data[,1+aa] < 1e-5) )
        {
          par(mar = c(6.5, 6.5, 2, 0.5) + 0.1, mgp = c(5, 1, 0))
          if(1==length(length(tplotmu)))
            plot( Data[ ,1], Data[ ,1+aa], xlab = "Sites", ylab = "Individual rates" ,yaxt='n',  xaxt='n' )    
          else
            plot( Data[ ,1], Data[ ,1+aa], xlab = "Sites", ylab = "Individual rates", main= sort(tplotmu)[aa] ,yaxt='n',  xaxt='n' )
          axis(2, at = axTicks(2), labels = scinot(axTicks(2), 3))
        }
        else
        {
          par(mar = c(6.5, 6.5, 2, 0.5) , mgp = c(5, 1, 0))
          if(1==length(length(tplotmu)))
            plot( Data[ ,1], Data[ ,1+aa], xlab = "Sites", ylab = "Individual rates" ,yaxt='n',  xaxt='n' )    
          else
            plot( Data[ ,1], Data[ ,1+aa], xlab = "Sites", ylab = "Individual rates", main= sort(tplotmu)[aa] , yaxt='n', xaxt='n' )
          axis(2, at = axTicks(2), labels = axTicks(2))
        }
        lablist<-as.vector(tempn)
        axis(1, at=tempn, labels = FALSE)
        if(NROW>100)
        {
          text(tempn, par("usr")[3]-1, labels = lablist, srt = 45, pos = 1, xpd = TRUE)
        }
        else
        {
          text(tempn, par("usr")[3]-1 , labels = lablist, srt = 0, pos = 1, xpd = TRUE)
        }
      }
    }
    options(scipen= 5)
    if(identical(plotloglik, TRUE))
    {
      Data <- data.frame(x= tempn, y= t(ultimate$loglik))
      
      plot(tempn , ultimate$loglik, xlab = "Sites", ylab= "Log Likelihoods", xaxt='n')
      axis(1, at=tempn) 
    }
    tlist <- list()
    dlist <- list()
    if(identical(mu, TRUE))
    {
      tplotmu <- sort(tempchars)
    }
    else
    {
      tplotmu <-NULL
    }
    for(bb in 1:NROW)
    {
      if(identical(mu,TRUE))
      {
        ttt <- numeric()
        trates <- ultimate$rates[ ,bb+1, drop=FALSE]
        frate <- numeric()
        for(cc in 1:nrow(trates))
          ttt<- c(ttt, trates[cc,])
        frates <- ttt
        nrates <- vector()
        for(cc in 1:(nrow(trates)))
          nrates<- c(nrates, c(paste(cc , sep="")))
      }
      else
      {
        ttt <- numeric()
        trates <- ultimate$rates[ ,(2*bb):(2*bb+1), drop=FALSE]
        frate <- numeric()
        for(cc in 1:nrow(trates))
          ttt<- c(ttt, trates[cc,])
        frates <- ttt
        nrates <- vector()
        for(cc in 1:(nrow(trates)))
          nrates<- c(nrates, c(paste(cc , sep="")), "std-err")
      }
      ttt <- numeric()
      tmu <- ultimate$mu[ ,(2*bb):(2*bb+1)]
      fmu <- numeric()
      for(cc in 1:nrow(tmu))
        ttt<- c(ttt, tmu[cc,])
      fmu <- ttt
      nmu <- vector()
      if(!identical(mu,FALSE))
        for(cc in 1:length(tplotmu))
          nmu<- c(nmu, sort(tplotmu)[cc], "std-err")
      ISSRPROB <- FALSE
      if(identical(rootprobability, TRUE))
      {
        ISSRPROB <-TRUE
        nrprob <- numeric()
        for(cc in 1:length(lvls))
          nrprob <-c(nrprob,  c( lvls[cc],"std-err") )
        frprob <-(ultimate$rprob)[bb,]
        ttt <- numeric()
        tsrprob <- ultimate$srprob[ ,(2*bb):(2*bb+1)]
        fsrprob <- numeric()
        for(cc in 1:nrow(tsrprob))
          ttt<- c(ttt, tsrprob[cc,])
        fsrprob <- ttt
        nsrprob <- vector()
        for(cc in 1:length( ultimate$srprob[,1]))
          nsrprob <-c(nsrprob, c( levels(ultimate$srprob[,1])[cc], "std-err"))
      }
      else
      {
        frprob <- ultimate$rprob[ ,1+(bb)]
        fsrprob <- 0
        nsrprob <-0
        nrprob <- lvls
      }
      flog <- ultimate$loglik[bb]
      nlog <- " "
      falpha <- 0
      ISALPHA <-FALSE
      if( !identical(alpha,FALSE))
        ISALPHA <- TRUE
      if(!identical(alpha, FALSE))
      {
        falpha <- ultimate$alpha[ ( 2*bb-1):(2*bb) ]
        nalpha <- c("","")
      }
      ISMU<-FALSE
      if(!identical(mu,FALSE))
        ISMU <-TRUE
      fmu <- unlist(fmu)
      tlist[[bb]] <- c(frates, rep(c( fmu[1:(2*length(tplotmu))]), ISMU),  frprob, rep(c( fsrprob), FALSE),  rep(c( falpha), ISALPHA),  flog)
      tlist[[bb]] <- data.frame(bb,  data.frame(t(as.matrix(tlist[[bb]] ))))
      if(identical(mu, TRUE))
        colnames(      tlist[[bb]]) <- c("Site","rate_index:1", nrates[-1],rep(c(paste("Mu:",nmu[1], sep= ""), nmu[-1] ),ISMU), paste("Prior_Prob:", nrprob[1] ,sep=""), 
                                         nrprob[-1], rep(c("   Parameter estimates:", nsrprob) , FALSE), rep(c("alpha: ", " "), ISALPHA), "Log Likelihood")
      else
        colnames(      tlist[[bb]]) <- c("Site", "rate_index:1", nrates[-1],rep(c( paste("Mu:",nmu[1], sep=""), nmu[-1] ),ISMU), paste("Prior_Prob:",nrprob[1], sep = ""), 
                                         nrprob[-1], rep(c("   Parameter estimates:", nsrprob) , FALSE), rep(c("alpha: ", " "), ISALPHA), "Log Likelihood")  
      
      dlist[[bb]] <-   data.frame(Site=bb, t(as.matrix(formatC(as.numeric(as.matrix(tlist[[bb]][,-1]) ), format='f', digits= 6))))
      if(identical(mu, TRUE))
        colnames(      dlist[[bb]]) <- c("Site","rate_index:1", nrates[-1],rep(c(paste("Mu:",nmu[1], sep= ""), nmu[-1] ),ISMU), paste("Prior_Prob:", nrprob[1] ,sep=""), 
                                         nrprob[-1], rep(c("   Parameter estimates:", nsrprob) , FALSE), rep(c("alpha: ", " "), ISALPHA), "Log Likelihood")
      else
        colnames(      dlist[[bb]]) <- c("Site", "rate_index:1", nrates[-1],rep(c( paste("Mu:",nmu[1], sep=""), nmu[-1] ),ISMU), paste("Prior_Prob:",nrprob[1], sep = ""), 
                                         nrprob[-1], rep(c("   Parameter estimates:", nsrprob) , FALSE), rep(c("alpha: ", " "), ISALPHA), "Log Likelihood")  
    }
    ufinal <- do.call(rbind, tlist)
    ufinal <- cbind(1:NROW, ufinal)
    ufinal <- as.data.frame(ufinal, row.names = NULL)
    names(tlist) <- paste("Site", 1:NROW)
    names(dlist) <- paste("Site", 1:NROW)
    ultimate$dlist <- dlist
    if(identical(mu, TRUE))
      colnames(ufinal) <- c("Sites", "rate_index:1", nrates[-1],rep(c(paste("Mu:",nmu[1], sep= ""), nmu[-1] ),ISMU), paste("Prior_Prob:", nrprob[1] ,sep=""), 
                            nrprob[-1], rep(c("   Parameter estimates:", nsrprob) , FALSE), rep(c("alpha: ", " "), ISALPHA), "Log Likelihood")
    else
      colnames(ufinal) <- c("Sites", "rate_index:1", nrates[-1],rep(c( paste("Mu:",nmu[1], sep=""), nmu[-1] ),ISMU), paste("Prior_Prob:",nrprob[1], sep = ""), 
                            nrprob[-1], rep(c("   Parameter estimates:", nsrprob) , FALSE), rep(c("alpha: ", " "), ISALPHA), "Log Likelihood")
    ufinal[ which(names(ufinal)!="std-err")][  suppressWarnings(is.na(   ufinal[ which(names(ufinal)!="std-err")    ])  ) ]    <- rep("", length(    ufinal[ which(names(ufinal)!="std-err")][  suppressWarnings(is.na(   ufinal[ which(names(ufinal)!="std-err")    ]) )  ]  ))    
    ultimate$total <- ufinal
    names(tlist) <- NULL
    ultimate$total <- tlist
    end <- proc.time() - start
    ultimate$time <- end[3]
    return(ultimate)
  }
  }
