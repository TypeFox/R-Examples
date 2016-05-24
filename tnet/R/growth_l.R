`growth_l` <-
function(net,perspective="actor",effects,window=NULL,binary=FALSE,nstrata=10,seed=NULL,regression=TRUE){
  # Check that an actor perspective is chosen
  if(perspective!="actor") 
    stop('not yet implemented');
  # Check if a window is used, and if triadic.closure.w.gm is part of the effect list. 
  if(!is.null(window) & "triadic.closure.w.gm"%in%effects) 
    warning('Issues when using a window and the geometric triadic closure function. Treat as experimental!');
  # If seed is set, formally set it
  if(!is.null(seed))
    set.seed(as.integer(seed))
  # Sample function with accurate for length=1
  sample.accurate <- function(x, ...) x[sample.int(length(x), ...)]

cat("1/3: Preparing data\n")
  if(is.null(attributes(net)$tnet))
    net <- as.tnet(net, type="longitudinal tnet")
  if(attributes(net)$tnet!="longitudinal tnet")
    stop("Network not loaded properly")
  # Remove reinforcement/weakening of ties in binary networks
  if(binary) 
    net <- net[!duplicated(net[,c("i","j","w")]),]
  # Add window to data, and remove timestamps (as.tnet orders the data already)
  if(!is.null(window)) 
    net <- add_window_l(net=net, window=window, remove.nodes=TRUE)
  net <- net[,c("i","j","w")]
  # Total number of nodes
  N <- max(c(net[,"i"],net[,"j"]))
  # Verify that nstrata are within boundaries
  if(!is.numeric(nstrata) | nstrata < 2 | nstrata > N | nstrata != as.integer(nstrata))
    stop("nstrata not properly set: must be equal or greater than 2, less than N, and an integer value")
  # Number of attributes/covariates included
  ax   <- as.integer(0);
  # Loop for all non-network measures
  all.net.effects <- c("indegree","instrength","reciprocity","reciprocityw","triadic.closure","triadic.closure.w.min","triadic.closure.w.gm","reinforcement")
  net.effects <- effects[effects %in% all.net.effects]
  non.net.effects <- effects[!effects %in% all.net.effects]
  for (t in non.net.effects) {
    # Make sure the effect has a same., simi., or dyad. prefix. 
    if(substr(t,1,5)!="same." & substr(t,1,5)!="simi." & substr(t,1,5)!="dyad.") 
      stop(paste("effect '", t, "' is not a valid effect", sep=""))
    #Load attribute
    attrx <- substr(t, 6,nchar(t))
    attrx <- eval(parse(text=attrx));
    # Dyad attributes
    if(substr(t,1,5)=="dyad.") {
      # Check dimensions
      if(nrow(attrx) != N | ncol(attrx) != N) 
        stop('node attributes (dimentions)');
      # Include dyad attributes as is
      amat <- attrx;
    # Nodal attribtes
    } else {
      # Check dimensions
      if(length(attrx) != N)
        stop('node attributes (dimentions)')
      # Check for missing data
      if(sum(is.na(attrx))>0)
        stop('node attributes (missing data)');
      # Check for numeric data
      attrx <- as.numeric(attrx)
      if(sum(is.na(attrx))>0)
        stop('node attributes (non-numeric data)');
      # Create dyadic term matrix of nodal attribute
      amat <- matrix(data=NaN, nrow=N, ncol=N)
      if(substr(t,1,5)=="same.") {
        rule <- "as.integer(attrx[i]==attrx[j])"
      } else {
        rule <- "1-(abs(attrx[i]-attrx[j])/r)"
        r <- max(attrx)-min(attrx)
      }
      for(i in 1:(N-1))
        for(j in (i+1):N)
          amat[i,j] <- eval(parse(text=rule))
      amat[lower.tri(amat)] <- t(amat)[lower.tri(amat)]
    }
    # Add name of effect to object
    attributes(amat)$nameeffect <- t
    # Increase the count of attribute variable by 1
    ax <- ax+1;
    # Rename the dyadic matrix as a.X
    eval(parse(text=paste("a.", ax, " <- amat;", sep="")))  
    # Clean up
    rm(amat,attrx) 
  }
  # Add sequence numbers to the edgelist and a logical column to signal whether a tie is part of the strata
  net <- cbind(seq=1:nrow(net), net, strata=(net[,"i"] != net[,"j"] & net[,"w"]==1));
  rownames(net) <- NULL;
  # Create regression table, add dependent variable
  o <- data.frame(tie=rep(1:nrow(net), each=nstrata), i=rep(net[,"i"], each=nstrata), j=0, w=rep(c(as.logical(TRUE), rep(as.logical(FALSE),(nstrata-1))), nrow(net)))
  # Create matrix to store the observed and control nodes' id
  strata <- matrix(data=0, nrow=nstrata, ncol=nrow(net))

cat("2/3: Calculating measures\n")
  # Add network measures to regression table
  # Create adjacency matrix
  d <- matrix(data=as.integer(0), nrow=N, ncol=N);
  diag(d) <- NA;
  # Create node list
  n <- data.frame(node=1:N, active=as.logical(FALSE))
  # Vector containing functions to calculate the network measures
  effectfunction <- NULL
  # Add to function vector updates to adjacency matrix
  updatefunctionba <- "d[i,j]<-1;"
  updatefunctionbs <- "d[i,j]<-0;"
  updatefunctionwa <- "d[i,j]<-d[i,j]+1;"
  updatefunctionws <- "d[i,j]<-d[i,j]-1;"
  # Add the included effects to function vector
  ## Reinforcement
  if("reinforcement" %in% effects) {
    reinforcement <- matrix(data=as.integer(0), nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "reinforcement[,t] <- as.integer(d[i,hs])")
  }
  ## In-degree
  if("indegree" %in% net.effects) {
    n <- cbind(n, indegree=0)
    indegree <- matrix(data=as.integer(0), nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "indegree[,t] <- as.integer(n[hs,'indegree'])")
    updatefunctionba <- c(updatefunctionba, "n[j,'indegree'] <- n[j,'indegree']+1")
    updatefunctionbs <- c(updatefunctionbs, "n[j,'indegree'] <- n[j,'indegree']-1")
  }
  ## In-strength
  if("instrength" %in% net.effects) {
    n <- cbind(n, instrength=0)
    instrength <- matrix(data=as.integer(0), nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "instrength[,t] <- as.integer(n[hs,'instrength'])")
    updatefunctionba <- c(updatefunctionba, "n[j,'instrength'] <- n[j,'instrength']+1")
    updatefunctionwa <- c(updatefunctionwa, "n[j,'instrength'] <- n[j,'instrength']+1")
    updatefunctionbs <- c(updatefunctionbs, "n[j,'instrength'] <- n[j,'instrength']-1")
    updatefunctionws <- c(updatefunctionws, "n[j,'instrength'] <- n[j,'instrength']-1")
  }
  ## Reciprocity (dummy)
  if("reciprocity" %in% net.effects) {
    reciprocity <- matrix(data=FALSE, nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "reciprocity[,t] <- as.logical(d[hs,i]>0)")
  }
  ## Reciprocity (w_hi)
  if("reciprocityw" %in% net.effects) {
    reciprocityw <- matrix(data=FALSE, nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "reciprocityw[,t] <- as.integer(d[hs,i])")
  }
  ## Triadic closure
  if("triadic.closure" %in% net.effects) {
    dc <- d;
    triadic.closure <- matrix(data=as.integer(0), nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "triadic.closure[,t] <- as.integer(dc[i,hs])")
    updatefunctionba <- c(updatefunctionba, "dc[which(d[,i]>0),j] <- dc[which(d[,i]>0),j]+1")
    updatefunctionbs <- c(updatefunctionbs, "dc[which(d[,i]>0),j] <- dc[which(d[,i]>0),j]-1") ##untested
  }
  ## Triadic closure weighted: triplet value = minimum
  if("triadic.closure.w.min" %in% net.effects) {
    dcwmin <- d;
    triadic.closure.w.min <- matrix(data=as.integer(0), nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "triadic.closure.w.min[,t] <- as.integer(dcwmin[i,hs])")
    updatefunctionb <- c(updatefunctionb, "dcwmin[which(d[,i]>0),j] <- dcwmin[which(d[,i]>0),j]+1")
    updatefunctionw <- c(updatefunctionw, "dcwmin[which(d[,i]>=d[i,j]),j] <- dcwmin[which(d[,i]>=d[i,j]),j]+1")
  }
  ## Triadic closure weighted: triplet value = geometric mean
  if("triadic.closure.w.gm" %in% net.effects) {
    dcwgm <- d;
    triadic.closure.w.gm <- matrix(data=as.integer(0), nrow=nstrata, ncol=nrow(net))
    effectfunction <- c(effectfunction, "triadic.closure.w.gm[,t] <- as.integer(dcwgm[i,hs])")
    updatefunctionb <- c(updatefunctionb, "dcwgm[which(d[,i]>0),j] <- dcwgm[which(d[,i]>0),j]+sqrt(d[which(d[,i]>0),i])")
    #Update function when a tie already exists between i and j
    ##1)  find the nodes connected to i (hs)
    ##2)  if j is in hs, remove it
    ##3)  if there are any hs, update the [hs,j]-dyads with
    ##4)  for each a in hs: 
    ##.1)  find the nodes that a is tied to, ls
    ##.2)  if j is in ls, remove it
    ##.3)  for each l in ls:
    ##..1)  add sqrt(d[a,l]*d[l,j]) to the [a,j]-dyad
    updatefunctionw <- c(updatefunctionw, 
      "hss <- which(d[,i]>0)
      if(length(which(hss==j))>0)
        hss <- hss[-which(hss==j)]
      for(a in hss) {
        out <- 0
        ls <- which(d[a,]>0); 
        if(length(which(ls==j))>0)
          ls <- ls[-which(ls==j)]
        for(l in ls)
          out <- out+sqrt(d[a,l]*d[l,j])
        dcwgm[a,j] <- out
      }")
  }
  # Add time indicator to data
  indicator <- rep(as.logical(TRUE),nrow(net))
  if(length(indicator)>50)
    indicator <- as.logical(c(rep(c(rep(FALSE, floor(length(indicator)/50)-1), 1),50), rep(FALSE, length(indicator)-length(rep(c(TRUE, rep(FALSE, floor(length(indicator)/50)-1)),50)))))
  # Display indicator bar
  cat(paste("0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n",
            "+----+----+----+----+----+----+----+----+----+----+\n|", sep=""));
  # Preform loop of all ties
  ## Extract tie information
  ## If s a.k.a. positive and non-self-loop
  ### h: Possible control nodes (not i or j; if binary analysis, not i's already established contacts
  ### If h is less than the number of control nodes, exclude tie from regression
  ### hs: Randomly sample h, and add observed receiver
  ### Calculate measures in function vector
  ### Run update vectors
  ## If s is a self-loop, set node as active or passive (w=1 or -1)
  ## Time indicator display
  for (t in 1:nrow(net)) {
    tie <- net[t,]
    i <- as.integer(tie[2])
    j <- as.integer(tie[3])
    w <- as.integer(tie[4])
    s <- as.logical(tie[5])
    if(s) {
      h <- n[n[,"active"],"node"]
      h <- h[!h %in% c(i,j)]
      if(binary)
        h <- h[!h %in% which(d[i,]==1)]
      if(length(h)<(nstrata-1)) {
        net[t,"strata"] <- FALSE
      } else {
        hs <- c(j, sample.accurate(h, size=(nstrata-1)))
        strata[,t] <- hs
        eval(parse(text=effectfunction))
      }
      if(d[i,j]==0) {
        eval(parse(text=updatefunctionba))
      } else {
        eval(parse(text=updatefunctionwa))
      }
    } else if(i != j) {
      if(d[i,j]==1) {
        eval(parse(text=updatefunctionbs))
      } else {
        eval(parse(text=updatefunctionws))
      }
    } else if(w==1) {
      n[i,"active"] <- TRUE
    } else {
      n[i,"active"] <- FALSE
    }
    ##Indicator
    if(indicator[t]) cat("|")
  }
  cat("\n")
  # Add receving nodes' id to regression table
  o[,"j"] <- as.vector(strata)
  # Add network measures to regression table
  for(t in effects[effects %in% net.effects])
    eval(parse(text=paste("o <- cbind(o, ", t, "=as.vector(", t, "))", sep="")))
  # Remove ties with missing control nodes
  o <- o[o[,"j"]!=0,]
  rownames(o) <- NULL
  
  # Add attribute terms to regression table
  if(ax>0) {
    for(t in 1:ax) {
      # Name of term
      tobj <- eval(parse(text=paste("attributes(a.",t, ")$nameeffect", sep="")))
      # If same(attribute), then add a logical column
      #if(substr(tobj,1,5)=="same.") {
      #  o <- eval(parse(text=paste("cbind(o, ", tobj, "=as.logical(a.", t, "[cbind(rep(net[,'i'], each=nstrata),as.vector(strata))]))", sep="")))
      # Otherwise add as-is
      #} else {
        o <- eval(parse(text=paste("cbind(o, ", tobj, "=as.vector(a.", t, "[cbind(rep(net[,'i'], each=nstrata),as.vector(strata))]))", sep="")))
      #}
    }
  }

  # Conduct regression or output regression table
  if(regression) {
    cat("3/3: Running regression\n");
    # Create regression formula
    formula <- "w ~"
    for(t in effects)
      formula <- paste(formula, t, "+")
    formula <- substr(formula, 1, (nchar(formula)-2))
    # Run regression
    reg <- eval(parse(text=paste("clogit(", formula, " + cluster(i) + strata(tie), data=o, method='approximate')", sep="")))
    #Return the regression object to user
    return(reg)
  } else {
    cat("3/3: Outputting data\n");
    #Return the regression table to user with variables
    return(o)
  }
}
