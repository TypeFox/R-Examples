#'@export
fix.techmox <-
function(pls, EXEV, X=NULL, signif=.05, size=.10, deep=2, tree=TRUE)
{ 
  # =========================== ARGUMENTS ======================================
  # pls: an object of class "plspm"
  # EXEV: data frame of segmentation variables
  # X: optional matrix or dataframe used when NO pls$data is available
  # signif: significance threshold of the F-test p-value
  # size: numeric value indicating the minimum number of elements inside a node
  # deep: an integer inidicating the tree's depth 
  # ============================================================================
  
  # ==================== Checking function arguments ===================
  if (class(pls)!="plspm")
    stop("\nSorry, 'pls' must be an object of class 'plspm'")
  if (!is.null(X)) # if X available
  {
    if (is.null(pls$data))
    {
      if (!is.matrix(X) && !is.data.frame(X))
        stop("\nSorry, 'X' must be a numeric matrix or data frame.")
      if (nrow(X)!=nrow(pls$latents))
        stop("\nOops, 'pls' and 'X' are incompatible. Different number of rows.")
      if (nrow(X)!=nrow(EXEV))
        stop("\nOops, 'X' and 'EXEV' are incompatible. Different number of rows")
    }
  } else { # if no X
    if (is.null(pls$data)) {
      stop("\nSorry, 'X' is missing. No dataset available.")
    } else {
      if (nrow(pls$data)!=nrow(EXEV)) 
        stop("\nOops, 'pls' and 'EXEV' are incompatible. Different number of rows")
    }
  }
  if (!is.data.frame(EXEV)) 
    stop("\nSorry, 'EXEV' must be a data frame containing factors")
  for (j in 1:ncol(EXEV))
    if (!is.factor(EXEV[,j]))
      stop("\nOne or more columns in 'EXEV' are not factors")    
  if (mode(signif)!="numeric" || length(signif)!=1 || signif<=0 || signif>=1)
  {
    cat("NOTICE: Invalid argument 'signif'. Default value 0.05 was used", "\n")
    signif <- 0.05
  }
  if (mode(size)!="numeric" || length(size)!=1 || size<=0)
  {
    cat("NOTICE: Invalid argument 'size'. Default value 0.10 was used", "\n")
    size <- 0.10
  }        
  if (mode(deep)!="numeric" || length(deep)!=1 || deep<1 || (deep%%1)!=0)
  {
    cat("NOTICE: Invalid argument 'deep'. Default value 2 was used", "\n")
    deep <- 2
  }           
  
  # ========================== INPUTS SETTING ==========================
  N <- nrow(EXEV)
  IDM <- pls$model$IDM
  blocks <- pls$model$blocks
  modes <- pls$model$modes
  scheme <- pls$model$scheme
  scaled <- pls$model$scaled
  tol <- pls$model$tol
  iter <- pls$model$iter
  outer <- pls$model$outer
  blocklist <- outer
  for (k in 1:length(blocks))
    blocklist[[k]] <- rep(k,blocks[k])
  blocklist <- unlist(blocklist)
  # data matrix DT
  if (!is.null(pls$data)) {
    DT <- pls$data
  } else {         
    # building data matrix 'DT'
    DT <- matrix(NA, nrow(pls$latents), sum(blocks))
    for (k in 1:nrow(IDM))
      DT[,which(blocklist==k)] <- as.matrix(X[,outer[[k]]])
    dimnames(DT) <- list(rownames(pls$latents), names(pls$out.weights))
  }
  list.nodes <- NULL
  FT.global <- NULL
  FT.partial <- NULL
  FT.labs <- NULL
  FT.otras <- NULL
  type.exev <- rep(0,ncol(EXEV))
  treat.exev <- rep("binary",ncol(EXEV))
  for (i in 1:length(type.exev))
  {
    type.exev[i] <- ifelse(is.ordered(EXEV[,i]), "ord", "nom")
    if (nlevels(EXEV[,i])>2)
      if (is.ordered(EXEV[,i]))  treat.exev[i]="ordinal"  else  treat.exev[i]="nominal"
  }
  df.exev <- data.frame(Nlevels=unlist(lapply(EXEV,nlevels)), 
                        Ordered=unlist(lapply(EXEV,is.ordered)), Treatment=treat.exev)
  
  # ========================= XEXELOA NELHUATL =========================
  elements <- 1:nrow(EXEV)   # nelhuatl elems
  elemnod <- rep(0, N)
  nv <- 0
  spli <- get_fix_xexeloa(pls, DT, EXEV, type.exev, elemnod, nv, size, "techmox")
  nodes <- 2:3# nodes 2 and 3
  parent <- c(1,1)# parent nodes
  nodes.level <- rep(1,2)# nodes level
  #    FT.global <- rbind(FT.global, spli$inner.test)# global
  FT.partial <- c(FT.partial, list(spli$inner.test$F.test)) # partial
  FT.labs <- c(FT.labs, spli$inner.test$exev)
  FT.otras <- c(FT.otras, list(spli$otras))# otras
  exvar <- rep(spli$inner.test[1], 2)
  best <- spli$inner.test$p.val     
  if (deep==1 || best>signif) {
    type.node=c("leaf","leaf")
    terminal.node <- c("yes", "yes")
  } else {
    type.node=c("node","node")  
    terminal.node <- c("no", "no")
  }
  lef.split <- spli$list.elems[[1]]# a.cuauh
  rig.split <- spli$list.elems[[2]]# b.cuauh
  list.nodes <- list(lef.split, rig.split)
  size.node <- c(length(lef.split), length(rig.split))
  elemnod[lef.split] <- 2# label elems a.cuauh
  elemnod[rig.split] <- 3# label elems b.cuauh   
  categories <- list(spli$categs[[1]], spli$categs[[2]]) 
  
  # ================ dibuix nelhuatl =====================
  #dev.new()
  op = par(mar = c(.4, .4, 1, 1.5))
  openplotmat()
  text(.25, .95, c("Fix Techmox"), cex=1, col="gray30")
  text(.75, .95, c("Supervised Tree"), cex=1, col="gray30")
  num.levels <- rep(1, deep+1)
  for (i in 1:deep)
    num.levels[i+1] <- 2^i
  elpos <- coordinates(num.levels) 
  fromto <- cbind(c(1,1), c(2,3))
  arrpos <- matrix(ncol=2, nrow=2)
  ## dibuix nelhuatl i ohtli
  arrpos[1,]<- straightarrow(to=elpos[2,], from=elpos[1,],
                             lty=2, lwd=1, lcol="grey", arr.pos=0.6, arr.length=0)
  arrpos[2,]<- straightarrow(to=elpos[3,], from=elpos[1,],
                             lty=2, lwd=1, lcol="grey", arr.pos=0.6, arr.length=0)
  textellipse(elpos[1,], 0.045, 0.045, lab=c("Root",N), box.col="#eeeeee",  
              col="#757575", lcol="#eeeeee", lwd=1, shadow.size=0, cex=0.8)    
  ## dibuix els cuauhmaitl
  for (i in 2:3)
  {
    nodlab <- c(paste("Node",nodes[i-1]), size.node[i-1])
    if (type.node[i-1]=="node") {
      textellipse(elpos[i,], 0.05, 0.03, lab=nodlab, box.col="#feb769",
                  lcol="#feb769", shadow.size=0, col="#555555", cex=.7)
    } else { # "leaf"
      textrect(elpos[i,], 0.045, 0.025, lab=nodlab, box.col="#93c4e5",
               lcol="#93c4e5", shadow.size=0, col="#555555", cex=.7)
    }
  }
  ## afegir seg vars i pvals
  x1 <- (arrpos[1,1] + arrpos[2,1]) / 2
  text(x1, arrpos[1,2], as.vector(unlist(exvar[2])), cex=.7, col="#2cb0a7")
  text(x1, arrpos[1,2], paste("pv.gm=",round(as.numeric(best),4),sep=""), 
       cex=.9*.7, col="#2cb0a7", pos=1) 
  ## afegir categs
  varcat <- NULL
  for (h1 in 1:length(categories))
    varcat <- rbind(varcat, paste(as.vector(categories[[h1]]),sep="",collapse="/"))    
  for (i in 1:2)
  {
    posi <- i+1        
    seg.cat <- as.character(varcat[i])
    seg.cat <- unlist(strsplit(seg.cat,"/"))
    if (posi%%2==0) {
      for (h in 1:length(seg.cat))
        text(arrpos[i,1]-0.03, arrpos[i,2]+h/55, seg.cat[h], cex=.7, col="#555555")
    } else {
      for (h in 1:length(seg.cat))
        text(arrpos[i,1]+0.03, arrpos[i,2]+h/55, seg.cat[h], cex=.7, col="#555555")
    }
  }
  
  ### nelhuatl xexeloa
  if (size<1)  size.limit=ceiling(size*N)  else  size.limit=size
  nivel <- 1   
  ### cuahuitl ohtli
  repeat
  {
    if (deep==1 || best>signif)
      break
    num.nods.niv <- 2^nivel   # nodos q ai n teoria
    theo.nods.niv <- (2^nivel):(2^(nivel+1)-1)  # nodos q ai n teoria
    nods.niv <- intersect(nodes, theo.nods.niv)   # nodos q ai
    for (nd in 1:length(nods.niv))   # p/c cuauhmaitl
    {
      nv <- nods.niv[nd]   # su etiketado
      indivs <- which(elemnod==nv)
      # cuahuitl aturat
      n.size <- length(indivs)
      if (n.size <= size.limit)# aturat
      {
        type.node[which(nodes==nv)] <- "leaf"     
        terminal.node[which(nodes==nv)] <- "yes"     
        next   
      } else# segueix
      {
        DT.set <- DT[indivs,]
        x.pls <- get_pls_basic(DT.set, IDM, blocks, pls$model$specs)
        # xexeloa
        spli <- get_fix_xexeloa(x.pls, DT.set, EXEV, type.exev, elemnod, nv, size, "techmox")    
        best <- spli$inner.test$p.val     
        # evaluate p-value significance stop condition
        if (best > signif)   # aturat
        {
          type.node[which(nodes==nv)] <- "leaf"     
          terminal.node[which(nodes==nv)] <- "yes"     
        } else   
        {
          #  FT.global <- rbind(FT.global, spli$inner.test)
          FT.partial <- c(FT.partial, list(spli$inner.test$F.test))
          FT.labs <- c(FT.labs, spli$inner.test$exev)
          FT.otras <- c(FT.otras, list(spli$otras))
          exvar <- c(exvar, rep(spli$inner.test[1],2))
          nodes <- c(nodes, 2*nv, 2*nv+1)   
          parent <- c(parent, rep(nv,2))
          nodes.level <- c(nodes.level, rep(nivel+1,2))
          type.node <- c(type.node, rep("node",2))   
          terminal.node <- c(terminal.node, rep("no",2))   
          lef.split <- spli$list.elems[[1]]
          rig.split <- spli$list.elems[[2]]
          ultim <- length(type.node)
          if (length(lef.split) <= size.limit)
          {
            type.node[ultim-1] <- "leaf"
            terminal.node[ultim-1] <- "yes"
          }
          if (length(rig.split) <= size.limit)
          {
            type.node[ultim] <- "leaf"
            terminal.node[ultim] <- "yes"
          }
          list.nodes <- c(list.nodes, list(lef.split, rig.split))
          siz.nd <- c(length(lef.split), length(rig.split))
          size.node <- c(size.node, siz.nd)
          elemnod[lef.split] <- 2*nv
          elemnod[rig.split] <- 2*nv+1
          cat <- length(categories)
          categories[[cat+1]] <- spli$categs[[1]]
          categories[[cat+2]] <- spli$categs[[2]]
          
          # =============================== DIBUIX ================================
          ## dibuixar ohtli
          fromto <- cbind(c(nv,nv), c(2*nv,2*nv+1))
          arrpos <- matrix(ncol=2, nrow=2)
          arrpos[1,]<- straightarrow(to=elpos[2*nv,], from=elpos[nv,],
                                     lty=2, lwd=1, lcol="grey", arr.pos=0.6, arr.length=0)
          arrpos[2,]<- straightarrow(to=elpos[2*nv+1,], from=elpos[nv,],
                                     lty=2, lwd=1, lcol="grey", arr.pos=0.6, arr.length=0)
          ## dibuix cuauhmaitl                        
          for (i in (ultim-1):(ultim))
          {                             
            parentlab <- c(paste("Node",nodes[nv-1]), size.node[nv-1])
            textellipse(elpos[nodes[nv-1],], 0.05, 0.03, lab=parentlab, box.col="#feb769",
                        lcol="#feb769", shadow.size=0, col="#555555", cex=.7)
            nodlab <- c(paste("Node",nodes[i]), size.node[i])
            if (type.node[i]=="node") {
              textellipse(elpos[nodes[i],], 0.05, 0.03, lab=nodlab, box.col="#feb769",
                          lcol="#feb769", shadow.size=0, col="#555555", cex=.7)
            } else { # "leaf"
              textrect(elpos[nodes[i],], 0.045, 0.025, lab=nodlab, box.col="#93c4e5",
                       lcol="#93c4e5", shadow.size=0, col="#555555", cex=.7)
            }
          }
          ## afegir exvars i pval
          x1 <- (arrpos[1,1] + arrpos[2,1]) / 2
          text(x1, 1.15*arrpos[1,2], as.vector(unlist(exvar[ultim])), cex=.7, col="#2cb0a7")
          text(x1, 1.15*arrpos[1,2], paste("pv.gm=",round(as.numeric(best),4),sep=""), 
               cex=.9*.7, col="#2cb0a7", pos=1) 
          ## afegir categs
          varcat <- NULL
          for (h1 in (ultim-1):(ultim))
            varcat <- rbind(varcat, paste(as.vector(categories[[h1]]),sep="",collapse="/"))    
          for (i in 1:2)
          {
            if (i==1)  posi=ultim  else  posi=ultim+1
            seg.cat <- as.character(varcat[i])
            seg.cat <- unlist(strsplit(seg.cat,"/"))
            if (posi%%2==0) {
              for (h in 1:length(seg.cat))
                text(arrpos[i,1]-0.03, arrpos[i,2]+h/55, seg.cat[h], cex=.7, col="#555555")
            } else {
              for (h in 1:length(seg.cat))
                text(arrpos[i,1]+0.03, arrpos[i,2]+h/55, seg.cat[h], cex=.7, col="#555555")
            }
          }
          # =========== s'acaba el dibuix
          
        }
      }
    }
    nivel <- nivel + 1 
    if (nivel >= deep) 
      break
    nods.nxt.niv <- (2^nivel):(2^(nivel+1)-1)# nodos q ai n teoria
    nnn <- which(nodes %in% nods.nxt.niv)  # (nnn: nods.nxt.niv)
    if (length(nnn)==0)   # aturat cullons!
      break
  }
  # reset par
  par(op)
  # end repeat
  cat("\n")
  cat("End Techmox Supervised", "\n")
  
  # cuauhmaitl terminal ultim xexeloa 
  nods.nxt.niv <- (2^nivel):(2^(nivel+1)-1)
  nnn <- which(nodes %in% nods.nxt.niv)
  if (length(nnn)>0)
    terminal.node[nnn] <- rep("yes",length(nnn))    
  excat <- NULL
  for (h1 in 1:length(categories))
    excat <- rbind(excat, paste(as.vector(categories[[h1]]),sep="",collapse="/"))    
  
  names(list.nodes) <- nodes
  nodes <- c(1, nodes)
  parent <- c(0, parent)
  nodes.level <- c(0, nodes.level)
  type.node <- c("root", type.node)
  terminal.node <- c("no", terminal.node)
  size.node <- c(N, size.node)
  Perc <- 100*(size.node/N)    
  Variable <- c(NA, as.vector(unlist(exvar)))
  Categ <- c(NA, excat)
  model <- list(mox="techmox", signif=signif, size=size, deep=deep, df.exev=df.exev)
  node.nums <- unique(parent)[-c(1,2)]
  node.labs <- c("Root", paste(rep("Node",length(node.nums)),node.nums,sep="") )   
  names(FT.otras) <- node.labs
  names(FT.partial) <- paste(node.labs, FT.labs, sep="_")
  #    FT.global[,5] <- format(FT.global[,5], scientific=FALSE)
  MOX <- data.frame( Node=nodes, Parent=parent, Depth=nodes.level, Type=type.node,  
                     Terminal=terminal.node, Size=size.node, Percent=Perc, 
                     Variable=Variable, Category=Categ )
  mox.resul <- list(MOX = MOX, 
                    FT = FT.partial, 
                    candidates = FT.otras, 
                    list.nodes = list.nodes, 
                    model = model)
  class(mox.resul) <- "treemox"
  if (tree) plot.treemox(mox.resul)
  mox.resul
}

