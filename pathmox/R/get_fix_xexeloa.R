#' @title Core algorithm for fixed version of pathmox and techmox
#' 
#' @description
#' Internal function. \code{get_fix_xexeloa} is called by \code{fix.pathmox}
#' and \code{fix.techmox}
#' 
#' @param pls object of class plspm
#' @param DT data table
#' @param EXEV data frame with segmentation variables
#' @param type.exev vector with types of categorical variables (ordinal, nominal)
#' @param elemnod element in node nv
#' @param nv number of node
#' @param size minimum size of elements inside a node
#' @param mox string indicating type of algorithm "pathmox" or "techmox"
#' @export
#' @keywords internal
get_fix_xexeloa <-
function(pls, DT, EXEV, type.exev, elemnod, nv, size, mox)
{
  # cuauhmaitl "nv" elems 
  N <- nrow(EXEV)
  if (nv==0) elems<-1:length(elemnod) else elems<-which(elemnod==nv)
  indelem <- cbind(1:length(elems), elems)# identificador
  E.node <- EXEV[elems,]# elems en cuauhmaitl actual      
  exevs <- apply(E.node, 2, function(x) nlevels(as.factor(x)))# vector fo exevs
  list.cata <- as.list(1:length(exevs))   # categs a.cuauhmaitl
  list.catb <- as.list(1:length(exevs))   # categs b.cuauhmaitl
  Ftest.global <- matrix(NA, length(exevs), 5)# + chingon global
  Ftest.partial <- as.list(1:length(exevs))# + chingon partial
  
  # parameters of plspm
  Y.lvs <- pls$scores# matrix of latent variables
  IDM <- pls$model$IDM# IDM matrix
  blocks <- pls$model$blocks
  modes <- pls$model$modes
  scheme <- pls$model$scheme
  scaled <- pls$model$scaled
  tol <- pls$model$tol
  iter <- pls$model$iter
  endo <- rowSums(IDM)
  endo[endo!=0] <- 1  
  n.endo <- sum(endo)
  
  # size limite
  if (size < 1)  size.limit=ceiling(size*N)  else  size.limit=size
  
  # for each exev
  for (i in 1:length(exevs))    
  {
    if (exevs[i] == 1)  # ya se uso
    {
      Ftest.global[i,4] <- 1   # sta bien chafa
      list.cata[[i]] <- NA# NAs a.cuauhmaitl
      list.catb[[i]] <- NA# NAs b.cuauhmaitl
      next   # next exev
    }  else   # se usa               
    {
      ### moxexeloa
      cat.exe <- unique(E.node[,i])  
      if (type.exev[i] == "ord") {
        bin.split <- .ordinal.split(sort(cat.exe))  
      } else { # nom
        if (exevs[i]==2) {
          bin.split <- as.list(cat.exe)
        } else { 
          bin.split <- .nominal.split(cat.exe)   # function of nominal splits
        }
      }
      ### sanchez-aluja moxexeloa
      Ftest.glo <- matrix(NA, length(bin.split[[1]]), 5)   # global
      Ftest.par <- as.list(1:length(bin.split[[1]]))      # partial
      for (aux in 1:length(bin.split[[1]]))   # p/c moxexeloa
      {
        split.a <- which(E.node[,i] %in% bin.split[[1]][[aux]])# a.cuauhmaitl
        split.b <- which(E.node[,i] %in% bin.split[[2]][[aux]])  # b.cuauhmaitl
        n.a <- length(split.a)
        n.b <- length(split.b)
        if (n.a<=size.limit || n.b<=size.limit)# cuahuitl aturat
        {
          Ftest.glo[aux,4] <- 1   # sta bien chafa
          Ftest.par[[aux]] <- matrix(1,1,5)# sta bien chafa
          next  
        } else# cuahuitl segueix
        {
          idm <- .sanchez.aluja(Y.lvs, IDM, split.a, split.b)
          Ftest.glo[aux,] <- idm[[1]]   # global
          Ftest.par[[aux]] <- idm[[2]]  # partial
        }
      }
      ### el gallo dels exev
      if (mox=="pathmox")
        min.p <- which(Ftest.glo[,4]==min(Ftest.glo[,4]))[1]  
      if (mox=="techmox") 
      {
        pvals.ave <- unlist(lapply(Ftest.par, function(x) exp(mean(log(x[,4]))) ))
        min.p <- which(pvals.ave == min(pvals.ave))[1]
      }
      Ftest.global[i,] <- Ftest.glo[min.p,]
      Ftest.partial[[i]] <- Ftest.par[[min.p]]
      list.cata[[i]] <- bin.split[[1]][[min.p]]
      list.catb[[i]] <- bin.split[[2]][[min.p]]
    }
  }
  ### els chafillas
  if (mox=="pathmox")
  {
    otras <- order(Ftest.global[,4]) 
    F.otras <- Ftest.global[otras,]
    F.otras[,4] <- round(F.otras[,4],6)
    F.otras[,5] <- round(F.otras[,5],4)
    colnames(F.otras) <- c("f.stat","df.num","df.den","p.val","sa/sb")
    cata.otras <- NULL
    catb.otras <- NULL
    for (i in 1:length(otras))
    {
      cata.otras <- rbind(cata.otras, paste(as.vector(list.cata[[otras[i]]]),sep="",collapse="/"))
      catb.otras <- rbind(catb.otras, paste(as.vector(list.catb[[otras[i]]]),sep="",collapse="/"))
    }
    test.otras <- data.frame(variable=as.vector(colnames(EXEV)[otras]), F.otras, 
                             categ.a=cata.otras, categ.b=catb.otras)
    test.otras <- test.otras[which(!is.na(F.otras[,1])),]
    test.otras[,5] <- format(test.otras[,5], scientific=FALSE)
    colnames(test.otras)[6] <- "sa/sb"
    who.otras <- otras[which(!is.na(F.otras[,1]))]
  }
  if (mox=="techmox")
  {
    pvals.ave <- unlist(lapply(Ftest.partial, function(x) exp(mean(log(x[,4]))) ))
    otras <- order(pvals.ave) 
    pvals.otras <- pvals.ave[otras]
    cata.otras <- NULL
    catb.otras <- NULL
    for (i in 1:length(otras))
    {
      cata.otras <- rbind(cata.otras, paste(as.vector(list.cata[[otras[i]]]),sep="",collapse="/"))
      catb.otras <- rbind(catb.otras, paste(as.vector(list.catb[[otras[i]]]),sep="",collapse="/"))
    }
    test.otras <- data.frame(variable=as.vector(colnames(EXEV)[otras]), pval.geomean=sort(pvals.ave),
                             categ.a=cata.otras, categ.b=catb.otras)
    test.otras <- test.otras[which(!is.na(Ftest.global[,1])),]
    test.otras[,2] <- format(test.otras[,2], scientific=FALSE)
    who.otras <- otras[which(pvals.otras!=1)]
  }
  
  ### el fijado         
  
  if (length(who.otras)>0)
  {
    candi.table <- cbind(Num=who.otras, test.otras)
    if (nv==0)
      cat("CANDIDATE SPLITS: ROOT NODE", "\n\n")
    if (nv>0)
      cat("CANDIDATE SPLITS: NODE", nv, "\n\n")
    print(candi.table)
    cat("\n")
    cat("Select a candidate split from column 'Num' and then press Enter:", "\n")
    fijo <- scan(what=integer(1), n=1)
    if (mode(fijo)!="numeric" || length(fijo)!=1 || fijo<1 || (fijo%%1)!=0) {
      cat("Invalid number of candidate split", "\n")
      stop() 
    }
    f.fijo <- Ftest.global[fijo,]
  } else {
    fijo <- which(Ftest.global[,4]==min(Ftest.global[,4]))[1]
    f.fijo <- Ftest.global[fijo,]
  }         
  cats.fijo <- list.cata[[fijo]]
  anticat.fijo <- list.catb[[fijo]]
  
  aux.a <- which(E.node[,fijo] %in% cats.fijo)
  aux.b <- which(E.node[,fijo] %in% anticat.fijo)
  part.a <- indelem[aux.a,2]
  part.b <- indelem[aux.b,2]
  list.elems <- list(part.a, part.b)
  
  ### resul
  if (mox=="pathmox")
  {
    test.fijo <- data.frame(exev=colnames(EXEV)[fijo], f.stat=f.fijo[1], df.num=f.fijo[2],
                            df.den=f.fijo[3], p.val=f.fijo[4], sa_sb=round(f.fijo[5],4))
    colnames(test.fijo)[6] <- "sa/sb"
  }
  if (mox=="techmox")
  {
    test.fijo <- list(exev=colnames(EXEV)[fijo], p.val=pvals.ave[fijo], 
                      F.test=Ftest.partial[[fijo]])
  }
  categs <- list(as.vector(cats.fijo), as.vector(anticat.fijo))
  resul <- list(inner.test=test.fijo, categs=categs, list.elems=list.elems, otras=test.otras)
  return(resul)
}
