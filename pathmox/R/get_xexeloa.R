#' @title Core algorithm for pathmox and techmox
#' 
#' @description
#' Internal function. \code{get_xexeloa} is called by \code{pathmox}
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
get_xexeloa <-
function(pls, DT, EXEV, type.exev, elemnod, nv, size, mox)
{
  # cuauhmaitl elems
  N = nrow(EXEV)
  if (nv == 0) {
    elems = 1:length(elemnod) 
  } else {
    elems = which(elemnod == nv)
  }
  # identificador
  indelem = cbind(1:length(elems), elems)
  # elems en cuauhmaitl actual 
  E.node = EXEV[elems,]      
  exevs = apply(E.node, 2, function(x) nlevels(as.factor(x)))
  # categs a.cuauhmaitl
  list.cata = as.list(1:length(exevs))
  # categs b.cuauhmaitl
  list.catb = as.list(1:length(exevs))
  # + chingon global
  Ftest.global = matrix(NA, length(exevs), 5)
  # + chingon partial
  Ftest.partial = as.list(1:length(exevs))
  
  # from plspm
  Y.lvs = pls$scores
  IDM = pls$model$IDM
  
  # size limite
  if (size < 1)  size.limit = ceiling(size*N)  else  size.limit = size
  
  # for each exev
  for (i in 1:length(exevs))    
  {
    if (exevs[i] == 1)  # ya se uso
    {
      Ftest.global[i,4] <- 1   # sta bien chafa
      Ftest.partial[[i]] <- matrix(1, 1, 5) # sta bien chafa
      list.cata[[i]] <- NA# NAs a.cuauhmaitl
      list.catb[[i]] <- NA# NAs b.cuauhmaitl
      next   
    }  else   # se usa               
    {
      ### moxexeloa
      cat.exe <- unique(E.node[,i])  
      if (type.exev[i] == "ord") { 
        bin.split <- get_ordinal_split(sort(cat.exe))  
      } else { # nom
        if (exevs[i] == 2) {  
          bin.split <- as.list(cat.exe)
        } else { 
          bin.split <- get_nominal_split(cat.exe)  
        }
      }
      ### sanchez-aluja moxexeloa
      Ftest.glo <- matrix(NA, length(bin.split[[1]]), 5)   # global
      Ftest.par <- as.list(1:length(bin.split[[1]]))      # partial
      for (aux in 1:length(bin.split[[1]]))   # p/c moxexeloa
      {
        # a.cuauhmaitl
        split.a <- which(E.node[,i] %in% bin.split[[1]][[aux]])
        # b.cuauhmaitl
        split.b <- which(E.node[,i] %in% bin.split[[2]][[aux]])
        n.a <- length(split.a)
        n.b <- length(split.b)
        # cuahuitl aturat
        if (n.a<=size.limit || n.b<=size.limit)
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
      ### el gallo dels exevs
      if (mox == "pathmox")
        min.p <- which(Ftest.glo[,4] == min(Ftest.glo[,4]))[1]  
      if (mox == "techmox") 
      {
        # geometric mean
        geom_mean = lapply(Ftest.par, function(x) exp(mean(log(x[,4]))) )
        pvals.ave = unlist(geom_mean)
        min.p <- which(pvals.ave == min(pvals.ave))[1]
      }
      Ftest.global[i,] <- Ftest.glo[min.p,]
      Ftest.partial[[i]] <- Ftest.par[[min.p]]
      list.cata[[i]] <- bin.split[[1]][[min.p]]
      list.catb[[i]] <- bin.split[[2]][[min.p]]
    }
  }
  
  ### els chafillas
  if (mox == "pathmox")
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
      cata.otras <- rbind(
        cata.otras, 
        paste(as.vector(list.cata[[otras[i]]]), sep="", collapse="/")
      )
      catb.otras <- rbind(
        catb.otras, 
        paste(as.vector(list.catb[[otras[i]]]), sep="", collapse="/")
      )
    }
    test.otras <- data.frame(variable = as.vector(colnames(EXEV)[otras]), 
                             F.otras, 
                             categ.a = cata.otras, 
                             categ.b = catb.otras)
    test.otras <- test.otras[which(!is.na(F.otras[,1])),]
    test.otras[,5] <- format(test.otras[,5], scientific = FALSE)
    colnames(test.otras)[6] <- "sa/sb"
  }
  if (mox == "techmox")
  {
    # geometric mean
    geom_mean = lapply(Ftest.partial, function(x) exp(mean(log(x[,4]))) )
    pvals.ave <- unlist(geom_mean)
    otras <- order(pvals.ave) 
    cata.otras <- NULL
    catb.otras <- NULL
    for (i in 1:length(otras))
    {
      cata.otras <- rbind(
        cata.otras, 
        paste(as.vector(list.cata[[otras[i]]]), sep="", collapse="/")
      )
      catb.otras <- rbind(
        catb.otras, 
        paste(as.vector(list.catb[[otras[i]]]), sep="", collapse="/")
      )
    }
    test.otras <- data.frame(variable = as.vector(colnames(EXEV)[otras]), 
                             pval.geomean = sort(pvals.ave),
                             categ.a = cata.otras, 
                             categ.b = catb.otras)
    test.otras <- test.otras[which(!is.na(Ftest.global[,1])),]
    test.otras[,2] <- format(test.otras[,2], scientific=FALSE)
  }
  
  ### el chingon
  geomean_chingon = lapply(Ftest.partial, function(x) exp(mean(log(x[,4]))))
  pvals.ave <- unlist(geomean_chingon)
  if (mox == "pathmox")
    optim <- which(Ftest.global[,4] == min(Ftest.global[,4]))[1]
  if (mox == "techmox")
    optim <- which(pvals.ave == min(pvals.ave))[1]
  f.optim <- Ftest.global[optim,]
  cats.optim <- list.cata[[optim]]
  anticat.optim <- list.catb[[optim]]
  aux.a <- which(E.node[,optim] %in% cats.optim)
  aux.b <- which(E.node[,optim] %in% anticat.optim)
  part.a <- indelem[aux.a,2]
  part.b <- indelem[aux.b,2]
  list.elems <- list(part.a, part.b)
  
  ### resul
  test.global <- data.frame(exev = colnames(EXEV)[optim], 
                            f.stat = f.optim[1], 
                            df.num = f.optim[2],
                            df.den = f.optim[3], 
                            p.val = f.optim[4], 
                            s1_s2 = round(f.optim[5],4))
  colnames(test.global)[6] <- "sa/sb"
  test.partial <- list(exev = colnames(EXEV)[optim], 
                       p.val = pvals.ave[optim], 
                       F.test = Ftest.partial[[optim]])
  categs <- list(as.vector(cats.optim), as.vector(anticat.optim))
  # output
  list(inner.global = test.global, 
       inner.partial = test.partial, 
       categs = categs, 
       otras = test.otras, 
       list.elems = list.elems)
}