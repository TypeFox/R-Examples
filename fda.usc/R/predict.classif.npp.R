predict.classif.npp<-function (object, newx= NULL, type = "class", ...) 
{
  #poner x en classif o lfdata??
  if (class(newx)=="list")
  new.fdataobj<-newx[[1]]

 else stop("newx not is a list object") 
 newlfdata<-newx
 if (is.null(new.fdataobj)) 
    return(object$group.est)
  isfdata <- is.fdata(new.fdataobj)
  object$group <- factor(object$group, levels = levels(object$group)[which(table(object$group) > 
                                                                             0)])
ny <- lev <- levels(object$group)
lenl<-length(newx) 
nn = nrow(new.fdataobj$data)
nam1<-names(newx)##################fdalta ver si es nulo que se creen  paste("var,ii)
y = object$y
h = object$h.opt
ny <- levels(y)
numg = nlevels(y)
C <- object$call
m <- object$m
Ker = object$Ker

if (object$C[[1]] == "classif.npp") {
  #(lfdata, lfdataref = lfdata, metric, par.metric = list(), 
  #          weights, method = "euclidean")
  par.ldata<-list()
  par.ldata$lfdataref<-object$x
  par.ldata$lfdata<-newx
  par.ldata$metric
  par.ldata$par.metric<-list()
  for (i in 1:lenl) par.ldata$par.metric[[nam1[i]]]=attributes(object$mdist)[[nam1[i]]]$par.metric
  #OJO DIVIDIR POR EL DSCALE CORREPONDIENTE

  par.ldata$weights<-attributes(object$mdist)$weights
  par.ldata$method<-attributes(object$mdist)$method
  nmdist<-do.call("metric.ldata",par.ldata)
  
  # print(nmdist[1:3,1:3])
  
  #print(object$mdist[1:3,1:3])
  #print(nn)
  object$par.S$cv<-FALSE
  object$par.S$tt <- nmdist  
    kmdist = object$type.S(nmdist, h = h, Ker = object$Ker,w=object$par.S$w,
                           cv =object$par.S$cv)
     
    pgrup = array(0, dim = c(numg, nn))
  l<- group.pred <- rep(0,length.out= nn)
    for (j in 1:numg) {
      grup = as.integer(y == lev[j])
      pgrup[j, ] <- kmdist %*% matrix(grup, ncol = 1)
    }
    group.pred <- factor(ny[apply(pgrup, 2, which.max)], 
                         levels = ny)
    pgrup <- t(pgrup)
    group.est <- numeric(nn)
    ty <- "S.KNN"
  if (ty == "S.KNN") {
    for (ii in 1:nn) {
      l = seq_along(pgrup[ii, ])[pgrup[ii, ] == max(pgrup[ii, 
                                                          ], na.rm = T)]      
      if (length(l) > 1) {
        abc<-which(nmdist[ii, ] == min(nmdist[ii, ], na.rm = T))
        
        group.est[ii] =y[abc]
      }
      else  group.est[ii] = ny[l[1]]
    }
    group.pred <-  factor(group.est,levels = ny)
  }
  }
  if (type == "class") 
    return(group.pred)
  if (type == "probs") {
    if (isfdata) 
      rownames(pgrup) <- rownames(new.fdataobj$data)
    else rownames(pgrup) <- rownames(new.fdataobj)
    colnames(pgrup) <- levels(object$group)
    return(list(group.pred = group.pred, prob.group = pgrup))
  }
  else stop("type argument should be one of class or probs")
}
############################################

