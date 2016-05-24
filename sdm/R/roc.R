# Author: Babak Naimi, naimi.b@gmail.com
# Date :  Feb. 2015
# Version 2.0
# Licence GPL v3
#-------
.roc <- function(o,p) {
  th <- as.vector(quantile(p,0:100/100))
  e <- matrix(nrow=length(th),ncol=2)
  colnames(e) <- c('sensitivity','1-specificity')
  
  for (i in seq_along(th)) {
    w <- which(p >= th[i])
    pt <- rep(0,length(p))
    pt[w] <- 1
    e[i,] <- .evaluate.cmx(.cmx(o,pt))[1:2]
  }
  e[,2] <- 1-e[,2]
  e <- rbind(e,c(0,0))
  e[,c(2,1)]
}
#------------

.getRoc <- function(x,p,...) {
  if (missing(p)) p <- NULL
  wtest <- NULL
  a <- c('species','method','replication','run','w','wtest')
  dot <- list(...)
  ndot <- names(dot)
  if (length(ndot) > 0) {
    ndot <- .pmatch(ndot,a)
    w <- !is.na(ndot)
    ndot <- ndot[w]
    dot <- dot[w]
    names(dot) <- ndot
    if (!'w' %in% ndot) {
      if (!is.null(p) && is.numeric(p) && length(p) == 1) dot[['w']] <- p
    }
    if ('wtest' %in% ndot) {
      wtest <- dot[['wtest']]
      dot <- dot[-which(ndot == 'wtest')]
      if (is.character(wtest))  wtest <- .pmatch(wtest,c("training","dep","indep"))
      else if (is.numeric(wtest) && any(wtest %in% 1:3)) wtest <- c("training","dep","indep")[wtest[wtest %in% 1:3]]
      
      wtest <- unlist(lapply(wtest,function(x) {
        if (x == 'dep') 'test.dep'
        else if (x == 'indep') 'test.indep'
        else 'training'})
      )
    } else {
      wtest <- names(x@models[[1]][[1]][[1]]@evaluation)
    }
    
  } else if (!is.null(p) && is.numeric(p) && length(p) == 1) {
    dot[['w']] <- p
    wtest <- names(x@models[[1]][[1]][[1]]@evaluation)
  }
  if (is.null(wtest)) stop('there is no evaluation in the model object!')
  names(wtest) <- wtest
  
  dot[['x']] <- x
  mi <- do.call(.getModel.info,dot)
  o <- list()
  for (sp in unique(mi[,2])) {
    for (m in unique(mi[,3])) {
      if (!is.na(unique(mi[,4])[1])) {
        for (r in unique(mi[,4])) {
          rr <- .getModel.info(x,species=sp,method=m,replication=r)[,1]
          names(rr) <- rr
          o[[sp]][[m]][[r]] <- lapply(rr,function(i) {
            lapply(wtest,function (j) {
              if (inherits(x@models[[sp]][[m]][[as.character(i)]]@evaluation[[j]],'sdmEvaluate')) {
                .roc(x@models[[sp]][[m]][[as.character(i)]]@evaluation[[j]]@observed,x@models[[sp]][[m]][[as.character(i)]]@evaluation[[j]]@predicted)
              } else NA
            })
          })
        }
      } else {
        rr <- .getModel.info(x,species=sp,method=m)[,1]
        names(rr) <- rr
        o[[sp]][[m]][['no_replication']] <- lapply(rr,function(i) {
          lapply(wtest,function (j) {
            if (inherits(x@models[[sp]][[m]][[as.character(i)]]@evaluation[[j]],'sdmEvaluate')) {
              .roc(x@models[[sp]][[m]][[as.character(i)]]@evaluation[[j]]@observed,x@models[[sp]][[m]][[as.character(i)]]@evaluation[[j]]@predicted)
            } else NA
          })
        })
      }
    }
  }
  o
}


#-----------
.plot.roc <- function(x,y=NULL,auc=NULL,smooth=TRUE,wtest=NULL,...) {
  dot <- list(...)
  ndot <- names(dot)
  
  if (!'xlab' %in% ndot) dot[['xlab']] <- '1-Specificity (false positive rate)'
  if (!'ylab' %in% ndot) dot[['ylab']] <- 'Sensitivity (true positive rate)'
  if (!'main' %in% ndot) dot[['main']] <-  'ROC plot'
  if (!'xlim' %in% ndot) dot[['xlim']] <- c(0,1)
  if (!'ylim' %in% ndot) dot[['ylim']] <- c(0,1)
  if (!'cex' %in% ndot) dot[['cex']] <- 1
  if ('lwd' %in% ndot) lwd <- dot[['lwd']]
  else dot[['lwd']] <- lwd <- 1.4
  if (!'type' %in% ndot) dot[['type']] <- 'l'
  
  
  if (!'cex' %in% ndot) dot[['cex']] <- 1
  
  if ('col' %in% ndot) {
    if (length(dot[['col']]) == 1) {
      if (as.character(dot[['col']]) %in% c('red',"#FF0000","2")) cl <- c(dot[['col']],'#FFC0CB','blue','#ADD8E6')
      else if (as.character(dot[['col']]) %in% c('blue',"#0000FF","4")) cl <- c(dot[['col']],'#ADD8E6','red','#FFC0CB')
      else cl <- c(dot[['col']],'#BEBEBE','red','#FFC0CB')
    } else if (length(dot[['col']]) == 2) {
        if (as.character(dot[['col']][1]) %in% c('red',"#FF0000","2")) cl <- c(dot[['col']],'blue','#ADD8E6')
        else if (as.character(dot[['col']][1]) %in% c('blue',"#0000FF","4")) cl <- c(dot[['col']],'red','#FFC0CB')
        else cl <- c(dot[['col']],'red','#FFC0CB')
      } else if (length(dot[['col']]) == 3) cl <- c(dot[['col']],'#BEBEBE')
      else cl <- dot[['col']] <- dot[['col']][1:4]
    } else cl <- c('red','#FFC0CB','blue','#ADD8E6')
  
  ym <- xm <- NULL
  if (inherits(x,'list')) {
    xm <- x[[1]]
    if (length(x) > 1) {
      for (i in 2:length(x)) xm <- xm + x[[i]]
      xm <- xm / length(x)
    }
    
    if (smooth) {
      for (i in seq_along(x)) {
        rr <- try(supsmu(x[[i]][,1],x[[i]][,2],bass=0),silent=TRUE)
        if (!inherits(rr,'try-error')) {
          if (rr[[2]][1] != 0)  {
            rr[[1]] <- c(0,rr[[1]])
            rr[[2]] <- c(0,rr[[2]])
          }
          x[[i]] <- rr
        }
      }
      rr <- try(supsmu(xm[,1],xm[,2],bass=0),silent=TRUE)
      if (!inherits(rr,'try-error')) {
        if (rr[[2]][1] != 0)  {
          rr[[1]] <- c(0,rr[[1]])
          rr[[2]] <- c(0,rr[[2]])
        }
        xm <- rr
      }
      rm(rr)
    }
    dot[['lwd']] <- lwd * 0.5
    dot[['col']] <- cl[2]
    dot[['x']] <- x[[1]]
    do.call(plot,dot)
    if (length(x) > 1) {
      for (i in 2:length(x)) {
        lines(x[[i]],col=dot[['col']][2],lwd=dot[['lwd']]*0.7)
      }
    }
    dot[['lwd']] <- lwd
    abline(a=0,b=1,lty=2,lwd=2)
    
  } else {
    dot[['col']] <- cl[1]
    dot[['x']] <- x
    do.call(plot,dot)
    #plot(x,type='l',...)
  }
  
  if (!is.null(y)) {
    if (inherits(y,'list')) {
      ym <- y[[1]]
      if (length(y) > 1) {
        for (i in 2:length(y)) ym <- ym + y[[i]]
        ym <- ym / length(y)
      }
      
      if (smooth) {
        for (i in seq_along(y)) {
          rr <- try(supsmu(y[[i]][,1],y[[i]][,2],bass=0),silent=TRUE)
          if (!inherits(rr,'try-error')) {
            if (rr[[2]][1] != 0)  {
              rr[[1]] <- c(0,rr[[1]])
              rr[[2]] <- c(0,rr[[2]])
            }
            y[[i]] <- rr
          }
        }
        rr <- try(supsmu(ym[,1],ym[,2],bass=0),silent=TRUE)
        if (!inherits(rr,'try-error')) {
          if (rr[[2]][1] != 0)  {
            rr[[1]] <- c(0,rr[[1]])
            rr[[2]] <- c(0,rr[[2]])
          }
          ym <- rr
        }
        rm(rr)
      }
      for (i in 1:length(y)) {
        lines(y[[i]],col=cl[4],lwd=dot[['lwd']]*0.5)
      }
    } else {
      lines(y,col=cl[3],lwd=dot[['lwd']]*0.5)
    }
  }
  
  abline(a=0,b=1,lty=2,lwd=2)
  if (!is.null(xm)) lines(xm,lwd=dot[['lwd']]*2,col=cl[1])
  if (!is.null(ym)) lines(ym,lwd=dot[['lwd']]*2,col=cl[3])
  
  if (!missing(auc) && !is.null(auc)) {
    if (length(auc) > 1) {
      if (!is.null(y)) {
        auc <- auc[1:2]
        if (!is.null(wtest) & length(wtest) == 2) n <- c(paste("AUC (",wtest[1],") = ",auc[1],sep=''),paste("AUC (",wtest[2],") = ",auc[2],sep=''))
        else n <- c(paste("AUC = ",auc[1],sep=''),paste("AUC = ",auc[2],sep=''))
      }
    } else n <- paste("AUC = ",auc,sep='')
    
    if (length(x) > 1) n <- paste('Mean ',n,sep='')
    
    if (length(auc) == 1) {
      legend(x="bottomright",legend=n,lty=1,col=cl[1],cex=dot[['cex']],lwd=dot[['lwd']])
    } else {
      legend(x="bottomright",legend=n,lty=c(1,1),col=cl[c(1,3)],cex=dot[['cex']],lwd=dot[['lwd']])
    }
    
  }
}

#----------
if (!isGeneric("roc")) {
  setGeneric("roc", function(x,p=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,wtest=NULL,smooth=FALSE,legend=TRUE,...)
    standardGeneric("roc"))
}


setMethod("roc", signature(x='vector',p='vector'),
          function(x,p,...) {
            dot <- list(...)
            ndot <- names(dot)
            l <- TRUE
            if ('legend' %in% ndot) l <- dot[['legend']]
            
            a <- c('species','method','replication','run','wtest','legend')
            if (any(ndot %in% a)) {
              dot <- dot[which(ndot %in% a)]
            }
            
            r <- getRoc(x,p)
            if (l) dot[['auc']] <- .auc(x,p)
            dot[['x']] <- r
            do.call(.plot.roc,dot)
          }
)

setMethod("roc", signature(x='sdmModels'),
          function(x,p=NULL,species=NULL,method=NULL,replication=NULL,run=NULL,wtest=NULL,smooth=FALSE,legend=TRUE,...) {
            op <- par(no.readonly = TRUE) 
            w <- NULL
            if (!is.null(p)) w <- p 
            #----------
            wt <- names(x@models[[1]][[1]][[1]]@evaluation)
            
            if (!is.null(wtest)) {
              if (is.numeric(wtest)) {
                if (any(wtest %in% 1:3)) wtest <- c("training","dep","indep")[wtest[wtest %in% 1:3]]
                else wtest <- names(x@models[[1]][[1]][[1]]@evaluation)
              } else wtest <- .pmatch(wtest,c("training","dep","indep"))
              
              wtest <- unlist(lapply(wtest,function(x) {
                if (x == 'dep') 'test.dep'
                else if (x == 'indep') 'test.indep'
                else 'training'
              }))
              wt <- wtest %in% wt
              if (any(wt)) wtest <- wtest[wt]
              else stop('wtest argument is not properly specified!')
              if (length(wtest) > 2) {
                wtest <- c('test.dep','test.indep')
                warning('wtest can have 2 items; training is excluded!')
              }
              
            } else wtest <- names(x@models[[1]][[1]][[1]]@evaluation)
            
            if (length(wtest) > 2) wtest <- c('test.dep','test.indep')
            if (is.null(wtest)) stop('there is no evaluation in the model object!')
            names(wtest) <- wtest
            #----------
            
            mi <- .getModel.info2(x,w=w,species=species,method=method,replication=replication,run=run)
            rownames(mi) <- as.character(mi[,1])
            
            if (length(wtest) == 2) mi <- mi[apply(mi[,wtest],1,function(x) all(x)),]
            else mi <- mi[mi[,wtest[1]],]
            
            if (nrow(mi) == 0) stop('no successfull evaluation exist in the model object!')
            
            species <- as.character(mi[1,2])
            m <- as.character(unique(mi[,3]))
            r <- unique(mi[,4])
            #-------
            
            dot <- list(...)
            ndot <- names(dot)
            
            if (is.na(r[1])) {
              np <- length(m)
              
              if (!'main' %in% ndot) Main <- 'ROC'
              else if (is.null(dot[['main']])) Main <- NULL
              else Main <- dot[['main']]
              
              if (np > 16) {
                warning(paste('Due to larger number of plots, only the first', w ,'models are used...'))
                par(mfrow=c(4,4))
              } else {
                w <- floor(sqrt(np))
                h <- ceiling(np/w)
                if (abs(w-h) > 1) {
                  w <- w + 1
                  h <- h - 1
                }
                
                par(mfrow=c(w,h))  
              }
              
              o <- auc <- mn <- list()
              
              for (mo in m) {
                mn[[mo]] <- paste(Main,' (',mo,')',sep='')
                rr <- mi[mi[,2] == species & mi[,3] == mo,1]
                names(rr) <- rr
                o[[mo]] <- lapply(wtest,function(j) {
                  lapply(rr,function (i) {
                    if (mi[as.character(i),j]) {
                      .roc(x@models[[species]][[mo]][[as.character(i)]]@evaluation[[j]]@observed,x@models[[species]][[mo]][[as.character(i)]]@evaluation[[j]]@predicted)
                    } else NA
                  })
                })
                
                auc[[mo]] <- lapply(wtest,function(j) {
                  lapply(rr,function (i) {
                    if (mi[as.character(i),j]) {
                      x@models[[species]][[mo]][[as.character(i)]]@evaluation[[j]]@statistics$AUC
                    } else NA
                  })
                })
              }
              #-------
              if (np > 10) z <- 0.04
              else if (np > 5 & np <= 10) z <- 0.06
              else z <- 0.13
              
              for (mo in m) {
                xx <- o[[mo]][[wtest[1]]][!is.na(o[[mo]][[wtest[1]]])]
                if (length(xx) > 0) {
                  dot[['x']] <- xx
                  if (length(wtest) == 2) {
                    dot[['y']] <- o[[mo]][[wtest[2]]][!is.na(o[[mo]][[wtest[2]]])]
                    if (legend) {
                      dot[['auc']] <- c(round(mean(unlist(auc[[mo]][[wtest[1]]]),na.rm=TRUE),3),round(mean(unlist(auc[[mo]][[wtest[2]]]),na.rm=TRUE),3))
                      dot[['wtest']] <- wtest
                    }
                  } else {
                    if (legend) {
                      dot[['auc']] <-round(mean(unlist(auc[[mo]][[wtest[1]]]),na.rm=TRUE),3)
                      dot[['wtest']] <- wtest
                    }
                  }
                  
                  if ('cex' %in% ndot) dot[['cex']] <- dot[['cex']] - np*z
                  else dot[['cex']] <- 1 - np*z
                  dot[['smooth']] <- smooth
                  
                  if (!is.null(Main) && np > 1) dot[['main']] <- mn[[mo]]
                  
                  do.call(.plot.roc,dot)
                }
              }
            } else {
              np <- length(r) * length (m)
              
              if (!'main' %in% ndot) Main <- 'ROC'
              else if (is.null(dot[['main']])) Main <- NULL
              else Main <- dot[['main']]
              
              
              if (length(r) == 3 & np > 12) {
                warning('Due to larger number of plots, only the first 4 models are used...')
                m <- m[1:3]
                par(mfrow=c(4,3))
              } else if (np > 16) {
                w <- floor(16 / length(r))
                m <- m[1:w]
                warning(paste('Due to larger number of plots, only the first', w ,'models are used...'))
                par(mfrow=c(4,4))
              } else {
                w <- floor(sqrt(np))
                h <- ceiling(np/w)
                if (abs(w-h) > 1) {
                  w <- w+1
                  h <- h-1
                }
                par(mfrow=c(w,h))  
              }
              
              o <- auc <- mn <- list()
              for (mo in m) {
                mn[[mo]] <- list()
                for (re in r) {
                  mn[[mo]][[re]] <- paste(Main,' (',mo,' - ',re,')',sep='')
                  rr <- mi[mi[,2] == species & mi[,3] == mo & mi[,4] == re,1]
                  names(rr) <- rr
                  o[[mo]][[re]] <- lapply(wtest,function(j) {
                    lapply(rr,function (i) {
                      if (mi[as.character(i),j]) {
                        .roc(x@models[[species]][[mo]][[as.character(i)]]@evaluation[[j]]@observed,x@models[[species]][[mo]][[as.character(i)]]@evaluation[[j]]@predicted)
                      } else NA
                    })
                  })
                  
                  auc[[mo]][[re]] <- lapply(wtest,function(j) {
                    lapply(rr,function (i) {
                      if (mi[as.character(i),j]) {
                        x@models[[species]][[mo]][[as.character(i)]]@evaluation[[j]]@statistics$AUC
                      } else NA
                    })
                  })
                }
              }
              
              if (np > 10) z <- 0.04
              else if (np > 5 & np <= 10) z <- 0.06
              else z <- 0.13
              
              
              for (mo in m) {
                for (re in r) {
                  xx <- o[[mo]][[re]][[wtest[1]]][!is.na(o[[mo]][[re]][[wtest[1]]])]
                  if (length(w) > 0) {
                    dot[['x']] <- xx
                    if (length(wtest) == 2) {
                      dot[['y']] <- o[[mo]][[re]][[wtest[2]]][!is.na(o[[mo]][[re]][[wtest[2]]])]
                      if (legend) {
                        dot[['auc']] <- c(round(mean(unlist(auc[[mo]][[re]][[wtest[1]]]),na.rm=TRUE),3),round(mean(unlist(auc[[mo]][[re]][[wtest[2]]]),na.rm=TRUE),3))
                        dot[['wtest']] <- wtest
                      }
                    } else {
                      if (legend) {
                        dot[['auc']] <-round(mean(unlist(auc[[mo]][[re]][[wtest[1]]]),na.rm=TRUE),3)
                        dot[['wtest']] <- wtest
                      }
                    }
                    
                    if ('cex' %in% ndot) dot[['cex']] <- dot[['cex']] - np*z
                    else dot[['cex']] <- 1 - np*z
                    dot[['smooth']] <- smooth
                    
                    if (!is.null(Main) && np > 1) dot[['main']] <- mn[[mo]][[re]]
                    
                    do.call(.plot.roc,dot)
                  }
                  
                }
              }
              
            }
            par(op)
          }
)

#---------



if (!isGeneric("getRoc")) {
  setGeneric("getRoc", function(x,p,...)
    standardGeneric("getRoc"))
}  

setMethod('getRoc', signature(x='vector',p='vector'),
          function(x, p) {
            x <- .occurrence(x)
            w <- which(!is.na(x) & !is.na(p))
            p <- p[w]; x <- x[w]
            if (length(x) != length(p)) stop('observed and predicted vectors should have the same length')
            if (!.isBinomial(x) || !all(sort(unique(x)) == c(0,1))) stop('observed data should have 1 and 0 records (i.e., presence-absence)...')
            .roc(x,p)
          }
)


setMethod('getRoc', signature(x='sdmModels'),
          function(x, p,...) {
            if (missing(p)) p <- NULL
            wtest <- NULL
            a <- c('species','method','replication','run','w','test')
            dot <- list(...)
            ndot <- names(dot)
            
            mi <- .getModel.info(x)
            u <- unique(mi[,2])
            m <- unique(mi[,3])
            r <- unique(mi[,4])
            
            if (length(ndot) > 0) {
              ndot <- .pmatch(ndot,a)
              w <- !is.na(ndot)
              ndot <- ndot[w]
              dot <- dot[w]
              names(dot) <- ndot
              if (!'w' %in% ndot) {
                if (!is.null(p) && is.numeric(p) && length(p) == 1) dot[['w']] <- p
              }
              if ('test' %in% ndot) {
                wtest <- dot[['test']]
                dot <- dot[-which(ndot == 'test')]
                if (is.character(wtest))  wtest <- .pmatch(wtest,c("training","dep","indep"))
                else if (is.numeric(wtest) && any(wtest %in% 1:3)) wtest <- c("training","dep","indep")[wtest[wtest %in% 1:3]]
                
                wtest <- unlist(lapply(wtest,function(x) {
                  if (x == 'dep') 'test.dep'
                  else if (x == 'indep') 'test.indep'
                  else 'training'
                }))
              } else {
                wtest <- names(x@models[[1]][[1]][[1]]@evaluation)
              }
              
              if ('species' %in% ndot) {
                if (length(dot[['species']]) > 1) {
                  dot[['species']] <- dot[['species']][1]
                  warning('only the first specified species is considered!')
                }
                if (is.numeric(dot[['species']])) {
                  if (length(u) <= dot[['species']]) dot[['species']] <- u[dot[['species']]]
                  else stop('The specified species is not recognised!')
                } else {
                  dot[['species']] <- .pmatch(dot[['species']],u)
                  if (is.na(dot[['species']])) stop('The specified species is not recognised!')
                }
              } else {
                if (length(u) > 1) stop('This object contains models for several species, the names of species should be specified!')
              }
            } else if (!is.null(p) && is.numeric(p) && length(p) == 1) {
              dot[['w']] <- p
              wtest <- names(x@models[[1]][[1]][[1]]@evaluation)
            } else {
              if (length(u) > 1 | length(m) > 1 | length(r) > 1 ) stop('This object contains several models; the model ID (w) or the info including species name, model name, and replication method should be specified!')
              else {
                dot[['method']] <- m
                if (!is.na(r)) {
                  r <- .replicate.methodFix(r)
                  dot[['replication']] <- r
                }
              }
            }
            if (is.null(wtest)) stop('there is no evaluation in the model object!')
            else {
              wtest <- wtest[wtest %in% names(x@models[[1]][[1]][[1]]@evaluation)]
              if (length(wtest) > 1) {
                if ('test.indep' %in% wtest) wtest <- 'test.indep'
                else if ('test.dep' %in% wtest) wtest <- 'test.dep'
                else wtest <- 'training'
              }
            }
            
            dot[['x']] <- x
            mi <- do.call(.getModel.info,dot)
            o <- list()
            sp <- as.character(mi[1,2])
            m <- as.character(mi[1,3])
            for (i in seq_along(mi[,1])) {
              o[[i]] <- .roc(x@models[[sp]][[m]][[as.character(mi[i,1])]]@evaluation[[wtest]]@observed,x@models[[sp]][[m]][[as.character(mi[i,1])]]@evaluation[[wtest]]@predicted)
            }
            o
          }
)


