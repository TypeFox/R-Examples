kparts <-
function(x, y, parts, maxiter = 50, trials = 3, nblind=FALSE,
                  trialprint = TRUE, iterprint = FALSE) {
  
  G = data.frame('g' = tapply(x, x, min)
                 ,'n' = tapply(y, x, length)
                 ,'x' = tapply(y, x, mean)
                 ,'c' = as.integer(NA)
                 ,'cdn' = as.integer(NA)
                 ,'cup' = as.integer(NA)
                 ,'idn' = as.integer(NA)
                 ,'iup' = as.integer(NA)
                 ,'sse' = as.numeric(NA)
                 ,'ssedn' = as.numeric(NA)
                 ,'sseup' = as.numeric(NA)
                 ,'movup' = as.integer(NA)
                 ,'movdn' = as.integer(NA))
  
  if (nblind==TRUE) {
    G$n = rep(1, nrow(G))
  }
  
  rownames(G) = seq.int(1, length(G$g))
  for(r in 1:trials) {
    repeat{
      m = sort(runif(parts ,min(x), max(x)))
      G$c = apply(abs(outer(G$g, m, "-")), 1, which.min)
      if (any(duplicated(m))==FALSE & length(m)==length(unique(G$c))) {break}
    }
    
    l = length(m)
    
    M = data.frame('c' = seq.int(1:l)
                   ,'x' = as.numeric(NA)
                   ,'csi' = as.integer(NA)
                   ,'cei' = as.integer(NA)
                   ,'idn' = as.integer(NA)
                   ,'iup' = as.integer(NA)
                   ,'sse' = as.numeric(NA)
                   ,'ssedn' = as.numeric(NA)
                   ,'sseup' = as.numeric(NA)
                   ,'nsi' = as.numeric(NA)
                   ,'nei' = as.numeric(NA)
                   ,'isJump' = as.integer(NA))
    
    M = rbind(M, c(max(M$c)+1, rep.int(NA, (length(M[1,]) - 1))))
    M = rbind(M, c(0, rep.int(NA, (length(M[1,]) - 1))))
    
    M$csi[1:l] = tapply(kp.rint(G), G$c, min)
    M$cei[1:l] = tapply(kp.rint(G), G$c, max)
    
    for (k in 1:maxiter) {
      mr = kp.mre(G, M)
      M = mr$M
      l = mr$l
      G$c = mr$c; G$cdn = mr$cdn; G$cup = mr$cup
      G$sse = mr$sse; G$ssedn = mr$ssedn; G$sseup = mr$sseup
      G$idn = mr$idn; G$iup = mr$iup
      
      if (exists('J')==FALSE) {
        J = sum(G$sse)
        if (trialprint==1) {
          print(paste('initial',0,J))
        }
      }
      
      # old code that actually worked
       cJump = max(which(M$sse[1:l]==min(M$sse[1:l])))
       
       for (i in M$csi[cJump]:M$cei[cJump]) {
         if ((i==M$csi[cJump] & i!=length(G$g)) | M$cei[cJump] == 1) {
           sseJump = sum(G$sseup[i:G$iup[i]])
         } else if (M$csi[cJump] == length(G$g)) {
           sseJump = sum(G$ssedn[G$idn[i]:i])
         } else if (i==M$cei[cJump] & i!=1 ) {
           sseJump = c(sseJump,(sum(G$ssedn[G$idn[i]:i])))
         } else {
           sseJump = c(sseJump,
                       sum(G$sseup[i:G$iup[i]]) 
                       + sum(G$ssedn[G$idn[(i-1)]:(i-1)]))
         }
       }
       
      
      if (M$csi[cJump] == length(G$g)) {
        #  8/19/2013 code change
        #   rowJump = M$csi[cJump] + which.min(sseJump)
        rowJump = M$csi[cJump]
      } else {
        rowJump = M$csi[cJump] + which.min(sseJump) - 1
      }
      
      if (cJump != l & M$cei[cJump] >= rowJump) {
        M$csi[(cJump+1)] = rowJump
      }
      
      if (cJump != 1 & cJump != l & M$csi[cJump] <= rowJump - 1) {
        M$cei[(cJump-1)] = rowJump - 1
      } else if (cJump!=1 & M$csi[cJump]<= rowJump - 1) {
        M$cei[(cJump-1)] = rowJump
      } else if (M$cei[cJump] == nrow(G)) {
        M$cei[cJump-1] = nrow(G)
      }
      
      M[cJump,] = c(max(M$c) + 1
                    ,rep.int(NA,length(M[i,]) - 2)
                    ,1)
      
      mr = kp.mre(G,M)
      M = mr$M
      l = mr$l
      G$c = mr$c; G$cdn = mr$cdn; G$cup = mr$cup
      G$sse = mr$sse; G$ssedn = mr$ssedn; G$sseup = mr$sseup
      G$idn = mr$idn; G$iup = mr$iup
      
      cJump = which(M$isJump==1)
      
      cMost = max(which(M$sse[1:l]==max(M$sse[kp.rint(M) <= l 
                                              & (M$cei - M$csi - 2) > 0])))
      
      for (i in (M$csi[cMost] + 1):M$cei[cMost]) {
        iu = G$iup[i] #max i upper cluster
        id = G$idn[(i-1)] # min i of lower cluster
        ii = i - 1 # max i of lower cluster
        xu = sum(G$x[i:iu]*G$n[i:iu])/sum(G$n[i:iu])
        xd = sum(G$x[id:ii]*G$n[id:ii])/sum(G$n[id:ii])
        
        if (i != M$csi[cMost]+1) {
          sseMost = c(sseMost,
                      sum((G$x[i:iu] - xu)^2 * G$n[i:iu])
                      + sum((G$x[id:ii] - xd)^2 * G$n[id:ii]))
        } else {
          sseMost = sum((G$x[i:iu] - xu)^2 * G$n[i:iu]) 
          + sum((G$x[id:ii] - xd)^2 * G$n[id:ii])
        }
      }
      
      cJumpNewsi = M$csi[cMost] + which.min(sseMost)
      
      M$csi[cJump] = cJumpNewsi
      M$cei[cJump] = M$cei[cMost]
      M$cei[cMost] = cJumpNewsi - 1
      
      M$isJump[cJump]=NA
      
      mr = kp.mre(G,M)
      M = mr$M
      l = mr$l
      G$c = mr$c; G$cdn = mr$cdn; G$cup = mr$cup
      G$sse = mr$sse; G$ssedn = mr$ssedn; G$sseup = mr$sseup
      G$idn = mr$idn; G$iup = mr$iup
      
      for (j in 1:maxiter) {
        
        for (i in 1:length(G$c)) {
          G$movdn[i] = ifelse(ifelse(is.na(G$ssedn[i])==TRUE
                                     ,Inf
                                     ,sum(G$ssedn[G$idn[i]:i]))<
                                sum(G$sse[G$idn[i]:i]),1,0)
          
          G$movup[i] = ifelse(ifelse(is.na(G$sseup[i])==TRUE
                                     ,Inf
                                     ,sum(G$sseup[i:G$iup[i]]))<
                                sum(G$sse[i:G$iup[i]]),1,0)
        }
        
        for (i in 1:l) {
          mniup = tapply(kp.rint(G[G$movup==1 & G$c==i,])
                         ,G$c[G$movup==1 & G$c==i],min)
          mxidn = tapply(kp.rint(G[G$movdn==1 & G$c==i,])
                         ,G$c[G$movdn==1 & G$c==i],max)
          
          if (any(mxidn)==TRUE) {
            M$idn[i] = mxidn
            M$ssedn[i] = sum(G$ssedn[G$idn[mxidn]:mxidn])
          }
          if (any(mniup)==TRUE) {
            M$iup[i] = mniup
            M$sseup[i] = sum(G$sseup[mniup:G$iup[mniup]])
          }
          
        }
        
        for (i in 1:l) {   
          
          if (is.na(M$idn[i])==FALSE 
             & is.na(M$iup[i])==FALSE 
             & M$idn[i]>=M$iup[i]) {
            if (M$sseup[i]<=M$ssedn[i]) {
              if (M$iup[i]==M$csi[i]) {
                M$idn[i] = NA
                M$ssedn[i] = NA
              } else {
                nBest = M$iup[i] - 1
                mxidn = tapply(kp.rint(G[G$movdn==1 & G$c==i & kp.rint(G)<=nBest,])
                               ,G$c[G$movdn==1 & G$c==i & kp.rint(G)<=nBest],max)
                if (any(mxidn)==TRUE) {
                  M$idn[i] = mxidn
                  M$ssedn[i] = sum(G$ssedn[G$idn[mxidn]:mxidn])
                } else {
                  M$idn[i] = NA
                  M$ssedn[i] = NA
                }
              }
            } else {
              if (M$idn[i]==M$cei[i]) {
                M$iup[i] = NA
                M$sseup[i] = NA
              } else {
                nBest = M$idn[i] + 1
                mniup = tapply(kp.rint(G[G$movup==1 & G$c==i & kp.rint(G)>=nBest,])
                               ,G$c[G$movup==1 & G$c==i & kp.rint(G)>=nBest],max)
                
                if (any(mniup)==TRUE) {
                  M$iup[i] = mniup
                  M$sseup[i] = sum(G$sseup[mniup:G$iup[mniup]])
                } else {
                  M$iup[i] = NA
                  M$sseup[i] = NA
                }
              }
            }
          }    
        }
        
        for (i in 1:l) {
          if (is.na(M$iup[i]+M$idn[i+1])==FALSE) {
            if (M$sseup[i]>=M$ssedn[i+1]) {
              M$sseup[i] = NA
              M$iup[i] = NA
            } else {
              M$ssedn[i+1] = NA
              M$idn[i+1] = NA
            }
          }
        }
        
        for (i in 1:l) {
          if (is.na(M$idn[i])==FALSE) {
            M$nei[i-1] = M$idn[i]
            M$nsi[i] = M$idn[i] + 1
          }
          if (is.na(M$iup[i])==FALSE) {
            M$nsi[i+1] = M$iup[i]
            M$nei[i] = M$iup[i]-1
          }
        }
        
        M$nsi = ifelse(is.na(M$nsi)==TRUE,M$csi,M$nsi)
        M$nei = ifelse(is.na(M$nei)==TRUE,M$cei,M$nei)
        
        for (i in 1:l) {
          if (M$nsi[i]>=M$nei[i]) {
            M[i,] = c(max(M$c)+1
                      ,rep.int(NA,length(M[i,])-2)
                      ,1)
          }
        }
        
        if (sum(ifelse(is.na(M$isJump)==TRUE,0,1))!=0) {
          M = M[order(M$csi,-is.na(M$isJump),M$c),]
          l = length(M$c) - (sum(ifelse(is.na(M$isJump)==TRUE,0,1))+2)
          rownames(M) = seq.int(1,length(M$c))
          
          for (i in 1:l) {
            if (i==1) {
              M$nsi[i]=1
            } else if (i==l) {
              M$nsi[i] = M$nei[(i-1)] + 1
              M$nei[i] = length(G$c)
            } else {
              M$nsi[i] = M$nei[(i-1)] + 1
            }
          }
          
          M$csi = M$nsi
          M$cei = M$nei
          
          mr = kp.mre(G,M)
          M = mr$M
          l = mr$l
          G$c = mr$c; G$cdn = mr$cdn; G$cup = mr$cup
          G$sse = mr$sse; G$ssedn = mr$ssedn; G$sseup = mr$sseup
          G$idn = mr$idn; G$iup = mr$iup
          
          for (i in 1:length(M$isJump[is.na(M$isJump)==FALSE])) {
            cMost = max(which(M$sse[1:l]==max(M$sse[kp.rint(M)<=l 
                                                    & (M$cei - M$csi - 2) > 0])))
            cJump = kp.rint(M[M$isJump==1 & is.na(M$isJump)==FALSE,])
            for (i in (M$csi[cMost]+1):M$cei[cMost]) {
              iu = G$iup[i] #max i upper cluster
              id = G$idn[(i-1)] # min i of lower cluster
              ii = i - 1 # max i of lower cluster
              xu = sum(G$x[i:iu]*G$n[i:iu])/sum(G$n[i:iu])
              xd = sum(G$x[id:ii]*G$n[id:ii])/sum(G$n[id:ii])
              
              if (i != M$csi[cMost]+1) {
                sseMost = c(sseMost,
                            sum((G$x[i:iu] - xu)^2 * G$n[i:iu])
                            + sum((G$x[id:ii] - xd)^2 * G$n[id:ii]))
              } else {
                sseMost = sum((G$x[i:iu] - xu)^2 * G$n[i:iu]) + 
                  sum((G$x[id:ii] - xd)^2 * G$n[id:ii])
              }
            }
            cJumpNewsi = M$csi[cMost] + which.min(sseMost)
            
            
            M$csi[cJump] = cJumpNewsi
            M$cei[cJump] = M$cei[cMost]
            M$cei[cMost] = cJumpNewsi - 1
            M$isJump[cJump]=NA
            mr = kp.mre(G,M)
            M = mr$M
            l = mr$l
            G$c = mr$c; G$cdn = mr$cdn; G$cup = mr$cup
            G$sse = mr$sse; G$ssedn = mr$ssedn; G$sseup = mr$sseup
            G$idn = mr$idn; G$iup = mr$iup
          }
        } else {
          M$csi = M$nsi
          M$cei = M$nei
        }
        
        mr = kp.mre(G,M)
        M = mr$M
        l = mr$l
        G$c = mr$c; G$cdn = mr$cdn; G$cup = mr$cup
        G$sse = mr$sse; G$ssedn = mr$ssedn; G$sseup = mr$sseup
        G$idn = mr$idn; G$iup = mr$iup
        J = c(J,sum(G$sse))
        if (iterprint==1) {
          print(paste('    iteration',j,J[length(J)]))
        }
        
        if (j>=2) {
          if (J[length(J)]==J[(length(J)-1)]) {
            break
          }
        }
      }
      J = c(J,sum(G$sse))
      if (iterprint==1) {
        print(paste(' Partition Jump',k,J[length(J)]))
      }
      
      if (k>=2) {
        if (J[length(J)]==J[(length(J)-1)]) {
          break
        }
      }
    }
    
    if (trialprint==TRUE) {
      print(paste('trial',r,J[length(J)]))
    }
    
    
    if (exists('best.J')==TRUE & 
          J[length(J)] < ifelse(exists('best.J')==FALSE,Inf,best.J)) {
      best.J = J[length(J)]
      
      best.M = data.frame('parts' = M$c[1:l],
                          'range' = paste(tapply(G$g,G$c,min),
                                          tapply(G$g,G$c,max),
                                          sep='-'))
      
      best.G = merge(data.frame('x' = G$g
                                ,'y' = G$x
                                ,'parts' = G$c),
                     best.M,by.x = 'parts',by.y = 'parts')
    } else if (exists('best.J')==FALSE) {
      best.J = J[length(J)]
      
      best.M = data.frame('parts' = M$c[1:l],
                          'range' = paste(tapply(G$g,G$c,min),
                                          tapply(G$g,G$c,max),
                                          sep='-'))
      
      best.G = merge(data.frame('x' = G$g
                                ,'y' = G$x
                                ,'parts' = G$c),
                     best.M,by.x = 'parts',by.y = 'parts')
    }
  }
  structure(list(partitions = best.M, data = best.G, ss = best.J),class = 'kparts')
}
