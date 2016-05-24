kp.mre <-
function(G,M) {
  M = M[order(M$csi,-is.na(M$isJump),M$c),]
  l = length(M$c) - (sum(ifelse(is.na(M$isJump)==TRUE,0,1))+2)
  rownames(M) = seq.int(1,length(M$c))
  M$c = seq.int(1,length(M$c))
  M$c[l+2] = 0
  M$isJump[!is.na(M$isJump)] = seq.int(1,length(M$isJump[!is.na(M$isJump)]))
  M$idn = NA; M$iup = NA; M$ssedn = NA; 
  M$sseup = NA; M$nsi = NA; M$nei = NA
  if (sum(is.na(M$isJump)) < length(M$c)) {
    M[is.na(M$isJump)==FALSE,][,2:(length(M[1,]) - 1)] = NA 
  }
  
  
  G$c = which(outer(kp.rint(G),M$csi[1:l],'>=') * 
                outer(kp.rint(G),M$cei[1:l],'<=')==1,2)[,2]
  
  G$cdn = G$c-1
  G$cup = G$c+1
  
  M$x[1:l] = tapply((G$x*G$n),G$c,sum)/tapply(G$n,G$c,sum)
  
  G$sse = ((G$x - merge(G,M,by.x = 'c',by.y ='c')$x.y)^2)*G$n
  G$ssedn = ((G$x - merge(G,M,by.x = 'cdn',by.y ='c')$x.y)^2)*G$n
  G$sseup = ((G$x - merge(G,M,by.x = 'cup',by.y ='c')$x.y)^2)*G$n
  
  G[,7:8] = merge(data.frame('c' = sort(unique(G$c))
                             ,'mn' = tapply(kp.rint(G),G$c,min)
                             ,'mx' = tapply(kp.rint(G),G$c,max))
                  ,G,by.x = 'c',by.y='c')[,2:3]
  
  M$sse[1:l] = tapply(G$sse,G$c,sum)
  
  structure(list(M = M,l=l,c=G$c,cdn=G$cdn,cup=G$cup,sse=G$sse,
                 ssedn=G$ssedn,sseup=G$sseup,iup=G$iup,idn=G$idn))
}
