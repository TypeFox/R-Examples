Rayfan <-
function(phaselist, h, model, deltalist = 1:17 * 20, minp = .5, plot = TRUE, add = TRUE, col = rep('black',length(phaselist)), verbose = FALSE, mirror = FALSE){
  if(.Device == 'null device'){Earthplot(model)}
  
  imodel = ImproveModel(model)$newmodel
  rp = imodel$rp
  if(!add){plot(NaN,type='n',xlim = 1.15 * c(-rp,rp),ylim = 1.15 * c(-rp,rp),ann=FALSE,axes=FALSE)}
  outp = NULL
  outdist = NULL
  outphase = NULL
  outt = NULL
  for(i in 1:length(phaselist)){
    phase = phaselist[i]
    if(verbose){print(paste(phase,': ',i,' of ',length(phaselist),sep=''))}
    
    if(mirror){
      ttinfo = Traveltime(phase,c(deltalist,360-deltalist),h,imodel)
    } else {
      ttinfo = Traveltime(phase,deltalist,h,imodel)
    }
    p = ttinfo$p
    p = p[!is.na(p) & p >= minp]
    x = FindDist4p(phase,h,imodel,takeoff = ttinfo$angles)    
    w = which(!is.na(x$dist) & is.finite(x$dist))
    for(j in w){
      if(round(x$dist[j],2) %in% round(deltalist,2)) {
        PolarPlot(x$segx[[j]], rp-x$segz[[j]], method=lines, degrees=TRUE, geographical=TRUE, col = col[i])
      }
      if(mirror & (360 - round(x$dist[j],2)) %in% round(deltalist,2)){
        PolarPlot(-x$segx[[j]], rp-x$segz[[j]], method=lines, degrees=TRUE, geographical=TRUE, col = col[i])
      }
    }
    outp = c(outp, p)
    outdist = c(outdist, x$dist[w])
    outphase = c(outphase, rep(phase, length(p)))
    outt = c(outt, ttinfo$t[w])
  }
  
  return(list(p = outp, t = outt, dist = outdist, phase = outphase))
}

