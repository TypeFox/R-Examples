DistSummary <-
function(delta,h,model,phaselist = 'default', prop='vp',image.col=heat.colors(500),n=200,...){
  if(phaselist == 'default'){
    phaselist = c(
      'PS',
      'SP',
      'P',
      'PP',
      'PPP',
      'PKP',
      'PKPPKP','PKIKP','PKIKPPKIKP','PKKP',
      'PcP',
      'PcPPcP',
      'PcS',
      'ScP',
      'PcS2',
      'ScP2',
      'PKiKP','S','SS','SSS','ScS','ScSScS','SKS','SKSSKS','SKKS','SKIKS'
      )
  }
  
  Earthplot(model,prop,image.col,n)
  out=Rayfan(phaselist,h,model,deltalist=delta, verbose=TRUE,mirror=TRUE,...)
  return(out)
}

