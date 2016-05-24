piFunction <-
function(fit, kparam, qparam){

  if(fit < qparam){
    newp = qparam + (1-qparam) * 
      sign((qparam - fit) / qparam)*abs((qparam - fit) / qparam)^(1/kparam)
  }else if(fit > qparam){
    newp = qparam + qparam * 
      sign((qparam - fit) / (1 - qparam))*abs((qparam - fit) / (1 - qparam))^(1/kparam)
  }else{
    newp = qparam
  }
  return(newp)
}
