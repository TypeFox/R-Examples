pbaAgain <-
function(previous, newx, k=NA){

  if(nrow(newx) > 1) stop('newx has more than one row')

  x = rbind(previous$input$x, previous$input$newx)
  tr = c(previous$input$tr, previous$results$newtr)

  if(is.na(k)){
    pbaResults = pba(x=x, tr=tr, newx=newx, k=previous$input$k, 
                     global=previous$input$global)
  }else{
    if(k == previous$input$k){
      pbaResults = pba(x=x, tr=tr, newx=newx, k=previous$input$k, 
                       global=previous$input$global)
    }else{
      ## warn if k is changing from previous allocation
      warning(paste('Balancing parameter k changed from',
                    round(previous$input$k,5),'to',round(k,5), sep=' '))
      pbaResults = pba(x=x, tr=tr, newx=newx, k=k, global=previous$input$global)
    }
  }
  return(pbaResults)
}
