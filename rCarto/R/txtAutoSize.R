txtAutoSize <-
function(txtCex,height){
  if(is.null(txtCex)==TRUE){
    txtCex<-height*1.2/10}
  return(txtCex)
}
