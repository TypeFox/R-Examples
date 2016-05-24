showVarExplained <- function(result, col=c('grey20','grey43','grey65')){
  
  old.par <- par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(par(old.par))
  
  l <- length(result$data)
  VarJoint = rep(0,l)
  for(i in 1:l) VarJoint[i] = norm(result$joint[[i]],type='F')^2/norm(result$data[[i]], type = 'F')^2
  VarIndiv = rep(0,l)
  for(i in 1:l) VarIndiv[i] = norm(result$individual[[i]],type='F')^2/norm(result$data[[i]], type = 'F')^2
  VarResid = 1-VarJoint-VarIndiv
  
  par(mar=c(5.1,4.1,4.1,0))
  layout(matrix(c(1,2),1,2),heights=c(5,5),widths=c(5,2))
  barplot(rbind(VarJoint,VarIndiv,VarResid),col = col,main = "Variation Explained",names.arg=names(result$data))
  par(mar=c(0,0,0,0))
  plot.new()
  legend(x=0.05,y=0.8,legend=c('Joint','Individual','Residual'), bty = "n",fill= col)
}