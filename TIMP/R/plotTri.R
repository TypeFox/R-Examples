"plotTri" <- function(xx) {
  if(dev.cur() != 1)
    dev.new()
  par(mfrow=c(3, 1))	
  svddatalist <-doSVD(xx, 2, 2)
  
  matplot(svddatalist$left, type = "l", 
          main = "Left sing. vectors alpha",  
          log = "x")
  
  matplot(t(svddatalist$right), type = "l", 
          main = "Right sing. vectors alpha")
  
  plot(1:length(svddatalist$values), 
       log10(svddatalist$values), xlab="",
       main = "Sing. values alpha", type = "b")

}
