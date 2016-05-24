"plotResids" <- function (multimodel,  multitheta,  plotoptions)
{
  if(dev.cur() != 1)
    dev.new()
  plotrow<-2
  plotcol<-2
  par(plotoptions@paropt)
  par(mfrow = c(plotrow, plotcol))
  par(mar=c(5.1,7.1,4.1,2.),cex.lab=1.6, cex.axis=1.6,cex.main=1.4, 
      mgp=c(4,1, 0))
  ## RESIDUALS 
  ## make a list of resid matrix
  m <- multimodel@modellist
  res <- multimodel@fit@resultlist
  model <- multimodel@modellist[[1]]
  increasing_x2 <- model@x2[2] > model@x2[1]
  residlist <- svdresidlist <- list() 
  for(i in 1:length(m)) {
    residuals <- matrix(nrow = m[[i]]@nt, ncol = m[[i]]@nl)
    for(j in 1:length(res[[i]]@resid)){ 
      if(m[[i]]@clpType == "x2")
        residuals[,j] <- res[[i]]@resid[[j]]
      else
        residuals[j,] <- res[[i]]@resid[[j]]
    }
    svdresidlist[[length(svdresidlist)+1]] <- doSVD(residuals,2,2) 
    residlist[[length(residlist)+1]] <- residuals 
  }
  ## matplot function with "log" option is not compatible with 
  ## neg. x values; do the below to avoid warning
  xpos <- model@x
  xpos[which(model@x<=0)]<-NA	
  if(increasing_x2) {
    limd<- max(  max(residlist[[1]]), abs(min(residlist[[1]]))) 
    if (! (any(diff(model@x) <= 0) || any(diff(model@x2) <= 0)))
      image(model@x, model@x2, residlist[[1]], 
            xlab = plotoptions@xlab, ylab = plotoptions@ylab,
            main = "Residuals Dataset 1", zlim=c(-limd,limd),
            col = diverge_hcl(40, h = c(0, 120), c = 60, l = c(45, 90),
              power = 1.2))
  }
  if(model@nt > 1 && model@nl > 1){ 
    matplot(xpos, svdresidlist[[1]]$left[,1], type = "l",
            main = "1st left sing. vec. residuals ", log = "x", xlab =
            plotoptions@xlab, col = 1, ylab="")
    if(length(m) > 1){
      for(i in 2:length(m)) {
        matlines(m[[i]]@x, svdresidlist[[i]]$left[,1],log="x", type ="l",
                 col = i)
      }      
    }
    matplot(model@x2, svdresidlist[[1]]$right[1,], type = "l",  
            main = "1st right sing. vec. residuals ", xlab = plotoptions@ylab, 
            col = 1, ylab="")
    if(length(m) > 1){
      for(i in 2:length(m)) {
        matlines(m[[i]]@x2, svdresidlist[[i]]$right[1,], type ="l", col = i)
      }
    }
    plot(1:length(svdresidlist[[1]]$values),
         log10(svdresidlist[[1]]$values), 
         main = "Sing. values residuals", type= "b",xlab="",ylab="")
    if(length(m) > 1){
      for(i in 2:length(m)) {
        lines(1:length(svdresidlist[[i]]$values),
              log10(svdresidlist[[i]]$values), xlab="",ylab="", 
              type = "b", col = i)
      }
    }
  }
  ## MAKE PS
  if(dev.interactive() && length(plotoptions@makeps) != 0) {
    if(plotoptions@output == "pdf")
      pdev <- pdf 
    else  pdev <- postscript
    dev.print(device=pdev, 
              file=paste(plotoptions@makeps, "_resids.", 
                plotoptions@output, sep=""))
  }
  
}
