plot.meta4diag = function(x, var.type="var1", add=FALSE, overlay.prior = TRUE, save = FALSE, width=5, height=5, ...){
  if(is.logical(save)){
    if(save){
      mainDir <- getwd()
      subDir <- "meta4diag_Plot"
      if (file.exists(subDir)){
        file.name = paste(subDir,"/",var.type,".pdf",sep="")
      } else {
        dir.create(file.path(mainDir, subDir))
        file.name = paste(subDir,"/",var.type,".pdf",sep="")
      }
      save.flag = TRUE
    }else{save.flag = FALSE}
  } else if(is.character(save)){
    name_temp = unlist(strsplit(basename(save), "[.]"))
    fileform = name_temp[length(name_temp)]
    if(fileform %in% c("pdf","eps","jpg","png")){
      save.flag = TRUE
      file.name = paste(subDir,"/",save,sep="")
    }else{
      save.flag = FALSE
      stop("Please give the correct file name!")
    }
  } else{stop("Argument \"save\" could be FALSE, TRUE or a file name.")}
  
  fixed.name = rownames(x$summary.fixed)
  hyper.name = c("var1", "var2", "rho")
  fullnames = c(fixed.name, hyper.name)
  if(!(var.type %in% fullnames)){
    stop(paste("Please give the correct \"type\" name, which should be ", paste(fullnames,collapse=", "),sep=""))
  }
  fixed.prior = data.frame(x=seq(-10,10,len=100),y=dnorm(seq(-10,10,len=100),mean=0,sd=sqrt(1000)))
  fullmarginals = append(x$marginals.fixed,x$marginals.hyperpar)
  ind = which(fullnames==var.type)
  if(var.type=="var1" || var.type=="var2"){
    marginals.plot = fullmarginals[[ind]]
  }else{
    marginals.plot = INLA::inla.smarginal(fullmarginals[[ind]])
  }
  
  
  xlabnames = c(fixed.name, rownames(x$summary.hyperpar))
  
  if(add){
    lines(marginals.plot, ...)
  }else{
    if(save.flag){
      if(fileform=="eps"){
        setEPS()
        postscript(file.name, width=width, height=height,...)
      }else if(fileform=="pdf"){
        pdf(file.name, width=width, height=height,...)
      }else if(fileform=="jpg"){
        jpeg(filename = file.name,
             width = width, height = height, units = "in")
      }else{
        png(filename = file.name,
            width = width, height = height, units = "in")
      }
    }
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    plot(marginals.plot, type="l", xlab=xlabnames[ind], ylab="",xaxs = "r",family="sans",xaxt="s",yaxt="s",bty="o",...)
    if(!x$misc$wishart.flag){
      if(overlay.prior){
        nom = length(fullnames)
        if(ind<=(nom-3)){
          lines(fixed.prior,lty=2,col="darkgray",...)
        }else if(ind==(nom-2)){
          lines(x$priors.density[[1]],lty=2,col="darkgray",...)
        }else if(ind==(nom-1)){
          lines(x$priors.density[[2]],lty=2,col="darkgray",...)
        }else if(ind==nom){
          lines(x$priors.density[[3]],lty=2,col="darkgray",...)
        }else{
          stop("Wrong var.type!")
        }
      }
    }
    if(save.flag){
      dev.off()
    }
  }
}