plot2D<-function(seerSet, write=TRUE,outDir="~/Results/plots",col="red") {
  with(seerSet, {
    #     if(!file.exists(fD)) stop(paste0(fD,"does not exist!")) 
    #     load(fD) # brings in dataframe D
    if(nchar(bfn)==0) {bfn="tmp"
    warning("no base file name, so files will go to subfolder tmp")}
    outDir=file.path(outDir,bfn)
    if(!file.exists(outDir))  {   print(paste("Creating directory",outDir))
      dir.create(outDir,recursive=TRUE)    }
    #     require(rgl)
    #     require(dplyr)
    secondS=levels(D$cancer)
    #     head(D,2)
    print("If rgl is being loaded for the first time in this R session, make the plot bigger and get a nice angle on it before hitting returns at the R command prompt.")
    for (i in secondS) {
      pd=D%>%filter(cancer==i)%>%select(year,age,incid,Eincid)
      #     head(pd)
      #     sapply(pd,class)
      (nyrs=length(yrs<-unique(pd$year)))
      (length(ages<-unique(pd$age))*nyrs)
      if (!i%in%c("SLL","CLL")) pd$incid[pd$incid==0]=0.002 # raise a notch from 0.001 to avoid plotting -3 if all above 0.002
       else pd$incid[pd$incid==0]=0.0002 
      with(pd,plot3d(age,year,log10(incid),xlab="",ylab="",zlab="")) #,sub=paste0(j," ",i,"s")))
      # with(pd,plot3d(age,year,log10(incid),xlab="Age",ylab="Year",zlab="Log Incidence",cex=2)) #,sub=paste0(j," ",i,"s")))
      with(pd,surface3d(x=ages,y=yrs,z=matrix(log10(Eincid),ncol=nyrs),alpha=0.8,col=col)) 
      # with(pd,surface3d(x=ages,y=yrs,z=matrix(log10(Eincid),ncol=nyrs),alpha=0.8,col=col,lit=FALSE)) #=>slightly smaller
      # with(pd,persp3d(x=ages,y=yrs,z=matrix(log10(Eincid),ncol=nyrs),alpha=0.8,col=col,lit=FALSE)) #also too big
      # with(pd,surface3d(x=ages,y=yrs,z=matrix(log10(Eincid),ncol=nyrs),
                        # back="cull",front="cull",alpha=0.8,col=col)) # nothing seems to help=>use full size pngs
      clear3d(type="lights")
      light3d(theta = -90, phi = 75) 
      cat(i,"  ...  hit return at command prompt to move on:\n")
      readline() # way to pause until input
      if (write) {
        rgl.snapshot(filename=(f<-paste0(outDir,"/",i,".png")),fmt="png", top=TRUE)
        cat("writing file:",f,"\n")
        # title3d(paste0(i,sex,"s", col='black', line=3))
#         rgl.postscript(filename=(f<-paste0(outDir,"/",i,".eps")),fmt="eps") #these are big
#         cat("writing file:",f,"\n")
#         savePlot(filename=(f<-paste0(outDir,"/",i,".tiff")),type="tiff") #never worked
#         cat("writing file:",f,"\n")
        
      }
    } #i loop over cancers
  }) # end with 
  
}  # plot2D function
