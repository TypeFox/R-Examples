#' @export
"simmodel"<-function(vpcpath,
                     nsim=200,
                     pstrat=NULL,
                     setseed=NULL,
                     pred.corr=NULL,
                     var.corr=FALSE,
                     pi=c(0.025, 0.5, 0.975),
                     pi.ci=c(0.025, 0.975),
                     bin.option=NULL,
                     bin.bound=NULL,
                     bin.center=NULL,
                     clean=FALSE,
                     hold=FALSE, 
                     ivar="t",
                     model.file="test.mdl", 
                     cols.file="cols1.txt", 
                     data="data1.txt")
{ rootwd<-getwd()
  
  if(is.null(pstrat) | (!is.null(pstrat) & length(pstrat)<4)){
    if(any(file.exists(paste(vpcpath,c("mpiNLME7.exe","NLME7.exe"),sep="/"))) &
         all(file.exists(paste(vpcpath,c(model.file, data, cols.file),sep="/")))){
      
      ######set seed########################################
      if(!is.null(setseed)){
          setseed = as.integer(setseed)
        if(is.integer(setseed) & setseed>0 & setseed<=32767){
          set.seed(setseed) # 1-32767
        }else{
          write("Error: Seed is an integer number 1-32767. Random seed assigned.","")
          setseed<-NULL
        }
      }
      if(is.null(setseed)){
        setseed<-floor(runif(1,1,32767))
      }
      
      ####prediction correction options######################
      if(!is.null(pred.corr)){
        if(grepl("prop", pred.corr,ignore.case=T)) predpc=ifelse(var.corr," /predvpc"," /predpc") else{
          if(grepl("add", pred.corr,ignore.case=T)) 
            predpc=ifelse(var.corr, " /predvpc /predpcadd"," /predpc /predpcadd") else
              {write("Error: Incorrect pred.corr option. Default option pred.corr=NULL is used instead.","")
               pred.corr=NULL}
        }
      }
      if(is.null(pred.corr)) predpc=""
      
      
      ####binning options######################
      if(!is.null(bin.option)){
        if(grepl("k",bin.option, ignore.case=T)){bin.option=" /predkmeans"}else{
          if(grepl("bound",bin.option,ignore.case=T)&!is.null(bin.bound)){
            bin.option=paste(" /predboundaries  \" ", paste(bin.bound,collapse=","), "\" " ,sep=" ")}else
            if(grepl("cent",bin.option,ignore.case=T)&!is.null(bin.center)){
              bin.option=paste(" /predcenters  \" ", paste(bin.center,collapse=","), "\" ",sep=" ")}else
              {write(paste("Error: Incorrect bin.option or bin.bound/bin.center option.", 
                     ":\nDefault option pred.corr=NULL is used instead", sep=" "),"")
               bin.option=NULL}
        }
      }
      if(is.null(bin.option)) bin.option=" /prednobin"
      
      
      setwd(vpcpath)
      write(paste("Model",vpcpath,":\nSimulation/VPC N =",nsim,
                  ifelse(is.null(pstrat),"",paste("stratified by",paste(pstrat,collapse=", ")))),"")
      cat(paste("Executing simulation ... "))
      flush.console()
            
      ###change pi.ci format
      pi=pi*100
      pi.ci=pi.ci*100
      
      ### Create predckargs.txt file
      if(is.null(pstrat)) pstrarg<-"" else 
      { write(pstrat,file="strata.txt")
        pstrarg<-paste("/pstrat",c(1:length(pstrat))," \"",pstrat,"\"",sep="")}
      simargs<-paste("/predn",nsim, "/predout predout.csv /pcseed",setseed, 
                     predpc, bin.option, " /predx ", ivar,
                     "/pcpi", paste(pi, collapse=","), " /pcpe", paste(pi.ci, collapse=","), pstrarg)
      nlmeargs<-paste("/m 6 /n 0 /o 6 /e 0",
                      "-xnp 0 -anagrad 0 -logtran 1 -xrestart 0 -xnorderagq 1 -xfocehess 1 -xverbose 1 -xstderr 0",
                      "-xlameth 1 -xlandig 7 -xlatol 0.01 -xblmeth 1 -xblndig 13 -xbltol 0.002 -sort -csv",
                      "cols1.txt data1.txt out0001.txt")
      
      write(c(simargs,nlmeargs),"predckargs.txt")
      write(ivar, "ivar.txt")
      
      shellrun<-ifelse(hold,"cmd.exe /q /k","cmd.exe /q /c")
      
      try(system(paste(shellrun,"call \"./NLME7.exe\" ",
                       simargs,nlmeargs,"1> err1.txt 2> err2.txt"),wait=TRUE,invisible=TRUE))
      
      # mpi option, currently not implemented.
      #ncpu<-as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))
      #if(!all(is.numeric(nodes) & as.integer(nodes)==nodes & nodes>0 & nodes<ncpu)) nodes<-ncpu
      #try(system(paste(shellrun,"call \"%PhoenixMPIDir%bin/mpiexec\" -n",nodes,"-localonly ./mpiNLME7.exe",
      #                 simargs,nlmeargs,"1> err1.txt 2> err2.txt"),wait=TRUE,invisible=TRUE))

      
      ### Delete executables
      if(clean){
        unlink(list.files(pattern="\\.exe$"))
      }
      
      setwd(rootwd)
      cat("Done!\n")
    }else{
      write("Error: Required files are missing. Execute phxnlme.R first.","")
    }
  }else{
    write("Error: A maximun of 3 stratification variables exceeded.","")
  }
}
