#' @export
"bootmodel"<-function(model=NULL,nodes=NULL,
                      method=5,
                      niter=1000,
                      nboot=500,
                      bstrat=NULL, 
                      setseed=NULL,
                      clean=TRUE,
                      hold=FALSE){
  
  # Read model specification for fit
  mod.spec = read.csv("model.spec.csv")
  mdl=as.character(mod.spec$model.file)
  dat=as.character(mod.spec$data)
  cols=as.character(mod.spec$cols.file)
  
  if(!(method %in% c(2:6))){
    write(paste("Error: Wrong estmation method number. Valid method numbers:\n",
                "\t2 IT2S-EM (Iterated 2-stage expectation-maximization)\n",
                "\t3 FOCE L-B (First-Order Conditional Estimation, Lindstrom-Bates)\n",
                "\t4 FO (First Order)\n",
                "\t5 FOCE-ELS (Extended Least Squares)\n",
                "\t6 Naive pooled",
                sep=""),"")
  }else{
    if(nboot>0 & nboot<10000){
      if(is.null(bstrat) | (!is.null(bstrat) & length(bstrat)<4)){
        if(any(file.exists(c("mpiNLME7.exe","NLME7.exe"))) &
             all(file.exists(c(mdl,dat,cols)))){
          
          if(!is.null(setseed)){
              setseed = as.integer(setseed)            
            if(is.integer(setseed) & setseed>0 & setseed<=32767){
              set.seed(setseed) # 1-32767
            }else{
              write("Error: Seed is an integer number 1-32767. Random seed assigned.","")
            }
          }

          write(paste("Model",model,":\nBootstrap N =",nboot,
                      ifelse(is.null(bstrat),"",paste("stratified by",paste(bstrat,collapse=", ")))),"")
          write("Executing bootstrap ... ","")
          flush.console()
          
          ncpu<-as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))
          if(!all(is.numeric(nodes) & as.integer(nodes)==nodes & nodes>0 & nodes<ncpu)) nodes<-ncpu
          
          ### Create nlmeargs.txt file
          bstrarg<-ifelse(!is.null(bstrat),paste("/bstrat",c(1:length(bstrat))," \"",bstrat,"\"",sep=""),"")
          nlmeargs<-paste("/m",method,"/n",niter,"/o 6 /e 0",
                          "-xnp 0 -anagrad 0 -logtran 1 -xrestart 0 -xnorderagq 1 -xfocehess 1 -xverbose 1 -xstderr 0",
                          "-xlameth 1 -xlandig 7 -xlatol 0.01 -xblmeth 1 -xblndig 13 -xbltol 0.002 -sort -csv",
                          cols, dat, "out.txt")
          
          nrand<-floor(runif(nboot,1,1000000)) # 32-bit positive integers
          for(i in 1:nboot){
            cat(paste("\b\b\b\b\b\b#",sprintf("%04i",i)))
            flush.console()
            j<-nrand[i]
            
            ### Update nlmeargs.txt file
            bootargs<-paste("/boot",j,"/bootsamp",i,"/boottry 1",bstrarg)
            
            if(i==1) write(c(bootargs,nlmeargs),"nlmeargs.txt")
            
            shellrun<-ifelse(hold,"cmd.exe /q /k","cmd.exe /q /c")
            
            try(system(paste(shellrun,"call \"./NLME7.exe\"",bootargs,nlmeargs,"1> err1.txt 2> err2.txt"),wait=TRUE,invisible=TRUE))
            
            ### Create output
            write(paste("\n### Bootstrap:",sprintf("%04i",i)),"out0002.txt",append=TRUE)
            if(file.exists("out.txt")) file.append("out0002.txt","out.txt")
          }
          if(file.exists("out.csv")) file.rename("out.csv","out0002.csv")
          
          write(c(bootargs,nlmeargs),"nlmeargs.txt")
          
          ### Delete executables
          if(clean){
            unlink(list.files(pattern="\\.exe$"))
          }
                            
          write("\nDone!","")
        }else{
          write("Error: Required files are missing. Check that .exe file is in folder.","")
        }
      }else{
        write("Error: A maximun of 3 stratification variables exceeded.","")
      }
    }else{
      write("Error: Number of bootstrap samples is outside allowed range of 1-9999.","")
    }
  }
}
