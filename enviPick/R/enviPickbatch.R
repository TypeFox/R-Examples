enviPickbatch <-
function(
                      folderin,
                      folderout,
                      MSlevel=c(1),
                      dmzgap=15,  
                      dmzdens=4,       
                      ppm=TRUE,
                      drtgap=500, 
                      drtsmall=20,
                      drtdens=250,
                      drtfill=10,
                      drttotal=200,
                      minpeak=4,
                      recurs=10,
                      weight=2,
                      SB=3,
                      SN=2,
                      minint=1E5,
                      maxint=1E7,
                      ended=2,
                      progbar=FALSE

      ){

      ##########################################################################
      # checks & setups ########################################################
      if(!is.integer(minpeak)){minpeak<-ceiling(minpeak);minpeak<-as.integer(minpeak)}
      if(minpeak<=1){stop("minpeak must be >1")}
      if(!is.loaded("agglom")){stop(".Call to agglom failed; aborted.")}
      if(!is.loaded("indexed")){stop(".Call to indexed failed; aborted.")}
      if(!is.loaded("densvec")){stop(".Call to densvec failed")}
      if(!is.loaded("pickpeak")){stop(".Call to pickpeak failed!")};
      if(!is.loaded("gapfill")){stop(".Call to gapfill failed!")};
      if(!is.loaded("picklist")){stop(".Call to picklist failed!")};
      if(!is.loaded("getEIC")){stop(".Call to getEIC failed")}      
      if(!is.logical(ppm)){stop("invalid ppm argument, TRUE or FALSE")}  
      if(!is.logical(progbar)){stop("invalid progbar argument, TRUE or FALSE")}       
      if(dmzdens<=0 || drtdens<=0){stop("invalid: drtdens or dmzdens <= 0")}
      if(drtsmall<=0 || drtsmall<=0){stop("drt must be >0!")};
      if(minint>=maxint){stop("Revise your minint & maxint settings!")}
      if(ended<1){stop("Wrong ended argument!")}    
      if(!is.integer(recurs)){recurs<-ceiling(recurs);recurs<-as.integer(recurs)}
      if(!is.integer(ended)){ended<-ceiling(ended);ended<-as.integer(ended)}
      ##########################################################################
      filed<-list.files(folderin);
      if(progbar==TRUE){    prog<-winProgressBar("Peak Picking",min=0,max=5);}
      for(n in 1:length(filed)){
        if(progbar==TRUE){    setWinProgressBar(prog, n, title = paste("Peak Picking - file #",n), label = NULL);}                      
			MSlist<-enviPickwrap(     
                      file.path(paste(folderin,"\\",filed[n],sep="")),
                      MSlevel=c(1),
                      dmzgap=5,  
                      dmzdens=4,       
                      ppm=TRUE,
                      drtgap=1000, 
                      drtsmall=20,
                      drtdens=250,
                      drtfill=10,
                      drttotal=200,
                      minpeak=4,
                      recurs=10,
                      weight=2,
                      SB=3,
                      SN=2,
                      minint=10E4,
                      maxint=10E6,
                      ended=2,
                      progbar=FALSE)              
			file_out<-paste(filed[n],"_picked",sep="")
			writePeaklist(MSlist,folderout,file_out)
      }
      if(progbar==TRUE){setWinProgressBar(prog, n, title = "Peak Picking - done.", label = NULL);
                        Sys.sleep(0.5);close(prog);}
      ##########################################################################                         

}
