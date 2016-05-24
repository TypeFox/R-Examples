enviPickwrap <-
function(             filepath.mzXML,
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
                      recurs=3,
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
      if(progbar==TRUE){    prog<-winProgressBar("Peak Picking",min=0,max=5);
                            setWinProgressBar(prog, 1, title = "Peak Picking - reading scans", label = NULL);}                      
      MSlist<-readMSdata(filepath.mzXML,MSlevel,progbar=FALSE,minRT=FALSE,maxRT=FALSE,minmz=FALSE,maxmz=FALSE);   
      if(progbar==TRUE){setWinProgressBar(prog, 2, title = "Peak Picking - Partition", label = NULL);}       
      MSlist<-mzagglom(MSlist,dmzgap,ppm,drtgap,minpeak,maxint,progbar=FALSE)     
      if(progbar==TRUE){setWinProgressBar(prog, 3, title = "Peak Picking - EIC cluster", label = NULL);}
      MSlist<-mzclust(MSlist,dmzdens,ppm,drtdens,minpeak,merged=TRUE,maxint,progbar=FALSE,from=FALSE,to=FALSE) 
      if(progbar==TRUE){setWinProgressBar(prog, 4, title = "Peak Picking - Peaks in EICs", label = NULL);}            
      MSlist<-mzpick(MSlist,minpeak,drtsmall,drtfill,drttotal,recurs,weight,SB,SN,minint,maxint,
                    ended,progbar=FALSE,from=FALSE,to=FALSE)
      if(progbar==TRUE){setWinProgressBar(prog, 5, title = "Peak Picking - done.", label = NULL);
                        Sys.sleep(0.5);close(prog);} 
      ##########################################################################
      return(MSlist);

}
