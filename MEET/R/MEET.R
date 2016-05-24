MEET <-function(TF,nameTF=NULL,seqin,alg=NULL,method,mode,org,vector=NULL,num_motif=NULL,len_motif=NULL,direction="f",threshold=0.1,order=NULL,model=NULL, position=c(501),mv=50,gapopen=-500,maxiters=16,gapextend=-2, optionsFile='.optionsmeet'){

 
      require("seqinr")
      #require("fields")
      require("MEET")
      require("methods")

      checkpaths <- function(executable){
        out <- switch(.Platform$OS.type,
               unix= system(paste('which',executable), intern = TRUE),
                      
               windows=)
      }

      executables <- c("clustalw","meme","mast","transfac2meme","muscle","MDScan-linux")
      names(executables) <-  c("call.clustalw","call.meme", "all.mast","call.transfac2meme","call.muscle", "call.MDscan")
      
# Initialitzation code, this checks whether MEME/MAST, clustalx and muscle is present in the system

      if ( is.na(file.info(optionsFile)$size) )
        {
          
          cat( "Configuration File (",optionsFile,") not found\n")
          cat ("\tChecking external dependencies\n")
          found <- sapply( executables, function(x) checkpaths(x))
          meetconf <- sapply(names(found), function(x) {
            if (length(found[[x]]) == 0 ) {
              done <- 0
              while(!done){
              cat(executables[x], ' not found, please enter new path, or ENTER for not available:\n')
              cat('\t ',executables[x],' path: ')
              npath <- readline()
              if ( length(npath) == 1 ) {done <- 1}
              if ( !is.na(file.info(npath)$size) )  {done <- 1}
            }
              npath
            } else {
              cat("\tFound",x," in ", found[[x]],"\n")
              found[[x]]
            }
          })
          write.table(as.data.frame(meetconf),file=optionsFile,sep='=',col.names=FALSE)
        } else {
          meetconftmp <- read.table(file=optionsFile,sep='=',header=FALSE)
          meetconf <- as.character( meetconftmp$V2)
          names(meetconf) <- as.character(meetconftmp$V1)
        }

      options(warn=-1)
    
     data(organism,package="MEET", envir=environment())
    
     organism<-get("organism")
    
      write.fasta <- get("write.fasta",pos="package:seqinr")
      read.fasta <- get("read.fasta",pos="package:seqinr")
     
      iicc<-vector(mode="list", length=28)
      
      names(iicc)<-c("mode", "method","nameTF","organism", "background","alignment","pvalue","parameters","model","Transcriptionfactor",
      "nummotif","lenmotif","direction","outsequence","position","missing","vector","gapopen","maxiters","gapextend","clustalw","muscle","meme","mast","transfac2meme","Mdscan","DNA", "parameterIdeal")
	
      iicc$mode=mode
      iicc$method=method
      iicc$nameTF=nameTF
      iicc$organism=org
      iicc$background= organism[org,]
      iicc$alignment=alg
      iicc$threshold<-threshold
      iicc$parameters<-order
      iicc$model<-model
      iicc$nummotif<-num_motif
      iicc$lenmotif<-len_motif
      iicc$direction<-direction
      iicc$outsequence<-sequence
      iicc$position<-position
      iicc$missing<-mv     
      iicc$vector<-vector
      iicc$gapopen<-gapopen
      iicc$maxiters<-maxiters
      iicc$gapextend<-gapextend
	  iicc$clustalw<-as.character(meetconf["call.clustalw"])
	  iicc$muscle<-as.character(meetconf["call.muscle"])
	  iicc$meme<-as.character(meetconf["call.meme"])
	  iicc$mast<-as.character(meetconf["call.mast"])
	  iicc$transfac2meme<-as.character(meetconf["call.transfac2meme"])
	  iicc$MDscan<-as.character(meetconf["call.MDscan"])
      iicc$parameterIdeal<-NULL
	
	if(iicc$direction!="b"){ iicc$DNA<-vector(mode="list", length=1)
      
		  if(iicc$direction=="f"){
			  iicc$DNA[[1]]<-read.fasta(file=seqin,forceDNAtolower=FALSE, set.attributes=FALSE)[[1]]
			  iicc$DNA[[1]][iicc$DNA[[1]]=="a"]="A"; iicc$DNA[[1]][iicc$DNA[[1]]=="t"]="T";
			  iicc$DNA[[1]][iicc$DNA[[1]]=="c"]="C"; iicc$DNA[[1]][iicc$DNA[[1]]=="g"]="G";
		  }else{
			  iicc$DNA[[1]]<-read.fasta(file=seqin,forceDNAtolower=FALSE, set.attributes=FALSE)[[1]]
			  iicc$DNA[[1]]<-rev(iicc$DNA[[1]])
			  iicc$DNA[[1]][iicc$DNA[[1]]=="a"]="A"; iicc$DNA[[1]][iicc$DNA[[1]]=="t"]="T";
			  iicc$DNA[[1]][iicc$DNA[[1]]=="c"]="C"; iicc$DNA[[1]][iicc$DNA[[1]]=="g"]="G";
			 }
		 
	  }else{
			  iicc$DNA<-vector(mode="list", length=2)
			  iicc$DNA[[1]]<-read.fasta(file=seqin,forceDNAtolower=FALSE, set.attributes=FALSE)[[1]]
		      iicc$DNA[[1]][iicc$DNA[[1]]=="a"]="A"; iicc$DNA[[1]][iicc$DNA[[1]]=="t"]="T";
		      iicc$DNA[[1]][iicc$DNA[[1]]=="c"]="C"; iicc$DNA[[1]][iicc$DNA[[1]]=="g"]="G";
			  
	          iicc$DNA[[2]][iicc$DNA[[1]]=="A"]="T";iicc$DNA[[2]][iicc$DNA[[1]]=="T"]="A";
			  iicc$DNA[[2]][iicc$DNA[[1]]=="C"]="G";iicc$DNA[[2]][iicc$DNA[[1]]=="G"]="C"
		  
			  iicc$DNA[[2]]<-rev(iicc$DNA[[2]])
	  }

      background<-iicc$background
     
      write.table(background, "bfile", col.names=FALSE, quote=FALSE)
      
      output<-vector(mode="list", length=3)
      names(output)<-c("Consensus", "Summary", "Results")

      if (is.null(iicc$nameTF)=="TRUE"){
             
                alignedSequences<-Alignment(TF,iicc)
                iicc$Transcriptionfactor<-alignedSequences
                output$Consensus<-CreateConsensus(alignedSequences,iicc,filein=TF )
                iicc$Consensus<-output$Consensus

            }else{
            
            
            
                iicc$Transcriptionfactor<-NULL
                iicc$Consensus<-NULL
                
        
        }   
          
	  output$Results<-switch(iicc$mode, "training"=ConstructModel(iicc,TF), "detection"=detection(iicc))
      output$Summary<-iicc
    return(output)
}

