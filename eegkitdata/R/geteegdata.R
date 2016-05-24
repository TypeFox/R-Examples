geteegdata <-
  function(indir,outdir=indir,cond=c("S1","S2m","S2n"),nt=NULL,
           filename="eegdata",filetype=c(".rda",".csv",".txt")){
    ###### Create data matrix from UCI EEG data
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: September 9, 2014

    ### initial checks
    filetype=filetype[1]
    if(is.na(match(filetype,c(".rda",".csv",".txt")))){stop("Incorrect filetype input.")}
    cond=cond[1]
    if(is.na(match(cond,c("S1","S2m","S2n")))){stop("Incorrect cond input.")}
    if(is.null(nt)==FALSE){nt=as.integer(nt); if(nt<1){stop("Incorrect nt input.")}}
    filename=as.character(filename)
    ix=file.info(indir)
    if(is.na(ix[1]) | ix[2]==FALSE){stop("Invalid input directory (indir is not a valid directory).")}
    ox=file.info(outdir)
    if(is.na(ox[1]) | ox[2]==FALSE){stop("Invalid output directory (outdir is not a valid directory).")}
    
    ### load list of file directories (and possibly untar)
    alldir=list.dirs(indir,full.names=FALSE)[-1]
    if(length(alldir)==0L){
      alldir=list.files(indir)
      lad=length(alldir)
      cat(paste("untarring",lad,"files...\n"))
      for(j in 1:length(alldir)){untar(paste(indir,alldir[j],sep=""),exdir=indir)}
      alldir=list.dirs(indir,full.names=FALSE)[-1]
      if(length(alldir)==0L){stop("Invalid input directory (no data directories).")}
    }
    
    ### load all data
    eegdata=NULL
    for(j in 1:length(alldir)){
      cat(paste("subject:",alldir[j],"\n"))
      egroup=strsplit(as.character(alldir[j]),"")[[1]][4]
      thedir=paste(indir,alldir[j],sep="") # directory
      flist=list.files(path=thedir)         # list files
      flen=length(flist)                    # number of files
      k=m=0; maxk=FALSE
      while(k<flen & maxk==FALSE){
        k=k+1
        fn=paste(thedir,flist[[k]],sep="/")
        einfo=scan(file=fn,what=character(),skip=3,nlines=1,quiet=TRUE)
        condition=paste(einfo[2],einfo[3],sep="")
        if(condition=="S1obj"){
          condition="S1"
        } else if(condition=="S2match"){
          condition="S2m"
        } else if(condition=="S2nomatch,"){
          condition="S2n"
        } 
        if(condition==cond){
          eegtab=read.table(file=fn) # load data
          colnames(eegtab)=c("trial","channel","time","voltage")
          eegdata=rbind(eegdata,data.frame(subject=alldir[j],group=egroup,condition,eegtab))
          if(is.null(nt)==FALSE){m=m+1; if(m==nt){maxk=TRUE}}
        }
      } # end while(k<flen & maxk==FALSE)
    } # end for(j in 1:length(alldir))
    
    ### save data
    if(filetype==".rda"){
      save(eegdata,file=paste(outdir,filename,filetype,sep=""))
      return(eegdata)
    } else if(filetype==".csv"){
      write.csv(eegdata,file=paste(outdir,filename,filetype,sep=""))
    } else {
      write.table(eegdata,file=paste(outdir,filename,filetype,sep=""))
    }
    
  }