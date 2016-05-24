reGenotyper <- function(phenotype, genotype,fileName="test",thres=0.9,optGT=TRUE,optGTplot=FALSE,
                        optGT.thres=0,permu=FALSE,n.permu=10,wls.score.permu=NULL,process=TRUE,t.thres=1.5, GT.ref=NULL){
    #phenotype: a matrix containing trait value, nPheno x nSamples
    #genotype: a matrix containing genotype information( elements being 1 or 0 ), nSamples x nMakers  
    #fileName is the output file name
    #thres: probability threshold for deciding a sample being wrongly labeled (default=0.9)
    #optGT=TRUE: the optimal genotype will be recovered
    
 
    # 1. Checking arguments... 
    if(process){
      cat("STEP 1: Checking arguments... \n" )  
      flush.console()
    }
    
    if( is.null(colnames(genotype) ) )
        colnames(genotype) <- paste( "Sample", 1:ncol(genotype), sep="")
    if( is.null(rownames(genotype) ) )
        rownames(genotype) <- paste( "mk", 1:nrow(genotype), sep="")
	  if( is.null(colnames(phenotype) ) )
        colnames(phenotype) <- paste( "Sample", 1:ncol(phenotype), sep="") 
    if( is.null(rownames(phenotype) ) )
        rownames(phenotype) <- paste( "pheno", 1:nrow(phenotype), sep="")
			
		
   
    if(("A" %in% unique(as.vector(as.matrix(genotype)))) | ("B" %in% unique(as.vector(as.matrix(genotype)))))  {
        temp                <- genotype
        temp[genotype=="A"] <- 1
        temp[genotype=="B"] <- 0
        if("H" %in% unique(as.vector(as.matrix(genotype)))) {
           temp[genotype=="H"] <- 0.5
        }
        temp                <- matrix(as.numeric(as.matrix(temp)),
                                          nrow=nrow(genotype))
        dimnames(temp)      <- dimnames(genotype)
        if(process){
            cat("The gentoypes: A, (H) and B have been coded as 0, 0.5 and 1 in reGenotyper \n" )  
          flush.console()
        }
    }
    
    gt.unique <- unique(as.vector(as.matrix(genotype)))
    if(length(gt.unique)==2) {  #RILs, input gentoype has elements of "1" and "2"
        if("2" %in% gt.unique & "1" %in% gt.unique) genotype[genotype==2] <- 0
        if(process){
            cat("  There were two possible gentoypes: 1 and 2 in input genotypes, which have been coded as 0,1 in reGenotyper \n" )  
          flush.console()
        }
    }
    if(length(gt.unique)==3) {  #F2 or human, input gentoype has elements of "0", 1" and "2"
        if("0" %in% gt.unique & "1" %in% gt.unique & "2" %in% gt.unique) {
          genotype[genotype==1] <- 0.5
          genotype[genotype==2] <- 1
          if(process){
            cat("  There are three possible gentoypes: 0,1,2, which have been coded as 0,0.5,1 in reGenotyper \n" )  
            flush.console()
          }
        }
    }
    
    if(process){
      cat( "\n")
      cat("STEP 2: Computing QTL profiles using original genotype data... \n" )
      flush.console() 
    }
    
    if(nrow(phenotype)>8000){    #make the phenotype smaller: select 3000 traits with large variance
      var.all   <- apply(phenotype,1,var)
      var.thres <- sort(var.all,decreasing=TRUE)[3000]  
      pheno     <- phenotype[which(var.all>=var.thres),]
      if(process){
          cat("  The phenotypes with largest variance were selected and used in reGenotyper \n" )  
          flush.console()
      }
    } else{
      pheno  <- phenotype
    }
    gt		 <- genotype
    N      <- ncol(gt)       
    if( file.exists(paste("t_mat0_",fileName,".Rdata", sep="")) )
    {
        # We load the data previously calculated
        if(process){
          cat( ">> We load the data previously calculated \n" )
          flush.console()
        }
        load(file=paste("t_mat0_",fileName,".Rdata", sep=""))
    }
    #   If not, we had to compute it
    else 
    {
       t.mat0 <- tMatFunction(pheno,gt,fileName) # matrix probes x marker
       if (!is.null(fileName)) {save(t.mat0, file=paste("t_mat0_",fileName,".Rdata", sep=""))}   
    }
    if(process){
      cat("\n")
      cat("STEP 3: Computing wrongly labeled sample score... \n" )
      flush.console()
    }
    if( file.exists(paste("delta_t_mat_allmk_list_",fileName,".Rdata", sep="") ))
    {
        # We load the data previously calculated
        if(process){
          cat( ">> We load the data previously calculated\n" )
          flush.console()
        }
        load(paste("delta_t_mat_allmk_list_",fileName,".Rdata", sep=""))
    }
    #   If not, we had to compute it
    else 
    {                                   
        delta.t.mat.allmk.list  <- .deltaTMatAllmk_list(gt,pheno,t.mat0,indSample=NULL,t.thres=t.thres,fileName=fileName)#a list with length=nMK, each element is a matrix nSample by nSigGene ,
        # each element matrix is deltaT for each sigGene(row) when  sample i (column) is flipped
        if (!is.null(fileName)) {save(delta.t.mat.allmk.list, file=paste("delta_t_mat_allmk_list_",fileName,".Rdata", sep=""))}
      
    }      
    
    deltaT.sampleByall <- NULL   #row:samples, all: all markers for sig genes at each marker
    for(k in 1:N){
      deltaT.sampleByall <-  rbind(deltaT.sampleByall ,unlist(lapply(delta.t.mat.allmk.list,function(x) x[,k])) )
    }
      
    #calculat ethe area of deltaT>0
    wls.score      <- apply( deltaT.sampleByall,1, .area.my)
    names(wls.score) <- colnames(gt)

    
    #get p value for wls score using permutation
    wls.pValue      <- NULL
    wls.ind         <- NULL
    wls.names       <- NULL
    if(is.null(wls.score.permu) & permu==FALSE ){
      
      if(process){ 
        cat( " Plotting WLS score... \n\n" )
        flush.console()
      }
      mycol         <- rep("black",length(wls.score))
      plot(wls.score, ylim=c(min(wls.score)*0.95, max(wls.score)*1.1),ylab="wls.score",xlab="sample",
           pch=19,cex.lab=1.2,cex.axis=1.2,col=mycol,main="")          
      
      cat("In order to detect wrongly labeled samples using wls.score, function \"wls.score.permu\" needs to be calculated by setting permu=TRUE. An alternative is to use the function \"permutation\", e.g.
           wls.score.permu <- permutation(phenotype,genotype,n.permu=1000,process=TRUE,fileName\"test\",t.thres=3)")
            
      return(list(wls.score=wls.score,wls.names=NULL, gt.opt=NULL,wls.pValue=NULL,wls.score.permu=NULL))
    
      } else{ # compute wls.score.permu or use input wls.score.permu to do following steps

        if(permu==TRUE) {
          if(process){ 
            cat( "\n STEP 4: Permutation... \n" )
            flush.console()
          }
          
          wls.score.permu      <- permutation(pheno,gt,n.permu,process,fileName,t.thres) 
        } else{
          if(process){ 
            cat( "\n")
            cat("STEP 4: Input variable \"wls.score.permu\" from earlier permutation is used ... \n" )
            flush.console()
          }
        }
      
        for (i in 1:length(wls.score)){
          wls.pValue <- c(wls.pValue, length(which( as.vector(wls.score.permu) < wls.score[i] ))/length(as.vector(wls.score.permu)))
        }     
        wls.ind   <- which (wls.pValue >=thres) 
        wls.names <- colnames(gt)[wls.ind]
          
        #plot permutation results         
        #pdf(file="permutation_plot.png")
        if(process){ 
          cat("\n")
          cat("STEP 5: Plotting... \n" )
          wls.output <- list(wls.score=wls.score,wls.names=wls.names, gt.opt=1, wls.pValue=wls.pValue,wls.score.permu=wls.score.permu,thres=thres)
          class(wls.output)[2] <- "wls"
          plot.wls(wls.output)
          flush.console()
        }
      
        #print the identified WLS
         if(process){
          if(length(wls.ind>0)){
            cat(paste("RESULTS:", length(wls.ind), "samples have been detected as wrongly labled samples.",sep=" "),"\n")
            cat("They are: ");cat(wls.names,sep="; ","\n")
            
            cat("with score of ");cat(wls.score[wls.ind],sep="; ","\n")
            cat("with probabilty of ");cat(wls.pValue[wls.ind],sep="; ","\n\n")
            flush.console()
          } else{      
            cat("Luckily, no wrongly labled sample has been found! \n\n" )
            flush.console()
          }
        }
        
      
        # step 5 get optimal gt 
        gt.opt            <- NULL
        if( optGT==TRUE & length(wls.ind)>0 ){
          if(process){
            cat("\n")
            cat("STEP 5: Recovering the optimal genotype for the detected wrongly labeled samples...\n " )
            cat(" >> Please find it in the 3rd element of the output list from reGenotyper function.\n " )           
            flush.console()
          }
          
          gt.opt.allSamples <- optimalGT(delta.t.mat.allmk.list,gt, gt.thres=optGT.thres,optGTplot=optGTplot)     
          gt.opt  <- gt.opt.allSamples[,wls.ind]
          if(length(wls.ind)==1) {
            gt.opt          <- as.matrix(gt.opt)
            colnames(gt.opt)<- colnames(gt.opt.allSamples)[wls.ind]
          }
          if (!is.null(fileName)) {
            save(gt.opt,gt.opt.allSamples, file=paste("gt_opt_",fileName,".Rdata", sep="")) 
          }
          
          #step 6 find best matched genotype for detected WLS
          if(process){
            cat("\n")
            cat("STEP 6: Finding the best matched genotype for the detected wrongly labeled samples... \n " )
            flush.console()
          }
          if(is.null(GT.ref)) GT.ref <- genotype
          cat(paste("[WLS name]", "[Top 5 best matched genotype name]", "[percentage of common genotype between them]", sep="\t\t"),sep="\n")
          for(i in 1:length(wls.ind)){
            gt.diff <- 100*apply(abs(gt.opt[,i]-GT.ref),2,function(x) sum(x, na.rm =TRUE))/nrow(GT.ref) 
            cat(colnames(gt)[wls.ind[i]],names(sort(gt.diff))[1:5],100-sort(gt.diff)[1:5], sep="\t" ,"\n")      
          }        
        }
        wls.output <- list(wls.score=wls.score,wls.names=wls.names, gt.opt=gt.opt,wls.pValue=wls.pValue,wls.score.permu=wls.score.permu,thres=thres)
        class(wls.output)[2] <- "wls"
        #print output
        if(process ){  
          print.wls(wls.output)          
          flush.console()
        }
        invisible(wls.output)     
      }

}