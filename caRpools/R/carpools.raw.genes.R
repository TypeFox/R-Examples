carpools.raw.genes=function(untreated.list,treated.list, genes=NULL, namecolumn=1, fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), do.plot=TRUE, log=FALSE, put.names=FALSE, type="foldchange", controls.target= NULL, controls.nontarget=NULL, sort=TRUE)
  {
  # TYPE can be
  # foldchange -> plots fold change without error bars
  # readcount -> plots readcount with error bars
  # z-score -> plots z-score normalized data
  
  # for Z-ratio: highest quartil -> plot a line as might be significant
  
  # reset PAR
  #par(mfrow=c(1,1))
  
  par(mar=c(5,4,4,2))
  par(oma=c(2,2,2,2))
  
  
  ## Start calculations for Foldchange plotting
  if(type=="foldchange")
  {
  stopifnot(length(genes)<=16)
  
  untreated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1  
       
      )
    ),
    1,
    mean)
  treated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    mean
  )
  untreated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    sd)
  treated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    sd
  )
  gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
  designs=treated.list[[1]][,namecolumn]
  
  
  dataset.combined <- data.frame( 
    untreated = untreated.list.mean,
    treated = treated.list.mean,
    untreated.sd = untreated.list.sd,
    treated.sd = treated.list.sd,
    genes=gene.names,  
    log2foldchange=log2(as.numeric(treated.list.mean)/as.numeric(untreated.list.mean)),  
    designs=as.character(designs),
    row.names = designs, 
    stringsAsFactors=FALSE
    
   
  )

# plot data
    if(do.plot){
   
      for(i in genes){   
        
        # Sort according to fold change
        dataset.combined.genes = dataset.combined[dataset.combined$genes==i,]
        if(identical(sort, TRUE))
        {
          dataset.combined.genes = dataset.combined.genes[order(dataset.combined.genes$log2foldchange, decreasing=TRUE),]
        }
       
        fc = dataset.combined.genes$log2foldchange
        #fc=log2(dataset.combined$treated[dataset.combined$genes==i]/dataset.combined$untreated[dataset.combined$genes==i])
        #str(dataset.combined.genes)
        col.v=rep(rgb(46,98,166, 255, maxColorValue=255),times=length(fc))
        col.v[fc>0]=rgb(217,35,35, 255, maxColorValue=255)
        a<-barplot(dataset.combined.genes$log2foldchange,
                  col=col.v,
                  main=paste("Foldchange:",i, sep=" "),
                  ylab="Normalized log2(foldchange)",
                  names.arg=as.character(dataset.combined.genes$designs),
                  las=2, cex.names=0.6)
        #str(if(put.names){rownames(dataset.combined.genes)}else{NULL})
      }
      #dev.off()
      #par(mfrow=c(1,1),las=1) 

    }
      # return data table for this gene
      res=do.call("rbind.data.frame",tapply(genes,genes,function(x) return(dataset.combined[dataset.combined$genes==x,]),simplify=FALSE))
      return(res)

  } #End if type==foldchange



## Start calculations for Foldchange plotting
if(type=="foldchange.vioplot")
{
  stopifnot(length(genes)<=16)
  
  
  untreated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1  
        
      )
    ),
    1,
    mean)
  treated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    mean
  )
  untreated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    sd)
  treated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    sd
  )
  gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
  designs=treated.list[[1]][,namecolumn]
  
  
  dataset.combined <- data.frame( 
    untreated = untreated.list.mean,
    treated = treated.list.mean,
    untreated.sd = untreated.list.sd,
    treated.sd = treated.list.sd,
    genes=gene.names,  
    log2foldchange=log2(as.numeric(treated.list.mean)/as.numeric(untreated.list.mean)),  
    row.names = designs, 
    stringsAsFactors=FALSE
    
    
  )

  # plot data
  if(do.plot){
    
    for(i in genes){   
      
      
      fc=log2(dataset.combined$treated[dataset.combined$genes==i]/dataset.combined$untreated[dataset.combined$genes==i])
      min.y = min(fc)
      max.y = max(fc)
      # first will be GENE to plot
      sg.all = log2(dataset.combined$treated/dataset.combined$untreated)
      #sg.all = apply(sg.all,1, function(z) { if(is.infinite(z[1]) || is.na(z[1])) { return(as.numeric(0)) } else { return(as.numeric(z[1])) } } )
      
      
      if(!is.null(controls.nontarget))
      {
        sg.nontarget = log2(dataset.combined$treated[dataset.combined$gene %in% controls.nontarget]/dataset.combined$untreated[dataset.combined$gene %in% controls.nontarget])
      }
      else
      {
        sg.nontarget=sg.all
      }
      if(!is.null(controls.target))
      {
        sg.target = log2(dataset.combined$treated[dataset.combined$gene %in% controls.target]/dataset.combined$untreated[dataset.combined$gene %in% controls.target])
      }
      else
      {
        sg.target=sg.all
      }
      
      sg.list = list(gene = fc, all = sg.all, nontarget = sg.nontarget, target = sg.target)
      
      # remove NA / infinite
      lapply(sg.list, function(z) {
                                         if(is.infinite(z) || is.na(z) || is.null(z)) { return(as.numeric(0)) } else { return(as.numeric(z)) } } )
      fc = sg.list$gene
      sg.all = sg.list$all
      sg.nontarget = sg.list$nontarget
      sg.target = sg.list$target
      
      #sg.all = lapply(sg.list, function(z) { if(is.infinite(z["all"]) || is.na(z["all"])) { return(as.numeric(0)) } else { return(as.numeric(z["all"])) } } )
      #sg.nontarget = lapply(sg.list, function(z) { if(is.infinite(z["nontarget"]) || is.na(z["nontarget"])) { return(as.numeric(0)) } else { return(as.numeric(z["nontarget"])) } } )
      #sg.target = lapply(sg.list, function(z) { if(is.infinite(z["target"]) || is.na(z["target"])) { return(as.numeric(0)) } else { return(as.numeric(z["target"])) } } )
      
      if(min(c(sg.target,sg.all,sg.nontarget)) < min.y) { min.y = min(c(sg.target,sg.all,sg.nontarget)) }
      if(max(c(sg.target,sg.all,sg.nontarget)) > max.y) { max.y = max(c(sg.target,sg.all,sg.nontarget)) }
      
      #par(mfrow=c(1,1))
      #par(mar=c(5,4,4,2))
      if(is.null(controls.nontarget))
      {
        if(is.null(controls.target))
        {
          carpools.vioplot(fc, sg.all, range=1.5, ylim=c(min.y,max.y), names=c(i,"all"), horizontal=FALSE,
                   col=c("orange","grey"), border="black", lty=1, lwd=1,
                   colMed="white", pchMed=16, add=FALSE, wex=1,
                   drawRect=TRUE)
          title(paste("Log2 Foldchange of sgRNAs for",i, "compared to sgRNA from all genes", sep=" ")) 
        }
        else
        {
          carpools.vioplot( fc, sg.all, sg.target, range=1.5, ylim=c(min.y,max.y), names=c(i,"all","positive"), horizontal=FALSE,
                   col=c("orange","grey","red"), border="black", lty=1, lwd=1,
                   colMed="white", pchMed=16, add=FALSE, wex=1,
                   drawRect=TRUE)
          title(paste("Log2 Foldchange of sgRNAs for",i, "compared to sgRNA from all genes and targeting controls", sep=" ")) 
          
        }
      }
      else if(is.null(controls.target))
      {
        carpools.vioplot( fc, sg.all, sg.nontarget, range=1.5, ylim=c(min.y,max.y), names=c(i,"all","non-targeting"), horizontal=FALSE,
                 col=c("orange","grey","blue"), border="black", lty=1, lwd=1,
                 colMed="white", pchMed=16, add=FALSE, wex=1,
                 drawRect=TRUE)
        title(paste("Log2 Foldchange of sgRNAs for",i, "compared to sgRNA from all genes and non-targeting controls", sep=" ")) 
      }
      else
      {
        carpools.vioplot( fc, sg.all, sg.nontarget, sg.target, range=1.5, ylim=c(min.y,max.y), names=c(i,"all","non-targeting","positive"), horizontal=FALSE,
                 col=c("orange","grey","blue","red"), border="black", lty=1, lwd=1,
                 colMed="white", pchMed=16, add=FALSE, wex=1,
                 drawRect=TRUE)
        title(paste("Log2 Foldchange of sgRNAs for",i, "compared to sgRNA from all, non-targeting and targeting controls", sep=" ")) 
      }
    
      
      
    }
  
    
    
    
    
    
    
    
  }
  # return data table for this gene
  res=do.call("rbind.data.frame",tapply(genes,genes,function(x) return(dataset.combined[dataset.combined$genes==x,]),simplify=FALSE))
  return(res)
  
} #End if type==foldchange.vioplot
  

  # plot readcount instead of foldchange data
  if(type=="readcount")
  {
    
    
    untreated.list.mean=apply(
      do.call(
        "cbind",
        lapply(
          untreated.list,
          function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))
        )
      ),
      1,
      mean)
    treated.list.mean=apply(
      do.call(
        "cbind",
        lapply(
          treated.list,
          function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))
        )
      ),
      1,
      mean
    )
    untreated.list.sd=apply(
      do.call(
        "cbind",
        lapply(
          untreated.list,
          function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))
        )
      ),
      1,
      sd)
    treated.list.sd=apply(
      do.call(
        "cbind",
        lapply(
          treated.list,
          function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))
        )
      ),
      1,
      sd
    )
    
    gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
    designs=treated.list[[1]][,namecolumn]
    
    
    dataset.combined <- data.frame( 
      untreated.readcount = as.numeric(untreated.list.mean),
      treated.readcount = as.numeric(treated.list.mean),
      untreated.sd = as.numeric(untreated.list.sd),
      treated.sd = as.numeric(treated.list.sd),
      genes=gene.names,  
      row.names = designs, 
      stringsAsFactors=FALSE

    )
    
    # check dataset SD for NA, INF and set to 0
    dataset.combined$untreated.sd = apply(dataset.combined, 1, function(z) {
      if(is.na(z["untreated.sd"]) || is.infinite(z["untreated.sd"]))
      {
        return(as.numeric(0))
      }
      else
      {
        return(as.numeric(z["untreated.sd"]))
      }
    })
    
    dataset.combined$treated.sd = apply(dataset.combined, 1, function(z) {
      if(is.na(z["treated.sd"]) || is.infinite(z["treated.sd"]))
      {
        return(as.numeric(0))
      }
      else
      {
        return(as.numeric(z["treated.sd"]))
      }
    })
    # plot data
    
    if(do.plot){
  
      # for each gene
      for(i in genes){ 
 
       # Treated Dataset
       data.plot = dataset.combined[dataset.combined$genes==i,]
           
            for (m in 1:nrow(data.plot))
           {

             if(m==1)
             {
               df.plot <- data.frame(
                 sgRNA = c(data.plot[m,"untreated.readcount"], data.plot[m, "treated.readcount"]),
                 stringsAsFactors=FALSE)
               sd.plot <- data.frame(
                 sgRNA = c(data.plot[m,"untreated.sd"], data.plot[m, "treated.sd"]),
                 stringsAsFactors=FALSE)
             }
             else
             {
             
               df.plot[,m] = as.numeric(c(data.plot[m,"untreated.readcount"], data.plot[m, "treated.readcount"])) 
               sd.plot[,m] = as.numeric(c(data.plot[m,"untreated.sd"], data.plot[m, "treated.sd"])) 
               
             }

           }
           colnames(df.plot) = rownames(data.plot)
           rownames(df.plot) = c("untreated","treated")
           # put in matrix
           df.plot = as.matrix(df.plot)
           sd.plot = as.matrix(sd.plot)
           
           maxy = max(df.plot + sd.plot)
           miny = min(df.plot - sd.plot)
           
           
           a<-barplot(df.plot,
                      col=c(rgb(46,98,166, 255, maxColorValue=255),rgb(217,35,35, 255, maxColorValue=255)),
                      main=paste("Readcount:",i, sep=" "),
                      ylab="Normalized Mean Readcount",
                      legend = rownames(df.plot),
                      ylim= c(0, maxy),
                      bty="n",
                      border=FALSE,
                      #ylim = c(0, maxy),
                      #names.arg=if(put.names){row.names(dataset.combined[dataset.combined$genes==i,])}else{NULL},
                      las=2, beside=TRUE,  cex.names=0.6)
           segments(a,df.plot-sd.plot,a,df.plot+sd.plot)
           segments(a-0.5,df.plot+sd.plot,a+0.5,df.plot+sd.plot)
           segments(a-0.5,df.plot-sd.plot,a+0.5,df.plot-sd.plot)
           
      }
      
    } else
    {
    # return data table for this gene
    res=do.call("rbind.data.frame",tapply(genes,genes,function(x) return(dataset.combined[dataset.combined$genes==x,]),simplify=FALSE))
    return(res)
    }
  }

## Start calculations for Z-Score plotting
if(type=="z-ratio")
{
  
  #Normalize each single dataset first with the norm.function
  untreated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    mean)
  treated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    mean
  )
  untreated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    sd)
  treated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))+1
      )
    ),
    1,
    sd
  )
  gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
  designs=treated.list[[1]][,namecolumn]
  
  # calculate final data.frame with Z-Score for each sgRNA
  # Z score is calculated for untreated and treated
  # then a z-ratio is calculated instead of a fold change
  zscore.untreated = (untreated.list.mean-mean(untreated.list.mean))/sd(untreated.list.mean)
  zscore.treated = (treated.list.mean-mean(treated.list.mean))/sd(treated.list.mean)
  
  #zscore.untreated=log2(zscore.untreated)
 # zscore.treated=log2(zscore.treated)
  
  # Z-Ratio according to 
  # Cheadle, C., Vawter, M. P., Freed, W. J., & Becker, K. G. (2003).
  # Analysis of microarray data using Z score transformation.
  zratio.initial = zscore.treated-zscore.untreated
  
  #zratio.initial = zscore.treated/zscore.untreated
  
  
  zratio = zratio.initial/sd(zratio.initial)
 
  # calculate highest quartile
  #quantiles = quantile(zratio)[4]
  dataset.combined <- data.frame( 
    genes=gene.names,
    untreated = untreated.list.mean,
    treated = treated.list.mean,
    untreated.sd = untreated.list.sd,
    treated.sd = treated.list.sd,
    z.score.untreated = zscore.untreated,
    z.score.treated  = zscore.treated,
    z.ratio=zratio/sd(zratio),  
    row.names = designs,
    designs = designs,
    col=NA,
    stringsAsFactors=FALSE  
  )
  
  
  # Color bars according to z-ratio
  dataset.combined$col = apply(dataset.combined, 1, function(i){
    if(as.numeric(i["z.ratio"]) < 0 )
    { return(rgb(46,98,166, 255, maxColorValue=255))
    }
    else
    {
      return(rgb(217,35,35, 255, maxColorValue=255))
    }
  })
   

  # plot data
  if(do.plot){
    

      for(i in genes){   
        
        # Sort according to fold change
        dataset.combined.genes = dataset.combined[dataset.combined$genes==i,]
        if(identical(sort, TRUE))
        {
          dataset.combined.genes = dataset.combined.genes[order(dataset.combined.genes$z.ratio, decreasing=TRUE),]
        }
        
        zr = dataset.combined.genes$z.ratio
        
        #zr=dataset.combined$z.ratio[dataset.combined$genes==i]
        col.v=dataset.combined.genes$col
        #col.v[zr>0]=rgb(217,35,35, 255, maxColorValue=255)
        
        a<-barplot(zr,
                   col=col.v,
                   main=paste("Z-Ratio:",i, sep=" "),
                   ylab="Normalized Z-Ratio",
                   names.arg=if(put.names){row.names(dataset.combined.genes)}else{NULL},
                   ylim = c(min(zr),max(zr)*1.2),
                   las=2, cex.names=0.6)
        }
    
 
  }
  # return data table for this gene
  res=do.call("rbind.data.frame",tapply(genes,genes,function(x) return(dataset.combined[dataset.combined$genes==x,]),simplify=FALSE))
  return(res)
  
} #End if type==z-score
 
# plot z-scores
if(type=="z-score")
{
  stopifnot(length(genes)<=16)
  
  untreated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (((x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))-mean(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))/sd(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))
      )
    ),
    1,
    mean)
  treated.list.mean=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (((x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))-mean(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))/sd(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))
      )
    ),
    1,
    mean
  )
  untreated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        untreated.list,
        function(x) (((x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))-mean(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))/sd(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))
      )
    ),
    1,
    sd)
  treated.list.sd=apply(
    do.call(
      "cbind",
      lapply(
        treated.list,
        function(x) (((x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn]))-mean(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))/sd(x[,fullmatchcolumn]/norm.function(x[,fullmatchcolumn])))
      )
    ),
    1,
    sd
  )
  
  gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
  designs=treated.list[[1]][,namecolumn]
  
  # calculate final data.frame with Z-Score for each sgRNA
  # Z score is calculated for untreated and treated
  # then a z-ratio is calculated instead of a fold change
  
  dataset.combined <- data.frame( 
    genes=gene.names,
    z.score.untreated = untreated.list.mean,
    z.score.treated = treated.list.mean,
    z.score.untreated.sd = untreated.list.sd,
    z.score.treated.sd = treated.list.sd,
    row.names = designs, 
    stringsAsFactors=FALSE
    
  )
  

  
  # plot data
  if(do.plot){
    
    # for each gene
    for(i in genes){ 
      
      
        data.plot = dataset.combined[dataset.combined$genes==i,]
        
        for (m in 1:nrow(data.plot))
        {
          if(m==1)
          {
            df.plot <- data.frame(
              sgRNA = c(data.plot[m,"z.score.untreated"], data.plot[m, "z.score.untreated"]),
              stringsAsFactors=FALSE)
            sd.plot <- data.frame(
              sgRNA = c(data.plot[m,"z.score.untreated.sd"], data.plot[m, "z.score.treated.sd"]),
              stringsAsFactors=FALSE)
            
          }
          else
          {
            df.plot[,m] = as.numeric(c(data.plot[m,"z.score.untreated"], data.plot[m, "z.score.treated"])) 
            sd.plot[,m] = as.numeric(c(data.plot[m,"z.score.untreated.sd"], data.plot[m, "z.score.treated.sd"])) 
            
          }
          
        }
        colnames(df.plot) = rownames(data.plot)
        rownames(df.plot) = c("untreated","treated")
        # put in matrix
        df.plot = as.matrix(df.plot)
        sd.plot = as.matrix(sd.plot)
        
        maxy = max(df.plot + sd.plot)
        miny = min(df.plot - sd.plot)
        

        
        
        a<-barplot(df.plot,
                   col=c(rgb(46,98,166, 255, maxColorValue=255),rgb(217,35,35, 255, maxColorValue=255)),
                   main=paste("Z-Score:",i, sep=" "),
                   ylab="Normalized Mean Z-Score",
                   legend = rownames(df.plot),
                   ylim= c(miny, maxy),
                   bty="n",
                   border=FALSE,
                   #ylim = c(0, maxy),
                   #names.arg=if(put.names){row.names(dataset.combined[dataset.combined$genes==i,])}else{NULL},
                   las=2, beside=TRUE,  cex.names=0.6)
        segments(a,df.plot-sd.plot,a,df.plot+sd.plot)
        segments(a-0.5,df.plot+sd.plot,a+0.5,df.plot+sd.plot)
        segments(a-0.5,df.plot-sd.plot,a+0.5,df.plot-sd.plot)
        
      
        
    }
   
    
    res=do.call("rbind.data.frame",tapply(genes,genes,function(x) return(dataset.combined[dataset.combined$genes==x,]),simplify=FALSE))
    return(res)
  } else
  {
    # return data table for this gene
    res=do.call("rbind.data.frame",tapply(genes,genes,function(x) return(dataset.combined[dataset.combined$genes==x,]),simplify=FALSE))
    return(res)
  }
}

par(mar=c(5,4,4,2))
par(oma=c(2,2,2,2))


}
