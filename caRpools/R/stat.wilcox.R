stat.wilcox=function(untreated.list=list(NULL, NULL),treated.list=list(NULL, NULL), namecolumn=1, fullmatchcolumn=2,normalize=TRUE,norm.fun=median, extractpattern=expression("^(.+?)_.+"), controls=NULL, control.picks=300, sorting=TRUE){
  # dataset can be list of designs with readcount
  # fullmatchreadcount must be in fullmatchcolumn
  ###### create dataframe with all datasets
  # get length of each dataset and make sure they are of same length
  #untreated.list= list(ctrl1,ctrl2)
  #treated.list=list(trail1,trail2)
  #print("Start data generation")
  non=lapply(
        c(untreated.list,treated.list),
        function(x) stopifnot(
                        identical(
                          x[,namecolumn],
                          treated.list[[1]][,namecolumn]
                          )
                      )
        )
  #print("get gene names")
  
  gene.names = sub(extractpattern,"\\1",treated.list[[1]][,namecolumn],perl=TRUE)  
  designs=treated.list[[1]][,namecolumn]
  untreated.list=apply(
                    do.call(
                      "cbind",
                      lapply(
                        untreated.list,
                        function(x)
                        {
                            x[,fullmatchcolumn] = apply(x, 1, function(v)
                              {
                              
                              if(as.numeric(v[fullmatchcolumn]) == 0 || is.na(v[fullmatchcolumn]))
                              {
                                return(as.numeric(1))
                              }
                              else
                              {
                                return(as.numeric(v[fullmatchcolumn]))
                              }
                            })
                          
                            if(normalize){ return(x[,fullmatchcolumn]/norm.fun(x[,fullmatchcolumn])+1)}else{return(x[,fullmatchcolumn])}
                          
                          
                        }
                        )
                      ),
                      1,
                      mean)
  #print("Untreated finished")
  treated.list=apply(
                  do.call(
                    "cbind",
                    lapply(
                      treated.list,
                      function(x){
                        
                        x[,fullmatchcolumn] = apply(x, 1, function(v)
                        {
                          if(as.numeric(v[fullmatchcolumn]) == 0 || is.na(v[fullmatchcolumn]))
                          {
                            
                            return(as.numeric(1))
                          }
                          else
                          {
                           
                            return(as.numeric(v[fullmatchcolumn]))
                          }
                        })
                        
                      if(normalize){ return(x[,fullmatchcolumn]/norm.fun(x[,fullmatchcolumn])+1)}else{ return(x[,fullmatchcolumn])}
                      }
                      )
                    ),
                    1,
                    mean
                  )
  #print("Treated finished")
  # Generate new data frame
  # 1.  UNTREATED meaned
  # 2.  TREATED meaned
  # 3.  designs
  # 4.  gene name
  # 5.  foldchange
  # 6.  wilcox p-value
  # get gene names for the dataframe

  dataset.combined <- data.frame( 
    untreated = as.numeric(untreated.list),
    treated = as.numeric(treated.list),
    genes=gene.names,
    foldchange=as.numeric(treated.list)/as.numeric(untreated.list),
    row.names = designs, 
    stringsAsFactors=FALSE)
  
  # remove data where entry is NaN, -Inf or +Inf
  # this will lead to NaN in future calculations
  #print("remove unwanted chars")
  for(i in 1: nrow(dataset.combined))
  {
    # check if either TREATED or UNTREATED Readcount is either NaN, NA, +Inf or -Inf
    # if this is the case, delete the line
    if(!is.finite(dataset.combined[i,4])){
      # set fold change to 1
      dataset.combined[i,4] = 1
    }else if(dataset.combined[i,4] ==0 ){
      dataset.combined[i,4] = 1
    }
  }
  
  ### Wilcoxon Rank sum test on DSS
  ## compare distribution of all sgRNAs for a gene to the group of negative controls (wilcox approach)
  
  # IF no non-targeting controls are set (controls = NULL), all sgRNAs will be used as control reference
  
  # make log2 foldchange
 # print("Get Controls")
  # get neg controls
  if(!is.null(controls))
  {
    # we take only the non-targeting sgrnas
    control.test = dataset.combined$foldchange[dataset.combined$genes==controls]

  }
  else
  {
    # non-targeting controls are NULL, so we take all sgRNA data in the dataset OR control.picks number of sgRNAs
    random.picked = dataset.combined[sample(nrow(dataset.combined), control.picks),"foldchange"]
    control.test = random.picked
    #control.test = dataset.combined$foldchange
  }
  
  #print("Do p-val calculation")
  # Windows machine?
  
    pvals=do.call(
      "rbind.data.frame",
      lapply(
        split(dataset.combined ,
               f = dataset.combined$genes ),
        FUN = function(x) 
        {
          #print("list")
          c(mean(x$untreated),
            mean(x$treated),
            mean(x$foldchange),
            wilcox.test(x$foldchange,
                        control.test,alternative = "t")$p.value
          )
        }
          
      )
    )
    
 #str(pvals)
 #print(pvals[1:10,])
#    # print("Parallel")
#     library("parallel")
#     pvals=do.call(
#       "rbind.data.frame",
#       mclapply(
#         split( dataset.combined ,
#                f = dataset.combined$genes ),
#         FUN = function(x) 
#         {
#         #  print("list")
#           c(mean(x$untreated),
#           mean(x$treated),
#           mean(x$foldchange),
#           wilcox.test(x$foldchange,
#                       control.test,alternative = "t")$p.value
#         )
#         }
#         , mc.cores = detectCores()/
#       )
#     )
  
  
  names(pvals)=c(
                "untreated",
                "treated",
                "foldchange",
                "p.value"
                 )
  #print("correct for multiple testing")
  # correct for multiple testing
  pvals$p.value = p.adjust(pvals$p.value, method = "BH", n = length(pvals$p.value))
  
 # print("Sorting")
  # Sort data frame according to pValue
  if(sorting)
  {
    return(pvals[order(pvals$p.value),])
  }
  else
  {
    return(pvals)
  }
}