YuGene = function(data.prop, progressBar = TRUE){
    
    #KA dded this bit to deal with negative values
    if(sum(data.prop < 0) >=1){
        warning('some negative values were set to 0')
        data.prop[data.prop < 0] = 0
    }
    # ---- end of addition
    
  CumProp= data.frame(matrix(nrow = nrow(data.prop), ncol = ncol(data.prop)) )
  rownames(CumProp) = rownames(data.prop)
  colnames(CumProp) = colnames(data.prop)
  if(progressBar == TRUE) pb <- txtProgressBar(style=3);
   
  for (col in 1:ncol(data.prop)) {                # for each sample (col)
	sample = data.prop[,col]    
	ind=duplicated(sample)
	b=sort(sample,index.return=TRUE,decreasing=TRUE)
	indice=!ind[b$ix] # cumprop on this indice set
	m = cbind(1:length(sample),sample,2)  
	byOrder = m[b$ix,] # sort from highest expression to lowest
	sampleTotal = sum(sample)
	cs = cumsum(byOrder[,2])/sampleTotal  # get cummulative sum of expression on the indice set
	byOrder[indice,3] = cs[indice]  # replace with cumulative values
	for(i in which(!indice))  byOrder[i,3]=byOrder[i-1,3]
	byIndex = byOrder[order(byOrder[,1]),] # re-sort by index

    CumProp[,col] = byIndex[,3] # replace with cum prop values
    if(progressBar == TRUE)
    {
        setTxtProgressBar(pb,col/ncol(data.prop)) # show progress
    }
  }
  cat("\n")# go to next line
  out=1-as.matrix(CumProp,nrow = nrow(data.prop), ncol = ncol(data.prop))
  class(out)="YuGene"
  invisible(out)
} # end function