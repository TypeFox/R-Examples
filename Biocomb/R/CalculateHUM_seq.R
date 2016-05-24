CalculateHUM_seq<-function(data,indexF,indexClass,indexLabel)
{
  indexL=NULL

  label=levels(data[,indexClass])
  if(is.null(label)) return

  for(i in 1:length(indexLabel))
  {
  indexL=c(indexL,which(label==indexLabel[i]))
  }


  output=NULL
  seqMax=NULL

  indexEach=NULL
  indexUnion=NULL
  for(i in 1:length(label))
  {
      vrem=which(data[,indexClass]==label[i])
      indexEach=c(indexEach,list(vrem))
      if(label[i]%in%indexLabel)
        indexUnion=union(indexUnion,vrem)
  }


  len=length(indexL)
  seq=gtools::permutations(len,len,1:len)

  for(i in 1:length(indexF))
  {
    s_data=NULL
    dataV=data[,indexF[i]]
    prodValue=1
    for (j in 1:length(indexLabel))
    {
    vrem=sort(dataV[indexEach[[indexL[j]]]])

    s_data=c(s_data,list(vrem))
    prodValue = prodValue*length(vrem)
    }
    #seq=sort(d_median,index.return=TRUE)

    #out=CalcGene(s_data,seq$ix, prodValue)

    #calculate the threshold values for plot of 2D ROC and 3D ROC
    thresholds <- sort(unique(dataV[indexUnion]))
    thresholds=(c(-Inf, thresholds) + c(thresholds, +Inf))/2

    out=CalcGene(s_data,seq,prodValue,thresholds)

    output=c(output,out$HUM)
    seqMax=cbind(seqMax,out$seq)
  }

  names(output)=names(data[,indexF,drop = FALSE])
  colnames(seqMax)=names(output)

  #return(output)
  return(list(HUM=output,seq=seqMax))
}
