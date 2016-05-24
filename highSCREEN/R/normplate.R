normplate = function(mainplate, datbefore, datafter, cmap, plate, triplicate, norm="bscore", poscont=NULL, negcont=NULL){
  cnamesb = c(as.character(cmap[,1]), as.character(cmap[,2]))
  cnamesa = c(as.character(cmap[,1]), as.character(cmap[,2]))

  if (!is.null(poscont)){
    posb = mean(unlist(datbefore[,c(1,ncol(datbefore))])[grep(poscont, cnamesb)])
    posa = mean(unlist(datafter[,c(1,ncol(datafter))])[grep(poscont, cnamesa)])
  }

  if (!is.null(negcont)){
    negb = mean(unlist(datbefore[,c(1,ncol(datbefore))])[grep(negcont, cnamesb)])
    nega = mean(unlist(datafter[,c(1,ncol(datafter))])[grep(negcont, cnamesa)])
  }

  if (norm == "bscore"){
    datbefore[,-c(1,ncol(datbefore))] = medpolish(datbefore[,-c(1,ncol(datbefore))])[["residuals"]]/mad(unlist(datbefore[,-c(1,ncol(datbefore))]))
    datafter[,-c(1,ncol(datafter))] = medpolish(datafter[,-c(1,ncol(datafter))])[["residuals"]]/mad(unlist(datafter[,-c(1,ncol(datafter))]))
  }
  else if (norm == "zscore"){
    datbefore[,-c(1,ncol(datbefore))] = (datbefore[,-c(1,ncol(datbefore))]-mean(unlist(datbefore[,-c(1,ncol(datbefore))])))/sd(unlist(datbefore[,-c(1,ncol(datbefore))]))
    datafter[,-c(1,ncol(datafter))] = (datafter[,-c(1,ncol(datafter))]-mean(unlist(datafter[,-c(1,ncol(datafter))])))/sd(unlist(datafter[,-c(1,ncol(datafter))]))
  }
  else if (norm == "cscore"){
    datbefore[,-c(1,ncol(datbefore))] = (datbefore[,-c(1,ncol(datbefore))]-negb)/(posb-negb)*100
    datafter[,-c(1,ncol(datafter))] = (datafter[,-c(1,ncol(datafter))]-nega)/(posa-nega)*100
  }
  else if (norm == "raw"){
    print("raw")
  }
  else
    stop ("unknown normalization.")

  dat = rbind(datbefore, datafter)

  dat = as.vector(t(dat))
  dat = data.frame(Time=c(rep("Before", length(dat)/2),rep("After", length(dat)/2)), Plate=plate, Triplicate=triplicate, Norm=norm, well=rep(c(as.vector(sapply(LETTERS[seq(1:8)],function(x){sapply(seq(1,12), function(y){paste(x,y,sep="")})}))),2), row=rep(c(as.vector(sapply(LETTERS[seq(1:8)],function(x){rep(x,12)}))),2), col=rep(c(as.vector(sapply(seq(1:8),function(x){seq(1,12)}))),2), score=c(dat[1:(length(dat)/2)], dat[(length(dat)/2+1):length(dat)]))

  MainPlate = rep(mainplate, nrow(dat))
  welltype = rep("Compound", nrow(dat))
  welltype[seq(1,nrow(dat),12)] = c(as.character(cmap[,1]), as.character(cmap[,1]))
  welltype[seq(12,nrow(dat),12)] = c(as.character(cmap[,2]), as.character(cmap[,2]))
  dat = data.frame(MainPlate, dat, welltype)

  return(dat)
}
