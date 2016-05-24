
gkmsvm_train = function (kernelfn, posfn, negfn, svmfnprfx,  Type="C-svc", C=1, shrinking=FALSE, ...){
#TODO: add bootstrapping and cv capabilities -- also autyomatic choise of C  . check if kernlab does that 
  

  #  library(seqinr)
  #  library(kernlab)
  #  library(utils)
  if (requireNamespace("seqinr", quietly = TRUE)&
      requireNamespace("utils", quietly = TRUE)&
      requireNamespace("kernlab", quietly = TRUE)){
        
      
    
  #  negfn= '/Users/mghandi/gkmsvm/test/testneg9.fa'
  #  posfn= '/Users/mghandi/gkmsvm/test/testpos9.fa'
  #  kernelfn= '/Users/mghandi/gkmsvm/test/test9kernel.txt'
    
    pos = seqinr::read.fasta(posfn)
    npos = length(unique(names(pos))) #length(pos)
    neg = seqinr::read.fasta(negfn)
    nneg = length(unique(names(neg))) #length(neg)
    nseq = npos+nneg; 
    
    if(length(which(duplicated(names(pos))))>0){
      print(paste("Error: duplicated sequence ID in", posfn))
      print(names(pos)[which(duplicated(names(pos)))])
      return;
    }
    if(length(which(duplicated(names(neg))))>0){
      print(paste("Error: duplicated sequence ID in", negfn))
      print(names(neg)[which(duplicated(names(neg)))])
      return;
    }
    
    mat <- data.matrix( utils::read.table(file=kernelfn, fill=TRUE, col.names=paste("V", 1:nseq)))
    mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
    rownames(mat)=colnames(mat)
    K <- kernlab::as.kernelMatrix(mat)
    y = c(rep(1, npos), rep(0, nneg)); names(y)=rownames(mat)
  
  #  svp <- ksvm(K, y, type="C-svc", C=1)
    svp <- kernlab::ksvm(K, y, type=Type, C=C, shrinking=shrinking, ...)
    
    seqnames = c(names(pos), names(neg))
    
    if(svp@nSV>0){
      alpha = unlist(svp@alpha )
      ii = unlist(svp@SVindex)
      jj = which(ii>npos); 
      alpha[jj]= -alpha[jj];
      
      utils::write.table(cbind(seqnames[ii], sprintf("%11.6e",alpha)),
                  file = paste(svmfnprfx, 'svalpha.out', sep='_'),
                  col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\t')
      
      svseqs = c(pos,neg)[ii]; 
      seqinr::write.fasta(svseqs, names(svseqs),  file.out= paste(svmfnprfx, 'svseq.fa', sep='_'))
    }
  }
}
