toBaseTen=function(x,m=128,prec=256,toFile=FALSE,file){
  A=apply(x,2,function(x){
    if (m>64){
    mpfr(sum(2 ^ (which(as.logical(rev(x))) - 1)),precBits=prec)
    }else{
      sum(2 ^ (which(as.logical(rev(x))) - 1))
    }
    return(dat)
  })
  if (toFile == TRUE){
    print("Saving to file...")  
    sapply(dat,function(x) write.table(format(x,scientific=FALSE),file = file,append = TRUE,col.names = F, row.names = F),simplify=TRUE)
    print("File save completed...")
  }  
  dat=0  
}