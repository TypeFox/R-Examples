ibdhap.reduce.states <-
function(qibd.filename, dat.filename, output.filename){


   #glean information from the par file
   par.file<-read.table(dat.filename, fill =TRUE, colClasses = "numeric")

   n.snps<- par.file[3,1]  #number of markers
   n.sets<-par.file[2,1]   #number of set of haplotypes
   
   temp<-FALSE
   

   for( iset in 1:n.sets){

      #make the text above the (sets of) marker information in both files the same
      line1<-scan(qibd.filename, skip=((n.snps+2)*(iset-1)), nlines=1, what='raw')
      write(line1, file=output.filename, ncolumns=length(line1), append=temp)
      temp<-TRUE
  
      line2<-scan(qibd.filename, skip=((n.snps+2)*(iset-1)+1), nlines=1, what='raw')
      write(line2, file=output.filename, ncolumns=length(line2), append=TRUE)

      rm(line1,line2)
      
      #read the qibd.filename, just the information for the markers in "iset"
      ibd.dat<-scan(qibd.filename,what = "numeric",skip = (2*(iset)+(n.snps)*(iset-1)), nlines = n.snps)
      ibd.dat<-t(matrix(as.numeric(unlist(ibd.dat)), nrow=17))




     reduced.dat <- t( apply(ibd.dat,1,sumcol) )

     reduced.dat <- cbind( ibd.dat[,1:2], reduced.dat)

     write.table( reduced.dat,file=output.filename,  append = TRUE, col.names=FALSE, row.names=FALSE)
     }
     
     return(paste("file", output.filename, "written in", getwd()) )
}

