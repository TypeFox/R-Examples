H.out <- function(tab,data.name=FALSE,filename){

  filename <- filename
          # A function that calculates the Hs Hs_est, Ht and Ht_est values for each locus separately.
          # The arguments of the function:

          # If data.name=TRUE, the name of the data file the result
          # table is saved at, can be exactly determined.  In this
          # case v is thre front part of the name and h the hind
          # part. The mean part will always be: H.values, so that the
          # file name will be as follows: v.H.values.h.txt. The parts
          # v and h have to be quoted, like v="FrontPart",
          # h="HindPart".

tab.pops <- split(tab,tab$population)

          # The table is splitted up into the data that belong to different
          # populations.           

allelefreq(tab)

          # From the table that contains the empirical data, the allelefrequencies
          # as well as the sample.sizes are calculated for each locus.
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are assigned to the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes.
  allelefrequency <- DEMEtics.env$allelefrequency

          allelefrequency2<-split(DEMEtics.env$allelefrequency,allelefrequency$locus)

                    # The data.frame allelefrequency is splitted according to the loci
                    # that have been examined.
                    


Hj.values1 <- lapply(allelefrequency2,Hj.values.calculation)
Hj.values <- do.call(rbind,Hj.values1)
Hj.values<-as.data.frame(Hj.values)
colnames(Hj.values)=c("Locus","Population","Hj.value")


  # Calculating non-biased Hj values according to Nei (1978)
sample.sizes <- DEMEtics.env$sample.sizes
Hj.values.sample.sizes <- cbind(Hj.values,sample.size=subset(sample.sizes$sample.size,sample.sizes$population%in%Hj.values$Population & sample.sizes$locus%in%Hj.values$Locus,drop=TRUE))

Hj.est.values <- apply(Hj.values.sample.sizes,1,function(x){Hj=as.numeric(x[3])
                                                                N=as.numeric(x[4])
                                                                (2*N/(2*N-1))*Hj})
Hj.est.table <- as.data.frame(cbind(Locus=as.character(Hj.values$Locus),Population=as.character(Hj.values$Population),Hj.est=as.numeric(Hj.est.values)))
Hj.means <- tapply(as.vector(as.numeric(as.vector(Hj.est.table$Hj.est))),as.character(Hj.est.table$Population),mean)
Hj.est.means <- as.data.frame(cbind(Population=as.character(rownames(Hj.means)),Hj.est.mean=as.vector(as.numeric(Hj.means))))
          


Hj.values2<-split(Hj.values,Hj.values$Locus)
          
                    # The data.frame Hj.values is split in order to separate the data for
                    # the different loci.


Hs.values1 <- lapply(Hj.values2,Hs.values.calculation)
Hs.values <- do.call(rbind,Hs.values1)

  Hs.values<-as.data.frame(Hs.values)

  colnames(Hs.values)=c("Hs.value","locus")
          
                    # The Hs values are combined in a data frame and the columns are
                    # named.
          
          Ht.values<-as.data.frame(Ht(DEMEtics.env$allelefrequency))

                       
          sample.sizes2<-split(sample.sizes,sample.sizes$locus)
          
                    # The sample size values are split according to the locus they belong
                    # to.

          # To be able to calculate Hs_est and Ht_est values per locus, some functions and sizes first have to be defined:



 Hs.values2 <- split(Hs.values,Hs.values$locus)
 Ht.values2 <- split(Ht.values,Ht.values$locus)


H.est.values1 <- mapply(Hs.Ht.est.calculation,Hs.values2,Ht.values2,sample.sizes2,SIMPLIFY=FALSE)
Hs.est.values1 <- lapply(H.est.values1,function(x) x[[1]])
Hs.est.values <- do.call(rbind,Hs.est.values1)

Ht.est.values1 <- lapply(H.est.values1,function(x) x[[2]])
Ht.est.values <- do.call(rbind,Ht.est.values1)
 


          Hs.est.values<-as.data.frame(Hs.est.values)
          colnames(Hs.est.values)=c("Hs.est.value","locus")

          Ht.est.values<-as.data.frame(Ht.est.values)
          colnames(Ht.est.values)=c("Ht.est.value","locus")

          H.output <- cbind(as.character(as.vector(Ht.est.values$locus)),as.numeric(as.vector(Hs.values$Hs.value)),as.numeric(as.vector(Hs.est.values$Hs.est.value)),as.numeric(as.vector(Ht.values$Ht.value)),as.numeric(as.vector(Ht.est.values$Ht.est.value)))

          colnames(H.output)=c("locus","Hs","Hs.est","Ht","Ht.est")

          H.output <- as.data.frame(H.output)
  
cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("+++++++++++++++++++++++++++++++++++++++++ HETEROZYGOSITIES +++++++++++++++++++++++++++++++++++++++++","\n") 
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")
  
print(H.output)

if (data.name==TRUE){

filename.output <- paste(DEMEtics.env$v,"_H.values_",sep="")
filename.output <- paste(filename.output,DEMEtics.env$h,sep="")
filename.output <- paste(filename.output,".txt",sep="")

filename2.output <- paste(DEMEtics.env$v,"_Hest_loci.values_",sep="")
filename2.output <- paste(filename2.output,DEMEtics.env$h,sep="")
filename2.output <- paste(filename2.output,".txt",sep="")

filename3.output <- paste(DEMEtics.env$v,"_Hest_mean.values_",sep="")
filename3.output <- paste(filename3.output,DEMEtics.env$h,sep="")
filename3.output <- paste(filename3.output,".txt",sep="")


}else{

filename.output <- paste("H.values",".txt",sep="")
filename.output <- paste(filename,".",filename.output,sep="")

  filename2.output <- paste("Hest_loci.values",".txt",sep="")
filename2.output <- paste(filename,".",filename2.output,sep="")

filename3.output <- paste("Hest_mean.values",".txt",sep="")
filename3.output <- paste(filename,".",filename3.output,sep="")}

          # The filename under which the table that contains the Hs,
          # Hs.est, Ht and Ht.est values for all the loci seperately,
          # is created

write.table(as.data.frame(as.matrix(H.output)),file=filename.output, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

write.table(as.data.frame(as.matrix(Hj.est.table)),file=filename2.output, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

write.table(as.data.frame(as.matrix(Hj.est.means)),file=filename3.output, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The table is saved in the working directory. Two or more
          # analysis that are carried under the same working
          # directory, are saved under the same file name. The data
          # are combined in the order as they have been analysed. No
          # data is lost due to be overwritten.

cat("----------------------------------------------------------------------------------------------------","\n")         
cat("\n","Hs, Hs.est, Ht and Ht.est values for each locus are saved in ","'",filename.output,"'","\n",sep="")
cat("\n","Bias corrected Heterozygosities (according to Nei (1978)) per population and locus are saved in ","'",filename2.output,"'","\n",sep="")
cat("\n","Bias corrected Heterozygosities (according to Nei (1978)) per population averaged over all loci are saved in ","'",filename3.output,"'","\n",sep="")


          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.

}
