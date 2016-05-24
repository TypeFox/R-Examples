calc <- function(allelefrequency,sample.sizes,x){

          # A function that calculates the mean empirical values over all loci and
          # the values (D, Gst) for all loci separately.
          # The arguments of the function:
          # allelefrequency = table with allelefrequencies as given by the function allelefreq().
          # sample.sizes = table with the sample sizes as given by the function allelefreq().
          # x defines whether D, Dest, Gst or Gst.est is calculated

 
          allelefrequency2<-split(allelefrequency,allelefrequency$locus)

                    # The data.frame allelefrequency is splitted according to the loci
                    # that have been examined.
                    
          loci.numbers<-length(allelefrequency2)
          
                    # The number of different loci.

Hj.values <- lapply(allelefrequency2,Calculate.Hj.values)
Hj.values <- do.call(rbind,Hj.values)      

          Hj.values<-as.data.frame(Hj.values)
          colnames(Hj.values)=c("Locus","Population","Hj.value")
          
          
                    # The data are ascribed to a data frame that is called Hj.values.

          Hj.values2<-split(Hj.values,Hj.values$Locus)
          
                    # The data.frame Hj.values is split in order to separate the data for
                    # the different loci.




Hs.one.locus <- sapply(Hj.values2,Hs.value.calculation)
Hs.values <- cbind(Hs.one.locus,names(Hj.values2))

          Hs.values<-as.data.frame(Hs.values)
          colnames(Hs.values)=c("Hs.value","locus")
          
                    # The Hs values are combined in a data frame and the columns are
                    # named.
          
          Ht.values<-as.data.frame(Ht(allelefrequency))
          
                    # The Ht.values are calculated.
          
          sample.sizes2<-split(sample.sizes,sample.sizes$locus)
          
                    # The sample size values are split according to the locus they belong
                    # to.


                                       # Replacing the following for-loop with a function


          Hs.values2 <- split(Hs.values,Hs.values$locus)
          Ht.values2 <- split(Ht.values,Ht.values$locus)

                                        # 2012-02-15 FROM here on something is wrong          
one.locus <- mapply(Gendiff.values,Hs.values2,Ht.values2,sample.sizes2,x)
values1 <- cbind(as.numeric(as.vector(one.locus)),names(sample.sizes2))
                                          
          values1<-as.data.frame(values1)
          colnames(values1)=c(x,"locus")
          
          over.loci<-mean(as.numeric(as.vector((values1[,1]))))
          
          values<-list(values1,over.loci)
          names(values)=c("values.for.loci","Mean.value")
          invisible(values)
          assign("values",values,pos = DEMEtics.env)

}
