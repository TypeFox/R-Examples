Bootstrapping.p <- function (tab,bt,x){
# x defines whether D, Dest, Gst or Gst.est is calculated
# tab is the table
# bt defines the number of bootstraps to be carried out
  
                                  
          # The input for this function has to be a table of the following format:
          
          #     individual population      fragment.length   locus
          # 1        B1.1     Borkum            323          L12
          # 2        B1.1     Borkum            266          L12
          # 3        B1.2     Borkum            325          L12
          # 4        B1.2     Borkum            274          L12
          # 5        B1.3     Borkum            266          L12
          # 6        B1.3     Borkum            323          L12
          # 7        B1.4     Borkum            325          L12 
          
          # The column names must be equal to those in this example. 
          # Certain tables of another format can be converted to this format
          # by the function 'inputformat()' that is included in this package.

          # bt defines the times of bootstrap-resampling.
allelefreq(tab)  
  
Empirical.values <- calc(DEMEtics.env$allelefrequency,DEMEtics.env$sample.sizes,x)

          # D.Chao values are calculated for each locus as well as
          # averaged over all loci from the empirical data table.

locus.empirical <- DEMEtics.env$values[[1]]
means.empirical <- DEMEtics.env$values[[2]]  
  

tab2 <- split(tab,tab$locus)

          # The table is split according to the several loci.
          


                    # For each locus, it has to be found out if all populations are in
                    # Hardy Weinberg Equilibrium.

                           
HWEs <- sapply(tab2,Locuswise.HWE)

          # zeros are transformed to FALSE, ones to TRUE.  
cat("\n","Bootstrapping is carried out ...","\n",sep="")                                      

          # User information.   


                                         
repetition <- seq(1:bt)
LocusMeans1 <- lapply(repetition,function(y) locus.means.calc(repetition.val=y,tab2,HWEs,x))
locus <- lapply(LocusMeans1,function(x) x[[1]])
means <- sapply(LocusMeans1,function(x) x[[2]])

assign("means",means,pos = DEMEtics.env)

          # The list "locus" and the vector "means" are both saved
          # according to their names in the workspace of R in order to be
          # available for further calculations.

loci <- do.call(rbind,locus)

          # The D or Gst values for the several loci are combined in a single
          # data frame.
          
assign("loci",loci,pos = DEMEtics.env) 

          # The data.frame "loci" is saved with its according name
          # in the workspace of R in order to be available for further calculations.         
          
}




