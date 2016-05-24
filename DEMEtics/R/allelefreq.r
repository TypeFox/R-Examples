
allelefreq <- function(tab){

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           tab <- Dest.Chao(), all.pops.Dest(), all.pops.Gst();

# Output:
#           sample.sizes, allelefrequency, List -> Workspace;
#------------------------------------------------------------------------------------------------------------------------------
  
          # A function that calculates the allelefrequencies and sample sizes per
          # locus from a table of the following format:
          
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
     
y <- tab

          # the argument, a table,  is assigned to the object 'y'.

y2 <- split(y,y$locus)

          # The data are splitted according to the loci that have been
          # examined.
sample.sizes1 <- lapply(y2,sample.size.calc)
sample.sizes <- do.call(rbind,sample.sizes1)
colnames(sample.sizes)=c("population","sample.size","locus")

allelefrequency1 <- lapply(y2,allelefreq.calcs)
allelefrequency <- do.call(rbind,allelefrequency1)
colnames(allelefrequency)= c("allele","number","population","locus","proportion")
rownames(allelefrequency)=seq(1:length(allelefrequency[,1]))   
         
List<-list(allelefrequency,sample.sizes)
names(List)=c("Table","Sample sizes")
invisible(List)

          # The matrix 'allelefrequency' is combined with the matrix 'sample.sizes'
          # in a list.
          
assign("sample.sizes",sample.sizes,pos = DEMEtics.env)
assign("allelefrequency",allelefrequency,pos = DEMEtics.env)
assign("List",List,pos = DEMEtics.env)

          # The objects 'sample.sizes' and 'allelefrequency'are ascribed to the
          # workspace in order to be available for further calculations.

}
