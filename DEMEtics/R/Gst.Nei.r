Gst.Nei <- function(filename, bias="correct", object=FALSE, format.table=TRUE, pm="pairwise",statistics="CI",bt=1000){
    
# pm = pairwise/overall
# statistics = p,CI,none (for "all" warn, that it will take long)
# bt = 1000 (number of bootstrapping steps)

filename <- filename
  
if (bias=="correct"){x="Gst.est"
                       }else{x="Gst"}

if (object==TRUE){

                  tab <- get(filename)
                  
          # The data table can be either an object in the R Workspace (a data table
          # that is already loaded in the Workspace)...
                           
                  }else{
                        tab <- read.table(filename, header=TRUE, sep="")}

          # Or a data table that is saved in a .txt-file and that has to be 
          # assigned to the workspace.

if (format.table==TRUE){
                        
                        inputformat(filename,object) 
                        tab <- read.table("Output-Inputformat.txt", header=TRUE, sep="")
                        
                        }
                       
          # If the argument 'format.table is set as true, the format of the
          # table is changed in that format, that is needed for further 
          # calculations.

          # This function will calculate the D values for every locus and
          # the mean D value over all loci for all the populations that 
          # have been sampled.
          
          # The table 'tab' has to be of the following format:
          
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
          
allelefreq(tab)

          # From the table that contains the empirical data, the allelefrequencies
          # as well as the sample.sizes are calculated for each locus.
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are assigned to the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes.         
actual.date <- as.Date(Sys.time())


filename.for.sample.sizes <- paste("sample.sizes.",actual.date,sep="")
filename.for.sample.sizes <- paste(filename.for.sample.sizes,".txt",sep="")
filename.for.sample.sizes <- paste(filename,".",filename.for.sample.sizes,sep="") 

          # The filename under which the table that contains the sample sizes
          # for the several loci, is created. It will encompass the actual
          # date at the end.

filename.for.allelefrequencies <- paste("allelefreq.",actual.date,sep="")
filename.for.allelefrequencies <- paste(filename.for.allelefrequencies,".txt",sep="")
filename.for.allelefrequencies <- paste(filename,".",filename.for.allelefrequencies,sep="")

          # The filename under which the table that contains the allelefrequencies
          # for the several loci, is created. It will encompass the actual
          # date at the end.

write.table(as.data.frame(as.matrix(DEMEtics.env$List[[1]])),file=filename.for.allelefrequencies, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(DEMEtics.env$List[[2]])),file=filename.for.sample.sizes, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The tables are saved in the working directory. Two or more analysis
          # that are carried under the same working directory, are saved in the
          # under the same file name. The data are combined in the order as they
          # have been analysed. No data is lost due to be overwritten.
cat("\n","====================================================================================================","\n",sep="")
cat("====================================================================================================","\n")
cat("======================================== Start of Analysis =========================================","\n")
cat("====================================================================================================","\n")
cat("====================================================================================================","\n")       

cat("\n","Table with allelefrequencies is saved in ","'",filename.for.allelefrequencies,"'","\n",sep="")
cat("\n","Per locus calculated sample sizes are saved in ","'",filename.for.sample.sizes,"'","\n",sep="")

H.out(tab,data.name=FALSE,filename)

if (pm=="pairwise"){pair.pops(tab,statistics,bt,x,filename)
}else {over.all.pops(tab,statistics,bt,x,filename)}



}
