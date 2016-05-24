inputformat <- function(filename,object)
{

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#          filename,object <- all.pops.Dest.Chao(), all.pops.Dest(), all.pops.Gst(), pair.pops.Gst(),
#                             pair.pops.Dest.Chao(), pair.pops.Dest();

# Output:
#          Output-Inputformat.txt -> working directory (file);
#------------------------------------------------------------------------------------------------------------------------------
  

     cat("The table format is transformed...")
     
      if (object==TRUE){

                  y <- get(filename)
                  
          # The data table can be either an object in the R Workspace (a data table
          # that is already loaded in the Workspace)...
                           
                  }else{
                        y <- read.table(filename, header=TRUE, sep="")
                        } 

          # Or a data table that is saved in a .txt-file and that has to be 
          # assigned to the workspace.      
     
     
     attach(y)                                                                #Column names attached
     on.exit(detach(y))
     names(y)                                                                 # And returned

     c <- split(y, y[,1])                                                     # Matrix is separated into data belonging to different samples
     p <- c                                                                   # Matrix p stays unchanged
     s <- c(length(y)-3)                                                      # Number 3 Digs columns (less one!)
     v <- colnames(y)                                                         # Column names placed in vector v
     w <- length(c)                                                           # Number of samples
     k <- length(y[1,])                                                       # Length of input-table
     a <- y[1,]                                                               # Starting point of output-table

     for (g in 1:w)                                                           # For each sample
     {
          c[[g]] <- rbind(p[[g]], y[1:s,])                                    # For each samle, rows are bound to c
          c[[g]][2:(s+1),1] <- c[[g]][1,1]                                    # 1.Column - values transferred
          c[[g]][2:(s+1),2] <- c[[g]][1,2]                                    # 2.Column - values transferred


        for (l in 1:s)
                {
                c[[g]][l+1,3] <- c[[g]][1,l+3]                                # Alleles are transcribed to column 3
                }

        for (f in 1:(s+1))
                {

                 c[[g]][f,4] <- v[f+2]                                        # Column 4 is filled with names of loci for each allele
                }

      }

      for (j in 1:w)
       {
        a <- rbind(a, c[[j]])                                                  # Table is recombined in a
       }
      a <- a[-1,]                                                              # The first row is deleted
      Output <- a[,-c(5:10000)]                                                # Superfluous columns are deleted

      for (f in 1:length(Output[,1]))
       {
        if (f%%2==1) {Output[f+1,4] <- Output[f,4]}
        if (f%%2==0) {}                                                        # Alleles are combined
       }


      colnames(Output) <- c("individual", "population", "fragment.length", "locus")   # Columns renamed
      rownames(Output) <- seq(1:length(Output[,1]))



write.table(Output, file = "Output-Inputformat.txt", append = FALSE, quote = FALSE, sep = " ", na = "NA", dec = ",", row.names = TRUE, col.names = TRUE)


print("Transformed Table saved as 'Output-Inputformat.txt'")
invisible(Output)

}
