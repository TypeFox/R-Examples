over.all.pops <- function (tab,statistics,bt,x,filename){

  filename <- filename
# x defines whether D, Dest, Gst or Gst.est is calculated
# tab is the table
# statistics either has the value p, CI or both
# bt defines the number of bootstraps to be carried out


allelefreq(tab)

          # From the table that contains the empirical data, the allelefrequencies
          # as well as the sample.sizes are calculated for each locus.
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are assigned to the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes.         
          

calc(DEMEtics.env$allelefrequency,DEMEtics.env$sample.sizes,x)

          # The D or Gst values for the several loci and the mean D or Gst value of 
          # all loci are calculated.
          # The results are available from the object 'values'.


loci=as.data.frame(DEMEtics.env$values[[1]])
colnames(loci)=c(paste(x,".locus",sep=""),"Locus")

          # A table is created that contains the D or Gst values for the several
          # loci 

mean=as.data.frame(DEMEtics.env$values[[2]])
colnames(mean)=paste(x,".mean",sep="")

          # A table is created that contains the mean D or Gst value over all
          # loci

all.pops=list(loci,mean)
names(all.pops)=c(paste(x,".loci.over.all.populations",sep=""),paste(x,".mean.over.all.populations",sep=""))

          # These two tables are combined in a single list.

assign("all.pops",all.pops,pos = DEMEtics.env)

          # The object 'all.pops', a list, is assigned to the R workspace.
cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("++++++++++++++++++++++++++++++++++++ RESULTS WITHOUT STATISTICS ++++++++++++++++++++++++++++++++++++","\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")      

print(all.pops)

          # The result tables are printed.

actual.date <- as.Date(Sys.time())

          # The actual date is determined.

filename.for.loci <- paste("overall.",x,".loci.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

          # The filename under which the table that contains the D or Gst values
          # for all the loci seperately, is created. It will encompass the actual
          # date at the end.

filename.for.mean <- paste("overall.",x,".mean.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

          # The filename under which the table that contains the mean D or Gst value
          # over all the loci, is created. It will encompass the actual
          # date at the end.
          # date at the end.

write.table(as.data.frame(as.matrix(all.pops[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(all.pops[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The tables are saved in the working directory. Two or more analysis
          # that are carried under the same working directory, are saved in the
          # under the same file name. The data are combined in the order as they
          # have been analysed. No data is lost due to be overwritten.
cat("----------------------------------------------------------------------------------------------------","\n")
cat("The ",x,".mean value was calculated as the arithmetic mean of the ",x,".loci values","\n",sep="")
cat("\n",x," values for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n",x," values over all loci are saved in ","\n","'",filename.for.mean,"'","\n",sep="")


################# Bootstrapping part




tab.pops <- split(tab,tab$population)

          # The table is splitted up into the data that belong to different
          # populations.           

cat("\n","====================================================================================================","\n",sep="")
cat("=============================== Bootstrapping for P-value Calculation ==============================","\n")
cat("====================================================================================================","\n","\n")
      


          # User information.
          
print(paste("Start of analysis: ",Sys.time(),sep=""))

          # The actual time is printed.

time1 <- Sys.time()

          # The time before the bootstrapping.

if (statistics=="p"||statistics=="all"){
cat("\n","WARNING: Depending on the size of your input data, the performance of your computer and the number of ","\n","bootstrap resamplings you have chosen, bootstrapping can take very long (hours to days).","\n",sep="")

loci <- Bootstrapping.p(tab,bt,x)

          # The bootstrap values are available from the objects 'loci'
          # and 'means'.    
          
time2 <- Sys.time()

          # The time after the bootstrapping.

bootstrap.time <- difftime(time2,time1,units="secs") 

          # Duration of the Bootstrapping analysis (in seconds).
            
bootstrap.time.output <- round((as.numeric(bootstrap.time/60)),2) 

          # Duration of the Bootstrapping analysis (in minutes), rounded (2 digits).

cat("\n","... Bootstrapping for p-value calculation is terminated.","\n",sep=" ")
print(paste("Duration of the bootstrapping analysis (min):",bootstrap.time.output,sep=" "))

          # User information
allelefreq(tab)

          # From the table that contains the empirical data, the allelefrequencies
          # as well as the sample.sizes are calculated for each locus.
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are assigned to the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes.         
          

calc(DEMEtics.env$allelefrequency,DEMEtics.env$sample.sizes,x)

          # The D or Gst values for the several loci and the mean D or Gst value of 
          # all loci are calculated.
          # The results are available from the object 'values'.
          

                 
          # From the bootstrapping data, the p-values are assigned to the
          # empirical D or Gst values.
          
loci2 <- split(loci,loci$locus)

          # The bootstrap values are separated in different tables according
          # to the loci they belong to.
values2 <- split(DEMEtics.env$values[[1]],DEMEtics.env$values[[1]]$locus)

p.values.1 <- mapply(p.values.loci.calc,loci2,values2,SIMPLIFY=FALSE)
p.values.loci <- do.call(c,p.values.1)

assign("p.values.loci",p.values.loci,pos = DEMEtics.env)
                          
bootstrapped.values=DEMEtics.env$means

          # The mean D or Gst values over all loci that have been obtained from
          # the bootstrapping.

empirical.value=DEMEtics.env$values[[2]]

          # The empirical mean D value over all loci. 

p.value.over.all=p.val(empirical.value,bootstrapped.values)

          # This function returns the p-value for the actual empirical D or Gst value
          # in the object 'p.value'.
          
p.value.over.all=round(DEMEtics.env$p.value,4)

          # This p.value is rounded up to 4 decimal places.

loci=cbind(DEMEtics.env$values[[1]],as.vector(as.numeric(p.values.loci)))
loci=as.data.frame(loci)
colnames(loci)=c(paste(x,".locus",sep=""),"Locus","P.value")

          # A table is created that contains the D or Gst values for the several
          # loci with the according p-values obtained
          # from the bootstrapping analysis.

mean=cbind(DEMEtics.env$values[[2]],p.value.over.all)
mean=as.data.frame(mean)
colnames(mean)=c(paste(x,".mean",sep=""),"P.value")

          # A table is created that contains the mean D or Gst value
          # over all loci with the according p-values obtained from the
          # bootstrapping analysis.

all.pops=list(loci,mean)
names(all.pops)=c(paste(x,".loci.over.all.populations",sep=""),paste(x,".mean.over.all.populations",sep=""))

          # These two tables are combined in a single list.

assign("all.pops",all.pops,pos = DEMEtics.env)

          # The object 'all.pops', a list, is assigned to the R workspace.

cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("+++++++++++++++++++++++++++++++++++++++ RESULTS WITH P-VALUES ++++++++++++++++++++++++++++++++++++++","\n")
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")    


print(all.pops)

          # The result tables are printed.

actual.date <- as.Date(Sys.time())

          # The actual date is determined.

filename.for.loci <- paste("overall.",x,".loci.p.",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

          # The filename under which the table that contains the D or Gst values
          # for all the loci seperately, is created. It will encompass the actual
          # date at the end.

filename.for.mean <- paste("overall.",x,".mean.p.",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

          # The filename under which the table that contains the mean D value
          # over all the loci, is created. It will encompass the actual
          # date at the end.



write.table(as.data.frame(as.matrix(all.pops[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(all.pops[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

cat("----------------------------------------------------------------------------------------------------","\n")         
print(paste("This analysis was finished on",Sys.time()))
cat("The ",x,".mean value was calculated as the arithmetic mean of the ",x,".loci values","\n",sep="")

cat("\n",x, " values over all populations for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n",x, " values averaged over all loci and populations are saved in ","\n","'",filename.for.mean,"'","\n",sep="")

          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.





}







if (statistics=="CI"||statistics=="all"){
cat("\n","====================================================================================================","\n",sep="")
cat("========================= Bootstrapping for Confidence Interval Estimation =========================","\n")
cat("====================================================================================================","\n")    


cat("\n","WARNING: Depending on the size of your input data, the performance of your computer and the number of ","\n","bootstrap resamplings you have chosen, bootstrapping can take very long (hours to days).","\n",sep="")
time1 <- Sys.time()
Bootstrapping.CI(tab,bt,x)

          # The bootstrap values are available from the objects 'loci'
          # and 'means'.    
          
time2 <- Sys.time()

          # The time after the bootstrapping.

bootstrap.time <- difftime(time2,time1,units="secs") 

          # Duration of the Bootstrapping analysis (in seconds).
            
bootstrap.time.output <- round((as.numeric(bootstrap.time/60)),2) 

          # Duration of the Bootstrapping analysis (in minutes), rounded (2 digits).

cat("\n","... Bootstrapping for confidence interval estimation is terminated.","\n",sep=" ")
print(paste("Duration of the bootstrapping analysis (min):",bootstrap.time.output,sep=" "))

          # User information


     
allelefreq(tab)

          # From the table that contains the empirical data, the allelefrequencies
          # as well as the sample.sizes are calculated for each locus.
          # The table that contains the allelefrequencies and the table that
          # lists the sample sizes are assigned to the R workspace in the 
          # object List, but also separately in the object allelefrequency
          # and the object sample.sizes.

calc(DEMEtics.env$allelefrequency,DEMEtics.env$sample.sizes,x)

          # The D values for the several loci and the mean D value of 
          # all loci are calculated.
          # The results are available from the object 'D.values'. 
                 
          # From the bootstrapping data, the p-values are assigned to the
          # empirical D values.
          
loci2 <- split(tab,tab$locus)

          # The bootstrap values are separated in different tables according
          # to the loci they belong to.
          
number.loci <- length(loci2)

          # The number of loci that have been studied.

values2 <- split(DEMEtics.env$values[[1]],DEMEtics.env$values[[1]]$locus)
empirical.value.loci1 <- lapply(values2, function(x) as.numeric(as.vector(x[,1])))
empirical.value.loci <- do.call(c,empirical.value.loci1)

# The empirical D or Gst values for the several loci are assigned to a
# vector.                          
                                                   


empirical.value=DEMEtics.env$values[[2]]

          # The empirical mean D value over all loci. 


loci=cbind(DEMEtics.env$values[[1]],DEMEtics.env$confidence.limits[[1]])
loci=as.data.frame(loci)
colnames(loci)=c(paste(x,".locus",sep=""),"Locus","Lower.0.95.CI","Upper.0.95.CI")

          # A table is created that contains the D or Gst values for the several
          # loci 

mean=cbind(DEMEtics.env$values[[2]],DEMEtics.env$confidence.limits[[2]])
mean=as.data.frame(mean)
colnames(mean)=c(paste(x,".mean",sep=""),"Lower.0.95.CI","Upper.0.95.CI")

          # A table is created that contains the mean D or Gst value over all
          # loci with the according confidence levels and the p-value obtained
          # from the bootstrapping analysis.

all.pops <- list(loci,mean)
names(all.pops) <- c(paste(x,".loci.over.all.populations",sep=""),paste(x,".mean.over.all.populations",sep=""))

          # These two tables are combined in a single list.

assign("all.pops",all.pops,pos = DEMEtics.env)

          # The object 'all.pops', a list, is assigned to the R workspace.

cat("\n","++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n",sep="")
cat("+++++++++++++++++++++++++++++ RESULTS WITH CONFIDENCE INTERVAL LIMITS ++++++++++++++++++++++++++++++","\n")  
cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++","\n","\n")    

print(all.pops)

          # The result tables are printed.

actual.date <- as.Date(Sys.time())

          # The actual date is determined.

filename.for.loci <- paste("overall.",x,".loci.ci",actual.date,sep="")
filename.for.loci <- paste(filename.for.loci,".txt",sep="")
filename.for.loci <- paste(filename,".",filename.for.loci,sep="")

          # The filename under which the table that contains the D values
          # for all the loci seperately, is created. It will encompass the actual
          # date at the end.

filename.for.mean <- paste("overall.",x,".mean.ci",actual.date,sep="")
filename.for.mean <- paste(filename.for.mean,".txt",sep="")
filename.for.mean <- paste(filename,".",filename.for.mean,sep="")

          # The filename under which the table that contains the mean D value
          # over all the loci, is created. It will encompass the actual
          # date at the end.

write.table(as.data.frame(as.matrix(all.pops[[1]])),file=filename.for.loci, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(as.data.frame(as.matrix(all.pops[[2]])),file=filename.for.mean, append = TRUE, quote = FALSE, sep = " ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

          # The tables are saved in the working directory. Two or more analysis
          # that are carried under the same working directory, are saved in the
          # under the same file name. The data are combined in the order as they
          # have been analysed. No data is lost due to be overwritten.

cat("----------------------------------------------------------------------------------------------------","\n")         
print(paste("This analysis was finished on",Sys.time()))
cat("The ",x,".mean value was calculated as the arithmetic mean of the ",x,".loci values","\n",sep="")

cat("\n",x, " values over all populations for every locus are saved in ","\n","'",filename.for.loci,"'","\n",sep="")
cat("\n", x, " values averaged over all loci and populations saved in ","\n","'",filename.for.mean,"'","\n",sep="")


          # User information about the end date of the analysis and the filenames
          # under which the several tables have been saved.



}
  

cat("\n","====================================================================================================","\n",sep="")
cat("====================================================================================================","\n")            
cat("========================================= End of Analysis ==========================================","\n")   
cat("====================================================================================================","\n")
cat("====================================================================================================","\n")    
}
