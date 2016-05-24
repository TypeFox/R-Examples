Bootstrapping.CI <- function (tab, bt, x){


tab <- tab[order(tab$population),]

#Variables
#------------------------------------------------------------------------------------------------------------------------------
# Input:
#           tab <- all.pops.Dest;
#           allelefrequency,sample.sizes <- allelefreq()
#           HWE <- Hardy.Weinberg();

# Passed:
#           tab2,l -> Hardy.Weinberg();
#           tab3 -> allelefreq();
#           allelefrequency,sample.sizes -> Dest.calc();

# Output:
#           Dest.means, Dest.locis, significance.levels -> Workspace;
#           names(tab2) -> screen (print);
#------------------------------------------------------------------------------------------------------------------------------  
                                  
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

          # Dest.Chao values are calculated for each locus as well as
          # averaged over all loci from the empirical data table.

locus.empirical <- DEMEtics.env$values[[1]]
means.empirical <- DEMEtics.env$values[[2]]  
  

tab2 <- split(tab,tab$locus)

          # The table is split according to the several loci.

HWEs <- sapply(tab2,Locuswise.HWE)

                    # For each locus, it has to be found out if all populations are in
                    # Hardy Weinberg Equilibrium.
          
cat("\n","Bootstrapping is carried out ...","\n",sep="")                                      

          # User information.                                                                                        
                                      
          
                                         
repetition <- seq(1:bt)
LocusMeans1 <- lapply(repetition,function(y) locus.means.calc(repetition.val=y,tab2,HWEs,x))
locus <- lapply(LocusMeans1,function(x) x[[1]])
means <- sapply(LocusMeans1,function(x) x[[2]])

assign("means",means,pos = DEMEtics.env)

          # The list "Dest.locus" and the vector "Dest.means" are both saved
          # according to their names in the workspace of R in order to be
          # available for further calculations.

lower.difference <- abs(mean(as.numeric(as.vector(means)),na.rm=TRUE)-quantile(as.numeric(as.vector(means)),.025,na.rm=TRUE))
upper.difference <- abs(mean(as.numeric(as.vector(means)),na.rm=TRUE)-quantile(as.numeric(as.vector(means)),.975,na.rm=TRUE))

critical.values.means <- cbind(means.empirical-lower.difference,means.empirical+upper.difference)


colnames(critical.values.means)<- c("0.95.conf.int.lower","0.95.conf.int.upper")

          # The percentile bootstrap confidence intervals are calculated.
          
loci <- do.call(rbind,locus)

          # The Dest values for the several loci are combined in a single
          # data frame.
          
assign("loci",loci,pos = DEMEtics.env) 

          # The data.frame "Dest.loci" is saved with its according name
          # in the workspace of R in order to be available for further calculations.         
          
loci2<-split(loci,loci$locus)     

          # This data frame is splitted so that the data that belong to the same
          # locus are separated from those that belong to a different locus.
locus.empirical2 <- split(locus.empirical,locus.empirical$locus)
differences1 <- mapply(critical.ci.values,loci2, locus.empirical2,SIMPLIFY=FALSE)
critical.value.loci <- do.call(rbind,differences1)

critical.value.loci<-as.data.frame(critical.value.loci)
colnames(critical.value.loci)<-c("0.95.conf.int.lower","0.95.conf.int.upper")


confidence.limits<-list(critical.value.loci,critical.values.means)
names(confidence.limits)<-c("confidence.limits.per.locus","confidence.limit.over.all.loci")
invisible(confidence.limits)
assign("confidence.limits",confidence.limits,pos = DEMEtics.env)

          # The end result of the bootstrapping is the list called 
          # 'confidence.limits' that contains the significance levels
          # per locus and over all loci.
          # This list is printed and assigned in the workspace  to be available 
          # for further calculations by the name 'confidence.limits'.
}




