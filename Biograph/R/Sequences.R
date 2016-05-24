Sequences <-
function (Bdata,mean_median) {
# TASK: From Bdata, determine number of different sequences
# Output : sequences_sorted 
#
# ntcase = number of different sequences
# case(ntcase) = sequences based on Bdata$path
# ns_case(ntcase) = number of states in pathway type i 
# case_ind[i]   =  path i belongs to pathtype case_ind
#    frequency of case_ind[i] gives number of subjects with pathtype case_ind 
# initiate the following vectors:
#          - pathway in class k
#          - number of sequences in class k
# initiate number of classes (ntcase)
# n_match = Bdata$path matches the n_match-th element of vector case
# mages = mean ages / median ages at transitions in different pathtypes
z <- check.par (Bdata)
if (missing(mean_median)) mean_median <- "median"
mean_median <- ifelse (mean_median=="mean","mean","median")
agetrans <- AgeTrans(Bdata)
ages <- agetrans$ages
ageentry <- agetrans$ageentry
agecens <- agetrans$agecens
ncase <- c(sort(table(Bdata$path),decreasing=TRUE))  # number of states in this path
case <- names(ncase)  # same as: case <- unique(Bdata$path)
case_ind <- match(Bdata$path,case)
ns_case <- nchar(case)

# Mean/median age at transition and at censoring (M_age_int), by path type
mages <- aggregate(ages,by=list(case_ind),mean_median,na.rm=TRUE)
mages <- round(mages,2)

mageentry <- round(aggregate(ageentry,by=list(case_ind),function(x){mean(x,na.rm=TRUE)}),2)
magecens <- round(aggregate(agecens,by=list(case_ind),mean_median,na.rm=TRUE),2)
sequences <-data.frame(cbind(case,ns_case,ncase,mageentry[,-1],magecens[,-1],mages[,-1]))
colnames(sequences) <- c("case","ns","ncase","M_age_entry","M_age_exit",paste("tr",1:(ncol(Bdata)-locpath(Bdata)),sep="")) # -1 for cohort added dec 09
# sort sequences: sorted dataframe = sequences2 and sequences_sorted
sequences2 <- sequences[order(sequences$ncase,decreasing=TRUE),]
sequences_sorted <- array(" ",c(nrow(sequences2),(ncol(Bdata)-locpath(Bdata))))     # number of transitions !! 
# Add to each mean/median age at transition, the destination state (state occupied after transition) 
for (i in 1:nrow(sequences))
{   for (j in 1:sequences2$ns[i]-1) {sequences_sorted[i,j] <- paste(sequences2[i,j+5],">",substr(sequences2$case[i],j+1,j+1),sep="")}
}
#  create dataframe 
sequences_sorted <- data.frame(sequences_sorted)
# Calculate share and cumulative share of pathway types
percentage <- 100*sequences2$ncase/sum(sequences2$ncase)
cumpercentage <-cumsum(percentage)
# Arrange the columns in dataframe
sequences_sorted <- cbind(sequences2[,3],round(percentage,2),round(cumpercentage,2),
    sequences2[,4],sequences2[,5],sequences2[,2],sequences2[,1],sequences_sorted)
colnames(sequences_sorted) <- c("ncase","%","cum%","M_age_entry","M_age_exit","ns","case",paste("tr",1:(max(sequences$ns)-1),sep=""))
rownames (sequences_sorted) <- c(1:nrow(sequences_sorted))
return(list( mean_median=mean_median,
             sequences=sequences_sorted))
}
