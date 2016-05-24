####  Analysis of AP_MS data: SAINT framed by Pre- and Postprocessing steps  - whole pipeline   ####
# Aim: Identify interaction partners of a bait protein and separation from the contaminants

#Steps:   - Normalization
#         - Filtering
#         - original SAINT run on (normalized)(filtered) data
#         - Permutation approach on SAINT: shuffle of controls and bait samples
#         - Westfall&Young algorithm to receive FWER controlled p-values for each protein interaction-candidate

# Assumptions:
#     - Input of only one bait protein; controls required; - replicates for each group required
#     - bait protein name(V2 in baittable) same for all bait replicates, same for all controls
#     - minimum number of replicates for bait and control experiment is 3 

#####  Input:
##    - Baittable as required for SAINT: classification of bait(T) and ctrls(C) including their names
##    - Interactiontable as required for SAINT
#           * important: a protein which was not detected in a sample receives a zero count
##    - Proteintable = Preyfile as required for SAINT
      ### these 3 files should be stored in the working directory !! ##

#     - choose a normalization method; default="none" (no normalization)
#     - Filtering TRUE/FALSE, set all parameters for filtering (method, cutoff, limit)

#Inputfiles example:
# file_baittable <- "baittable.txt"               # read Baittable  
# file_inttable <- "LBcl_IntSaint.txt"            # read Interactiontable


#Output: - in case of normalization: normed interaction-table (txt file)
#        - in case of filtering: filtered (normed) interaction-table (txt file)
#        - Saint output: "unique_interaction" file on the original data (normalized: "_orig" ;filtered: "_orgF")
#        - permutation matrix perm.avgp, perm.maxp (Rdata) = scores for each protein from all permutation runs
#        - overall result "WY_Result.csv": Saint interaction output + WY adjusted p-values for each candidate



saint_permF  <- function(file_baittable, file_inttable, prottable, norm = c("none", "sumtotal", "upperquartile", "DESeq", "TMM", "quantile"),Filter=TRUE, filter.method= c("IQR", "overallVar", "noVar"), var.cutoff=NA, limit=0, intern.norm=FALSE, saint.options="2000 10000 0 1 0") {
#require(multtest)
#require(gtools)

norm <- match.arg(norm)
saint.options <- match.arg(saint.options)

baittable <- read.table(file_baittable)            # data input
inttable<- read.table(file_inttable)

baittable$V3 <- as.character(baittable$V3)
baittable$V2 <- as.character(baittable$V2)
baittable$V1 <- as.character(baittable$V1)  
inttable$V1 <- as.character(inttable$V1)
inttable$V2 <- as.character(inttable$V2)


if(all(baittable$V3%in%c("C","T"))==FALSE)  {
stop("Define bait and control samples only by C and T in the Baittable \n")    }

if (setequal(names(table(baittable$V2)), names(table(inttable$V2)))==FALSE | setequal(names(table(baittable$V1)), names(table(inttable$V1)))==FALSE ) {
stop("Names of bait and control samples need to be consistent in the Bait- and Interactiontable.  \n")    }


baittab <- baittable[order(baittable$V3),]       # baittab sorted: controls first followed by bait samples
inttable.raw <- inttable



                               #####    NORMALIZATION    #####
if (norm!="none")  { 
inttab.mat <- int2mat(inttable) 
norm.out <- norm.inttable (inttab.mat, baittab, norm)        
inttable <- mat2int(norm.out[[1]], baittab)
baitname <- as.character(unique(baittab$V2[grep("T", baittab$V3)]))
normSaint_filename <- paste(norm, paste(baitname, "IntSaint.txt",sep="_"),sep="_")
write.table(inttable, file=normSaint_filename, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
}


                              #######     FILTERING     #######
if(Filter==TRUE) {  
filter.method <- match.arg(filter.method)              
intmat <- int2mat(inttable)                                   
intmat <- varFilter(intmat, baittab, func=filter.method, var.cutoff, limit)        # FILTER FUNCTION
inttable <- mat2int(intmat, baittab)    # new filtered inttable with selected proteins
                                        
if(intern.norm==TRUE) {                 # repeated normalization on filtered interaction table
del.prots <- setdiff(unique(inttable.raw$V3), unique(inttable$V3)) 
del.pos <- unlist(lapply(del.prots, function(x){which(inttable.raw$V3==x)} ))
intmat.f <- int2mat(inttable.raw[-del.pos,] )
norm2.out <- norm.inttable (intmat.f, baittab, norm)    
inttable <- mat2int(norm2.out[[1]], baittab)
}

write.table(inttable, file="Inttable_filtered.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

# 1.SAINT run on original data on a filtered protein-set:
saint.command <- paste("saint-spc-ctrl", "Inttable_filtered.txt", prottable, file_baittable, saint.options,sep=" ") 
system("rm -r LOG")                       
system("rm -r MAPPING")          # deletion of former folders
system("rm -r MCMC")
system("rm -r RESULT")
system(saint.command)            # linux run
#SAINT output "unique_interactions":
filename.org <- "./RESULT/unique_interactions"
orig.filter <- new.saintoutput(filename.org)
write.table(orig.filter, file="Unique_Interactions_orgF.txt", quote=FALSE, sep="\t", row.names=FALSE)

proteins <- as.character(orig.filter$Prey)
}
##  NO Filtering:
else {      
    if (norm=="none")  {                 # 1.SAINT run on original data (unnormalized + no filter) 
    saint.command2 <- paste("saint-spc-ctrl", file_inttable, prottable, file_baittable, saint.options,sep=" ") 
    system("rm -r LOG")                       
    system("rm -r MAPPING")         
    system("rm -r MCMC")
    system("rm -r RESULT")
    system(saint.command2)   }        
    if (norm!="none")  {                # 1.SAINT run on original data (normalized+ no filter) 
    saint.command.norm <- paste("saint-spc-ctrl", normSaint_filename, prottable, file_baittable, saint.options,sep=" ") 
    system("rm -r LOG")                       
    system("rm -r MAPPING")          
    system("rm -r MCMC")
    system("rm -r RESULT")
    system(saint.command.norm)    }
  
  #SAINT output "unique_interactions":
  filename.org <- "./RESULT/unique_interactions"
  orig.out <- new.saintoutput(filename.org)
  write.table(orig.out, file="Unique_Interactions_Orig.txt", quote=FALSE, sep="\t", row.names=FALSE)
  proteins <- as.character(orig.out$Prey)
}


                            #########    Permutation     ###########

ctrl <- as.character(baittab$V1[grep("C", baittab$V3)])               # control replicates
ctrl.sup <- unique(as.character(baittab$V2[grep("C", baittab$V3)]))   #control sup-name
baits <- as.character(baittab$V1[grep("T", baittab$V3)])              # bait replicates
bait.sup <- unique(as.character(baittab$V2[grep("T", baittab$V3)]))   #superior bait name

# permutation setup:
shuf.vec <- baittab$V3
classlabel <- ifelse(shuf.vec=="C",0,1)
perms <- mt.sample.label(classlabel,test="t",fixed.seed.sampling="n",B=0)
perms2 <- apply(perms, 2, function(x){ifelse(x==0,"C","T")} )
perms.baitname <- apply( perms, 2, function(x) {ifelse(x==0,baittab$V2[1],baittab$V2[length(baittab$V2)]) } )

# permutation table to store the resulting scores from each permutation:
perm.avgp <- as.data.frame(matrix(data=NA, nrow=length(proteins), ncol=(dim(perms)-1) ) )
rownames(perm.avgp) <- proteins
perm.maxp <- as.data.frame(matrix(data=NA, nrow=length(proteins), ncol=(dim(perms)-1) ) )
rownames(perm.maxp) <- proteins

#### Each permutation: create new baittable + interactiontable -> new Saint run -> store resulting scores
#1.row = original data setup (skipped)

for ( j in 2:nrow(perms2)) {                          
cat(j, ". Run started \n", sep="")
system("rm -r LOG")
system("rm -r MAPPING")                   # delete folders for next permutation run 
system("rm -r MCMC")
system("rm -r RESULT")

baittab.perm <- baittab
baittab.perm$V3 <- perms2[j,]
baittab.perm$V2 <- perms.baitname[j,]     # permuted baittable

new.ctrls <- baittab.perm$V1[grep("C",baittab.perm$V3)]
new.ctrl <- new.ctrls [which(new.ctrls %in% ctrl ==FALSE) ]   # bait turned into control 
new.baits <- baittab.perm$V1[grep("T",baittab.perm$V3)]
new.bait <- new.baits [which(new.baits %in% baits ==FALSE) ]  # control turned into bait 

perm.inttable <- inttable                                     # permuted inttable
newctrl.pos <- unlist(sapply(new.ctrl, function(x){which(perm.inttable$V1==x)}))  
perm.inttable$V2[newctrl.pos] <- ctrl.sup                                         
newbait.pos <- unlist(sapply(new.bait, function(x){which(perm.inttable$V1==x)}))
perm.inttable$V2[newbait.pos] <- bait.sup


print(table(perm.inttable$V1, perm.inttable$V2))              # check permutation

write.table(baittab.perm, file="Baittable_perm.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
write.table(perm.inttable, file="Inttable_perm.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

### SAINT run on permutated data
saint.commandp1 <- paste("saint-spc-ctrl", "Inttable_perm.txt", prottable, "Baittable_perm.txt",saint.options,sep=" ") 
system(saint.commandp1)

# SAINT output "unique_interactions"
filename <- "./RESULT/unique_interactions"
saintout.perm <- try(new.saintoutput(filename))
while (class(saintout.perm)=="try-error")         # case of file generation failure
      { system("rm -r LOG")                       # do it again
        system("rm -r MAPPING")          
        system("rm -r MCMC")
        system("rm -r RESULT")
        saint.commandp1 <- paste("saint-spc-ctrl", "Inttable_perm.txt", prottable, "Baittable_perm.txt", saint.options,sep=" ") 
	 system(saint.commandp1) 
        saintout.perm <- try(new.saintoutput(filename))
      }
saintout.perm$AvgP <- as.numeric(as.vector(saintout.perm$AvgP))
saintout.perm$MaxP <- as.numeric(as.vector(saintout.perm$MaxP))
saintout.perm$Prey <- as.character(saintout.perm$Prey)

# extract scores from output "unique-interactions"
protpos.perm <- match(proteins, saintout.perm$Prey) 
perm.avgp[,(j-1)] <- saintout.perm$AvgP[protpos.perm]
perm.maxp[,(j-1)] <- saintout.perm$MaxP[protpos.perm]

}          
# permutation end

save(perm.avgp, file="perm_avgP.RData")
save(perm.maxp, file="perm_maxP.RData")

system("rm -r LOG")                       
system("rm -r MAPPING")          
system("rm -r MCMC")
system("rm -r RESULT")
                          #######     Generation p-values    #######

if(Filter==TRUE) {orig <- orig.filter} else  {orig <- orig.out}

orig$AvgP <- as.numeric(as.vector(orig$AvgP))
set.seed(12345)
pp <- permute(c(1:length(proteins)))
wy <- WY.permalg(orig$AvgP[pp], perm.avgp[pp,])
tor <- rev(order(abs(orig$AvgP)))
write.csv2(cbind(orig[tor,] ,"WY_counter"=wy[match(orig[tor,]$Prey,rownames(wy)),1],"WY_P.adjust"=wy[match(orig[tor,]$Prey,rownames(wy)),2]) ,file="WY_Result.csv")

}

# end 


