setClass("GENOME", representation(

BIG.BIAL     =   "list", # length 1 for ff objects
SLIDE.POS    =   "list",
big.data     =   "logical",
gff.info     =   "logical",
snp.data     =   "logical",
basepath     =   "character",
project      =   "character",
populations  =   "list",        # populations
poppairs     =   "vector",      # population pairs
outgroup     =   "vector",      # outgroup
#region.names    =   "character",   # region.names
region.names =   "character",
feature.names=   "character",
genelength   =   "numeric",     # number of genes that where calculated
n.sites      =   "numeric",
n.sites2     =   "numeric", # important for SNP data
n.biallelic.sites = "numeric",
n.gaps       =   "numeric",
n.unknowns   =   "numeric",
n.valid.sites       =   "numeric",
n.polyallelic.sites =   "numeric",
trans.transv.ratio  =   "numeric",
keep.start.pos      =   "numeric", 

## GFF Infos 

Coding.region ="numeric",
UTR.region    ="numeric",
Intron.region ="numeric",
Exon.region   ="numeric",
Gene.region   ="numeric",

######## Pops for each Test Module

Pop_Neutrality  =  "list",
Pop_FSTN        =  "list",
Pop_FSTH        =  "list",
Pop_Linkage     =  "list",
Pop_Slide       =  "list", 
Pop_MK          =  "list",
Pop_Detail      =  "list",
Pop_Recomb      =  "list",
Pop_Sweeps      =  "list",

##################################
FSTNLISTE              =   "list",
nucleotide.F_ST        =   "matrix",      # FSTALL for Nucleotides # FSTN
nucleotide.F_ST2       =   "matrix",      # from calc_hawhfsth
nuc.diversity.between  =   "matrix",      # Nucleotide Diversity between populations # PIA
nuc.diversity.within   =   "matrix",      # Nucleotide Diversity within  populations # PIW
nuc.F_ST.pairwise      =   "matrix",      # FST FOR POP-PAIRS (Nucleotides) # FSTNPAIR
nuc.F_ST.vs.all        =   "matrix",      # nuc
n.haplotypes           =   "matrix",      # number of haplotypes in each locus # nh
hap.diversity.within   =   "matrix",      # haplotype diversity within populations # hapw
hap.diversity.between  =   "matrix",      # hapa
Pi                     =   "matrix",      # Nucleotide Diversity
PIA_nei                =   "matrix",
haplotype.counts       =   "matrix",      # copies of each haplotype #sfreqh
haplotype.F_ST         =   "matrix",      # FSTALL for Haplotypes # FSTH
hap.F_ST.pairwise      =   "matrix",
Nei.G_ST.pairwise      =   "matrix",
hap.F_ST.vs.all        =   "matrix", 
Nei.G_ST               =   "matrix",      # GST from Nei # GST
Hudson.G_ST            =   "matrix",      # GST from Hudson (Nei) # GSTH
Hudson.H_ST            =   "matrix",      # HST Hudson # HST
Hudson.K_ST            =   "matrix",      # KST Hudson # KST

Hudson.Snn             =   "matrix",      # F_ST.stats.2 module
Phi_ST                 =   "matrix",      # F_ST.stats.2 module

hap.pair.F_ST          =   "matrix",       # FST FOR POP-PAIRS (Haplotypes) # FSTHPAIR 
#Ps          =   "matrix",      # Synonymous Sites (within the population)
#Pn          =   "matrix",      # Nonsynonymous Sites (within the population)
#Ds          =   "matrix",      # Syn Position between population pairs      
#Dn          =   "matrix",      # NonSyn Positions of poppairs 
MKT                    =   "matrix",      # Mcdonald Kreitman Test
Tajima.D               =   "matrix",      # TAJIMA D VALUES #  TajD
SLIDE                  =   "matrix",      # Sliding Window Tajima
Fay.Wu.H               =   "matrix",      # Fay Wu normalized
Zeng.E                 =   "matrix",      # E Zeng
theta_Tajima           =   "matrix",      # thetaT (Tajima) # thetaT
theta_Watterson        =   "matrix",      # Watterson # thetaS
theta_Fu.Li            =   "matrix", # thetaFL
theta_Achaz.Watterson  =   "matrix", # thetaSA
theta_Achaz.Tajima     =   "matrix", # thetaTA
theta_Fay.Wu           =   "matrix", # thetaFW
theta_Zeng             =   "matrix", # thetaL

Fu.Li.F                =   "matrix",      # Fu Li F* # FuLi_F
Fu.Li.D     	       =   "matrix",      # Fu Li D* # FuLi_D
Yach                   =   "matrix",      # Achaz
n.segregating.sites    =   "matrix",      # Segregating Sites # S
Rozas.R_2              =   "matrix",      # R2 Values
Fu.F_S                 =   "matrix",      # Fu's FS # FS
Strobeck.S             =   "matrix",

Kelly.Z_nS             =   "matrix", # Zns
Rozas.ZZ               =   "matrix",
Rozas.ZA               =   "matrix",
Wall.B                 =   "matrix",
Wall.Q                 =   "matrix",
mult.Linkage           =   "matrix",

RM                     =   "matrix",
CL                     =   "matrix",
CLmax                  =   "matrix",
CLR                    =   "matrix",

MDSD="matrix",
MDG1="matrix",
MDG2="matrix",

D		       =   "matrix",              # introgression slots
f		       =   "matrix",              
jack.knife             =   "logical",
missing.freqs          =   "matrix",

genes                  =   "list",                # a list of statistics objects
region.data            =   "region.data",         # list of class GEN
region.stats           =   "region.stats"         # list of class DATA


))

### Get Summary DATA #######################

 setGeneric("get.sum.data", function(object) standardGeneric("get.sum.data"))
 setMethod("get.sum.data", "GENOME",
 function(object){

res           <- cbind(object@n.sites,object@n.biallelic.sites,object@n.gaps,object@n.unknowns,object@n.valid.sites,object@n.polyallelic.sites,object@trans.transv.ratio)
rownames(res) <- object@region.names
colnames(res) <- c("n.sites","n.biallelic.sites","n.gaps","n.unknowns","n.valid.sites","n.polyallelic.sites","trans.transv.ratio")

 return(res)
 })


#############################################
### Show ####################################

setMethod("show", "GENOME",
 function(object){
 cat("-----\n")
 cat("Modules:\n")
 cat("-----\n")
 out <- data.frame(Calculation=c("readData","neutrality.stats","linkage.stats","recomb.stats","F_ST.stats","diversity.stats","sweeps.stats","MKT","detail.stats","MS","--------------","set.populations","sliding.window.transform","splitting.data","show.slots","get.status"),         
 
                Description=c("Reading data","Neutrality tests","Linkage disequilibrium","Recombination","Fixation index","Diversities","Selective sweeps","McDonald-Kreitman test","Mixed statistics","Coalescent simulation",
"-----------","Defines the populations","Sliding window","Splits the data","?provided slots?","Status of calculations"),
                
                Get.the.Result = c("get.sum.data","get.neutrality","get.linkage","get.recomb","get.F_ST,get.diversity","get.diversity","get.sweeps","get.MKT","get.detail","@","-------------","","","","","")

                # ,

                #Calculated = c(TRUE,object@Pop_Neutrality$calculated,object@Pop_Linkage$calculated,object@Pop_Recomb #$calculated,object@Pop_FSTH$calculated,object@Pop_MK$calculated,
#object@Pop_Detail$calculated,NA,NA,NA,object@Pop_Slide$calculated,NA,NA)
                )
  print(out)

# cat("\n")
#  cat("--------- Help ----------------")
#  cat("\n")
#  cat("usage(GENOME.class)\n")
#  cat("show.slots(GENOME.class)\n")
#  cat("get.status(GENOME.class)\n")
#  cat("-------------------------------\n") 
#  cat("Note ! GENOME.class ist an object of class::GENOME returned from readData. \n")

})


#### showSLots ######
setGeneric("show.slots", function(object) standardGeneric("show.slots"))
setMethod("show.slots", "GENOME",
 function(object){
 cat("-----\n")
 cat("SLOTS:\n")
 cat("-----\n")
 out <- data.frame(Slots=c("project","populations","outgroup","region.names","n.sites","Pi","nucleotide.F_ST","nuc.F_ST.pairwise",
"haplotype.F_ST","Nei.G_ST","Hudson.G_ST","Hudson.H_ST","hap.F_ST.pairwise",
"Tajima.D","Fu.Li.F","Fu.Li.D","Fu.F_S","Strobeck","n.segragting.sites","Rozas.R_2",
"Kelly.Z_nS","Rozas.ZA","Rozas.ZZ","Wall.B","Wall.Q","RM","CL","CLmax",
"MDSD","MDG1","MDG2","Coding.region","Gene.region","UTR.region","Intron.region","Exon.region","region.stats","region.data"),         
 
                Description=c("project name","list of populations","outgroup vector","region names","length of each region","Nucleotide Diversity (Nei)","F_ST (Nucleotides)",
                "Pairwise F_ST (Nucleotides)","F_ST (Haplotypes)","G_ST from Nei 1973","G_ST from Hudson 1992","H_ST from Hudson 1992","Pairwise F_ST (Haplotypes)","Tajima's D","Fu & Li's F*","Fu & Li's D*","Fu & Li's F_S","Strobeck","Segregating Sites","Ramos & Rozas R_2","Linkage Disequilibrium","Rozas ZA statistic","Rozas ZZ statistic","Walls B","Walls Q","Hudson & Kaplan RM","Composite Likelihood (Nielson)","Max Composite Likelihood (Nielson)","mismatch distribution","mismatch distribution","mismatch distribution","length of Coding region","length of Gene region","length of UTR region","length of Intron region","length of Exon region",
                "  detail stats (class::region.stats)"," summary data (class::region.data)"),
                
                Module = c("Data","Data","Data","Data","Data","FST","FST","FST","FST","FST","FST","FST","FST","FST","Neutrality","Neutrality","Neutrality",
"Neutrality","Neutrality","Neutrality","Linkage","Linkage","Linkage","Linkage","Linkage","Recomb","Sweeps","Sweeps","Detail","Detail","Detail",
"Gff/Gtf","Gff/Gtf","Gff/Gtf","Gff/Gtf","Gff/Gtf","!!!","!!!"))
  
  print(out)
 })
########################################################################
 
setGeneric("usage", function(object) standardGeneric("usage"))
setMethod("usage", "GENOME",
 function(object){

 out <- data.frame(Functions=c("readData","neutrality.stats","linkage.stats","F_ST.stats","MKT","detail.stats","MS","------------","sliding.window.transform","set.populations"),         
 
                   Usage=c("GENOME.class <- readData(path)","GENOME.class <- neutrality.stats(GENOME.class)","GENOME.class <- linkage.stats(GENOME.class)",
"GENOME.class <- F_ST.stats(GENOME.class)","GENOME.class <- MKT(GENOME.class)","GENOME.class <- detail.stats(GENOME.class)","MS.class <- MS(GENOME.class)","---------------------------------------","slide.GENOME.class <- sliding.window.transform(GENOME.class)","GENOME.class <- set.populations(GENOME.class,list(pop1,pop2,...,popN))"))
                


  print(out)
 })


#### overview ##########################################################
#setGeneric("overview", function(object) standardGeneric("overview"))
#setMethod("overview", "GENOME",

# function(object){
 
 # Summary Data
# datasumset             <- matrix(,3,1)
# colnames(datasumset)   <- "Data.Summary"
 #rownames(datasumset)   <- c("Name of the Genome","Number of genes","Total number of sites")
 #datasumset[1,]         <- object@genome
# datasumset[2,]         <- object@genelength
# datasumset[3,]         <- sum(object@n.sites)

# Summary Statistics
 #npops <- length(object@Pop_Neutrality$Populations)
# out   <- vector("list",npops)
 
#for(xx in 1:npops){
# out[[xx]] <- summary(get.neutrality.stats(object)[[xx]]) # beim Auruf von get.neutrality.stats gibts bei R check Probleme !
#}

 #fst  <- summary(getFST(object)) 

#return(list(Summary.Data=datasumset,Summary.Neutrality=out,Summary.FST=fst))
 
#})
 

# Get Status -----------------------------------------------------------

setGeneric("get.status", function(object) standardGeneric("get.status"))
 setMethod("get.status", "GENOME",
 function(object){
 
 cat("########################\n")
 cat("Calculated Statistics: \n")
 cat("########################\n")
 cat("\n")
 cat("########################\n")
 cat("Neutrality::neutrality.stats \n")
 cat("------------------ \n")
 
 cat("Populations:\n")
 print(object@Pop_Neutrality$Populations)
 cat("Sites:",object@Pop_Neutrality$sites,"\n")
 #cat("Calculated:",object@Pop_Neutrality$calculated,"\n")
 
 cat("########################\n")
 #cat("\n")
 cat("F_ST::F_ST stats \n")
 cat("------------------ \n")
 
 #cat("Populations:\n")
 #print(object@Pop_FSTN$Populations)
 #cat("Sites:",object@Pop_FSTN$sites,"\n")
 #cat("Calculated:",object@Pop_FSTN$calculated,"\n")
 
 # cat("########################\n")
 #cat("\n")
 #cat("Haplotype Module: (popFSTH()) \n")
 #cat("------------------ \n")
 
 cat("Populations:\n")
 print(object@Pop_FSTH$Populations)
 cat("Sites:",object@Pop_FSTH$sites,"\n")
 
 #cat("Calculated:",object@Pop_FSTH$calculated,"\n")
 
 cat("########################\n") 
 #cat("\n")
 cat("McDonald & Kreitman::MKT \n")
 cat("------------------ \n")
 
 cat("Populations:\n")
 print(object@Pop_MK$Populations)
 cat("\n")
 #cat("Calculated:",object@Pop_MK$calculated,"\n")
 
 #cat("########################\n") 
 #cat("\n")
 #cat("Sliding Window::sliding.window \n")
 #cat("------------------ \n")
 
 #cat("Populations:\n")
 #print(object@Pop_Slide$Populations)
 #cat("Calculated:",object@Pop_Slide$calculated,"\n")
 
 cat("########################\n") 
 #cat("\n")
 cat("Linkage Disequilibrium::linkage.stats \n")
 cat("------------------ \n")
 
 cat("Populations:\n")
 print(object@Pop_Linkage$Populations)
 cat("Sites:",object@Pop_Linkage$sites,"\n")
 #cat("Calculated:",object@Pop_Linkage$calculated,"\n")
 #cat("########################\n")



 cat("########################\n") 
 #cat("\n")
 cat("Mismatch Distribution::detail.stats \n")
 cat("------------------ \n")
 
 cat("Populations:\n")
 print(object@Pop_Detail$Populations)
 #cat("Calculated:",object@Pop_Detail$calculated,"\n")
 cat("########################\n")


 }) 
 

## Set Population
#--------------------------------------------
setGeneric("set.populations", function(object,new.populations=FALSE, diploid=FALSE,triploid=FALSE,tetraploid=FALSE) standardGeneric("set.populations"))
 setMethod("set.populations", "GENOME",
 function(object,new.populations,diploid,triploid,tetraploid){

npops              <- length(new.populations)
change             <- object@region.data
populations        <- vector("list",npops)
populations2       <- vector("list",npops)

# if diploid add individuals.2
if(diploid){
 for (yy in 1:npops) {
  dottwo                <- paste(new.populations[[yy]],".2", sep="")
  new.populations[[yy]] <- c(new.populations[[yy]],dottwo) 
 } 
}
if(triploid){
 for (yy in 1:npops) {
  dottwo                <- paste(new.populations[[yy]],".2", sep="")
  dotthree              <- paste(new.populations[[yy]],".3", sep="")
  new.populations[[yy]] <- c(new.populations[[yy]],dottwo,dotthree) 
 } 
}
if(tetraploid){
 for (yy in 1:npops) {
  dottwo                <- paste(new.populations[[yy]],".2", sep="")
  dotthree              <- paste(new.populations[[yy]],".3", sep="")
  dottfour              <- paste(new.populations[[yy]],".4", sep="")
  new.populations[[yy]] <- c(new.populations[[yy]],dottwo,dotthree,dottfour) 
 } 
}

# End of diploid

object@populations <- new.populations

# Init

XX_popmissing   <- vector("list",length(object@region.names))
XX_populations  <- vector("list",length(object@region.names))
XX_populations2 <- vector("list",length(object@region.names))

############################################################
progr <- progressBar()
############################################################


for(xx in 1:length(object@region.names)){


     bial <- object@region.data@biallelic.matrix[[xx]] # popGetBial(object,xx) # if it does not fit into RAM
     
     if(length(bial)==0){next}   
  
      for(yy in 1:npops){

           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]]  <- match(new.populations[[yy]],rownames(bial))
              naids              <- which(!is.na(populations[[yy]]))
              populations[[yy]]  <- populations[[yy]][naids]
              populations2[[yy]] <- rownames(bial)[populations[[yy]]]

           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }      
       }
           
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   

       # if(length(populations)==0){next} # Keine Population vorhanden

     XX_populations[[xx]]  <- populations 
     XX_populations2[[xx]] <- populations2 

     XX_popmissing[[xx]]   <- popmissing
      
 # PROGRESS #######################################################
    progr <- progressBar(xx,length(object@region.names), progr)
 ###################################################################

}

change@populations  <- XX_populations
change@populations2 <- XX_populations2
change@popmissing   <- XX_popmissing
object@region.data  <- change
return(object)

})



# ---------------------------------------------------------------------- 
# Get MS: TAJ S R2 FuLi_F FuLi_D
# ----------------------------------------------------------------------
 
 setGeneric("getMS", function(object,gene=FALSE,neutrality=TRUE,linkage=FALSE,F_ST=FALSE) standardGeneric("getMS"))
 setMethod("getMS", "GENOME",
 function(object,gene=FALSE,neutrality,linkage,F_ST){


if(neutrality){
   
  if(!object@Pop_Neutrality$calculated){stop("Statistics have to be calculated first !")}

  if(gene[1]==FALSE){
  
  res  <- matrix(,object@genelength,11)
  pops <- vector("list",length(object@Pop_Neutrality$Populations))
  for(xx in 1:length(object@Pop_Neutrality$Populations)){
     res[,1]       <- object@Tajima.D[,xx]
     res[,2]       <- object@n.segregating.sites[,xx]
     res[,3]       <- object@Rozas.R_2[,xx]
     res[,4]       <- object@Fu.Li.F[,xx]
     res[,5]       <- object@Fu.Li.D[,xx] 
     res[,6]       <- object@Fu.F_S[,xx]
     res[,7]       <- object@Strobeck.S[,xx]
     res[,8]       <- object@Fay.Wu.H[,xx]
     res[,9]       <- object@Zeng.E[,xx]
     res[,10]      <- object@theta_Tajima[,xx]
     res[,11]      <- object@theta_Watterson[,xx]
     colnames(res) <- c("Tajima.D","n.segregating.sites","Rozas.R_2","Fu.Li.F","Fu.Li.D","Fu.F_S","Strobeck.S","Fay.Wu.H","Zeng.E","theta_Tajima","theta_Watterson")
     rownames(res) <- object@region.names
     pops[[xx]]    <- res 
  }
 pops           <- as.matrix(pops)
 rownames(pops) <- paste("pop",1:length(object@Pop_Neutrality$Populations))
# colnames(pops) <- "Genome"
 return(pops)
 
 }else{               # gene is is defined
 
  res  <- matrix(,1,11)
  pops <- vector("list",length(object@Pop_Neutrality$Populations))
  for(xx in 1:length(object@Pop_Neutrality$Populations)){
     res[,1]       <- object@Tajima.D[gene,xx]
     res[,2]       <- object@n.segregating.sites[gene,xx]
     res[,3]       <- object@Rozas.R_2[gene,xx]
     res[,4]       <- object@Fu.Li.F[gene,xx]
     res[,5]       <- object@Fu.Li.D[gene,xx] 
     res[,6]       <- object@Fu.F_S[gene,xx]
     res[,7]       <- object@Strobeck.S[gene,xx]
     res[,8]       <- object@Fay.Wu.H[gene,xx]
     res[,9]       <- object@Zeng.E[gene,xx]
     res[,10]      <- object@theta_Tajima[gene,xx]
     res[,11]      <- object@theta_Watterson[gene,xx]
     colnames(res) <- c("Tajima.D","n.segregating.sites","Rozas.R_2","Fu.Li.F","Fu.Li.D","Fu.F_S","Strobeck.S","Fay.Wu.H","Zeng.E","theta_Tajima","theta_Watterson")
     rownames(res) <- object@region.names[gene]
     pops[[xx]]    <- res 
  }
 pops           <- as.matrix(pops)
 rownames(pops) <- paste("pop",1:length(object@Pop_Neutrality$Populations))
# colnames(pops) <- "Genome"
 return(pops)
 
 }
 
}# end if neutrality


if(linkage){

if(!object@Pop_Linkage$calculated){stop("Statistics have to be calculated first !")}

  if(gene[1]==FALSE){
  
  res  <- matrix(,object@genelength,5)
  pops <- vector("list",length(object@Pop_Linkage$Populations))
  for(xx in 1:length(object@Pop_Linkage$Populations)){
     res[,1]       <- object@Wall.B[,xx]
     res[,2]       <- object@Wall.Q[,xx]
     res[,3]       <- object@Rozas.ZA[,xx]
     res[,4]       <- object@Rozas.ZZ[,xx]
     res[,5]       <- object@Kelly.Z_nS[,xx] 
     colnames(res) <- c("Wall.B","Wall.Q","Rozas.ZA","Rozas.ZZ","Kelly.Z_nS")
     rownames(res) <- object@region.names
     pops[[xx]]    <- res 
  }

 pops           <- as.matrix(pops)
 rownames(pops) <- paste("pop",1:length(object@Pop_Linkage$Populations))
# colnames(pops) <- "Genome"
 return(pops)
 
 }else{               # gene is is defined
 
  res  <- matrix(,1,5)
  pops <- vector("list",length(object@Pop_Linkage$Populations))
  for(xx in 1:length(object@Pop_Linkage$Populations)){
     res[,1]       <- object@Wall.B[gene,xx]
     res[,2]       <- object@Wall.Q[gene,xx]
     res[,3]       <- object@Rozas.ZA[gene,xx]
     res[,4]       <- object@Rozas.ZZ[gene,xx]
     res[,5]       <- object@Kelly.Z_nS[gene,xx] 
     colnames(res) <- c("Wall.B","Wall.Q","Rozas.ZA","Rozas.ZZ","Kelly.Z_nS")
     rownames(res) <- object@region.names[gene]
     pops[[xx]]    <- res 
  }

 pops           <- as.matrix(pops)
 rownames(pops) <- paste("pop",1:length(object@Pop_Linkage$Populations))
# colnames(pops) <- "Genome"
 return(pops)
 
 }
}# end if linkage

if(F_ST){

  if(length(object@Hudson.Snn)==0){object@Hudson.Snn <- matrix(,object@genelength,1)} # damit Snn nicht unbedingt berechnet werden muss


if(!object@Pop_FSTH$calculated){stop("Statistics have to be calculated first !")}

  if(gene[1]==FALSE){
  
  res  <- matrix(,object@genelength,6)
  pops <- vector("list",length(object@Pop_FSTH$Populations))
  for(xx in 1:length(object@Pop_FSTH$Populations)){
     res[,1]       <- object@hap.diversity.within[,xx]
     res[,2]       <- object@Pi[,xx] #*object@n.valid.sites
     res[,3]       <- object@haplotype.F_ST
     res[,4]       <- object@nucleotide.F_ST
     res[,5]       <- object@Nei.G_ST 
     res[,6]       <- object@Hudson.Snn
     colnames(res) <- c("hap.diversity.within","Pi","haplotype.F_ST","nucleotide.F_ST","Nei.G_ST","Hudson.Snn")
     rownames(res) <- object@region.names
     pops[[xx]]    <- res 
  }

 pops           <- as.matrix(pops)
 rownames(pops) <- paste("pop",1:length(object@Pop_FSTH$Populations))
# colnames(pops) <- "Genome"
 return(pops)
 
 }else{               # gene is is defined
 
  res  <- matrix(,1,6)
  pops <- vector("list",length(object@Pop_FSTH$Populations))
  for(xx in 1:length(object@Pop_FSTH$Populations)){
     res[,1]       <- object@hap.diversity.within[gene,xx]
     res[,2]       <- object@Pi[gene,xx]*object@n.valid.sites[gene]
     res[,3]       <- object@haplotype.F_ST[gene,1]
     res[,4]       <- object@nucleotide.F_ST[gene,1]
     res[,5]       <- object@Nei.G_ST[gene,1] 
     res[,6]       <- object@Hudson.Snn[gene,1]
     colnames(res) <- c("hap.diversity.within","Pi","haplotype.F_ST","nucleotide.F_ST","Nei.G_ST","Hudson.Snn")
     rownames(res) <- object@region.names[gene]
     pops[[xx]]    <- res 
  }

 pops           <- as.matrix(pops)
 rownames(pops) <- paste("pop",1:length(object@Pop_FSTH$Populations))
# colnames(pops) <- "Genome"
 return(pops)
 
 }
}# end if linkage

 })
 
 # -----------------------------------------------------------------------
 # Get Neutrality
 # -----------------------------------------------------------------------
 
 setGeneric("get.neutrality", function(object,theta=FALSE,stats=TRUE) standardGeneric("get.neutrality"))
 setMethod("get.neutrality", "GENOME",
 
 function(object,theta,stats){
 
 if(!object@Pop_Neutrality$calculated){stop("Statistics have to be calculated first !")}

 if(stats){
  res   <- matrix(,object@genelength,9)
  pops1 <- vector("list",length(object@Pop_Neutrality$Populations))
  
  for(xx in 1:length(object@Pop_Neutrality$Populations)){
     res[,1]       <- object@Tajima.D[,xx]
     res[,2]       <- object@n.segregating.sites[,xx]
     res[,3]       <- object@Rozas.R_2[,xx]
     res[,4]       <- object@Fu.Li.F[,xx]
     res[,5]       <- object@Fu.Li.D[,xx]
     res[,6]       <- object@Fu.F_S[,xx] 
     res[,7]       <- object@Fay.Wu.H[,xx]
     res[,8]       <- object@Zeng.E[,xx]
     res[,9]       <- object@Strobeck.S[,xx]
     colnames(res) <- c("Tajima.D","n.segregating.sites","Rozas.R_2","Fu.Li.F","Fu.Li.D","Fu.F_S","Fay.Wu.H","Zeng.E","Strobeck.S")
     rownames(res) <- object@region.names
     pops1[[xx]]   <- res 
  }
 
 pops1            <- as.matrix(pops1)
 rownames(pops1)  <- paste("pop",1:length(object@Pop_Neutrality$Populations))
 colnames(pops1)   <- "neurality stats"
 if(!theta){
    return(pops1)
 }
 
 }
 
 if(theta){

  res   <- matrix(,object@genelength,7)
  pops2 <- vector("list",length(object@Pop_Neutrality$Populations))


 for(xx in 1:length(object@Pop_Neutrality$Populations)){
    
     res[,1]       <- object@theta_Tajima[,xx]
     res[,2]       <- object@theta_Watterson[,xx]
     res[,3]       <- object@theta_Fu.Li[,xx]
     res[,4]       <- object@theta_Achaz.Watterson[,xx]
     res[,5]       <- object@theta_Fay.Wu[,xx]
     res[,6]       <- object@theta_Zeng[,xx] 
     res[,7]       <- object@theta_Achaz.Tajima[,xx]
     
     #res[,8]       <- object@thetaS[,xx]
     colnames(res)  <- c("theta_Tajima","theta_Watterson","theta_Fu.Li","theta_Achaz.Watterson","theta_Fay.Wu","theta_Zeng","theta_Achaz.Tajima")
     rownames(res)  <- object@region.names
     pops2[[xx]]    <- res 
  
}
 
 pops2            <- as.matrix(pops2)
 rownames(pops2)  <- paste("pop",1:length(object@Pop_Neutrality$Populations))
 colnames(pops2)  <- "theta values"
 
 if(!stats){
    return(pops2)
 }
  
 }

 if(stats & theta){
 
  npops  <- length(object@Pop_Neutrality$Populations)
  allpop <- vector("list",npops)
 
 for(xx in 1: npops){
    
  allpop[[xx]] <- cbind(pops1[[xx]],pops2[[xx]])
 
 }
   allpop             <- as.matrix(allpop)
   rownames(allpop)   <- paste("pop",1:length(object@Pop_Neutrality$Populations))
   colnames(allpop)   <- "neutrality stats + theta values"
   return(allpop)
 }
 
 })

### GetLinkage

setGeneric("get.linkage", function(object) standardGeneric("get.linkage"))
 setMethod("get.linkage", "GENOME",
 
 function(object){

if(!object@Pop_Linkage$calculated){stop("Statistics have to be calculated first !")}

 res  <- matrix(,object@genelength,5)

  pops <- vector("list",length(object@Pop_Linkage$Populations))
  
  for(xx in 1:length(object@Pop_Linkage$Populations)){
     res[,1]       <- object@Wall.B[,xx]
     res[,2]       <- object@Wall.Q[,xx]
     res[,3]       <- object@Rozas.ZA[,xx]
     res[,4]       <- object@Rozas.ZZ[,xx]
     res[,5]       <- object@Kelly.Z_nS[,xx]

     colnames(res) <- c("Wall.B","Wall.Q","Rozas.ZA","Rozas.ZZ","Kelly.Z_nS")
     rownames(res) <- object@region.names
     pops[[xx]] <- res 
  }
 
 pops <- as.matrix(pops)
 rownames(pops)  <- paste("pop",1:length(object@Pop_Linkage$Populations))
 colnames(pops)  <- "Linkage Disequilibrium"
 
 return(pops)
 
 
 return(res)
 }) 
 
#######################################
### get.diversity

setGeneric("get.diversity", function(object,between=FALSE) standardGeneric("get.diversity"))
 setMethod("get.diversity", "GENOME",
 
 function(object,between){

if(!object@Pop_FSTH$calculated){stop("F_ST have to be calculated first !")}


if(!between){

  res  <- matrix(,object@genelength,5)

  pops <- vector("list",length(object@Pop_FSTH$Populations))
  
  for(xx in 1:length(object@Pop_FSTH$Populations)){

     res[,1]       <- object@nuc.diversity.within[,xx]
     res[,2]       <- object@hap.diversity.within[,xx]
     if(length(object@Pi)!=0){
     res[,3]       <- object@Pi[,xx]
     }
     if(length(object@hap.F_ST.vs.all)!=0){
     res[,4]       <- object@hap.F_ST.vs.all[,xx]
     }	
     if(length(object@nuc.F_ST.vs.all)!=0){
     res[,5]       <- object@nuc.F_ST.vs.all[,xx]
     }

     colnames(res) <- c("nuc.diversity.within","hap.diversity.within","Pi","hap.F_ST.vs.all","nuc.F_ST.vs.all")
     rownames(res) <- object@region.names
     pops[[xx]] <- res 
  }
 
 pops <- as.matrix(pops)
 rownames(pops)  <- paste("pop",1:length(object@Pop_FSTH$Populations))
 colnames(pops)  <- "Diversity"
 return(pops)
 
}else{

 res          <- vector("list",2)
 res[[1]]     <- t(object@hap.diversity.between)
 res[[2]]     <- t(object@nuc.diversity.between)

res           <- as.matrix(res)
rownames(res) <- c("hap.diversity.between","nuc.diversity.between")
colnames(res) <- "Diversity"
 
return(res)

} 
 

}) 
 
####################################### 
# -----------------------------------------------------------------------
# GetBayes
# ----------------------------------------------------------------------
#######################################

setGeneric("getBayes", function(object,snps=FALSE) standardGeneric("getBayes"))
 setMethod("getBayes", "GENOME",
 

function(object,snps){

if(!snps){
  
if(!object@Pop_FSTH$calculated){
stop("First you have to calculate the F_ST.stats module")
}


   npops          <- length(object@Pop_FSTH$Populations)
   frequencies    <- lapply(object@region.stats@haplotype.counts,function(x){if(length(x)>0){return(x)}else{return(NULL)}})
   haps           <- lapply(frequencies,function(x){return(dim(x)[2])})
   
   maxhaps        <- max(unlist(haps))
   
 berendlist   <- vector("list",npops)
 berendmatrix <- matrix(NA,object@genelength,maxhaps)

 for(xx in 1:npops){
    for(yy in 1:object@genelength){
        if(length(frequencies[[yy]])>1){
         berendmatrix[yy,1:haps[[yy]]]  <-  frequencies[[yy]][xx,]
        }
    }
   berendlist[[xx]] <- berendmatrix
 } 

# Calculate hapcounts and popualtion_size

  for(xx in 1:npops){
     pop_size     <- apply(berendlist[[xx]],1,function(vv){
	         	return(sum(vv,na.rm=TRUE))		
		  })	
   hapcounts     <- apply(berendlist[[xx]],1,function(vv){
		        return(sum(!is.na(vv)))		
		  })	
  berendlist[[xx]]           <- cbind(pop_size,hapcounts,berendlist[[xx]])
  rownames(berendlist[[xx]]) <- object@region.names
  }


 return(list(LISTE=berendlist,FUNC=NULL))

}#end of change


if(snps){

GLOBAL <- new.env()
bayes_IN <- vector("list",length(object@populations))
 
 for(xx in 1:length(object@populations)){

   GLOBAL$loci <- 1 
   sss         <- sapply(object@region.data@biallelic.matrix,function(bial){

         bial <- popGetBial(object, GLOBAL$loci)
         if(length(bial)!=0 & length(object@populations)==length(object@region.data@populations[[GLOBAL$loci]])){         
	   region_pop        <- object@region.data@populations[[GLOBAL$loci]]
           sub_bial          <- bial[region_pop[[xx]],,drop=FALSE]	  
	   nullen            <- colSums(sub_bial==0,na.rm=TRUE)   
           einsen            <- colSums(sub_bial==1,na.rm=TRUE)	
           rueck      	     <- rbind(nullen,einsen) 
           GLOBAL$loci       <- GLOBAL$loci + 1
           return(rueck)
	
         }else{
	   GLOBAL$loci       <- GLOBAL$loci + 1
           return(NULL)
	 } 
	 
          })
   matt                     <- matrix(unlist(sss),ncol=2,byrow=2)
   bayes_IN[[xx]]           <- cbind(rowSums(matt),2,matt)
   colnames(bayes_IN[[xx]]) <- c("size","haps","nullen","einsen")
  }

## Create the sets 1)
GLOBAL$loci  <- 1
func  <- sapply(object@region.data@biallelic.sites,function(x){
        if(length(object@populations)==length(object@region.data@populations[[GLOBAL$loci]])){
	 q    <- rep(GLOBAL$loci,length(object@region.data@biallelic.sites[[GLOBAL$loci]]))
	 GLOBAL$loci <- GLOBAL$loci + 1; 
	 return(q)

        }else{
	   GLOBAL$loci       <- GLOBAL$loci + 1
           return(NULL)
        } 

})

rm(GLOBAL) # remove environment

## Create sets 2)
 FUNC  <- unlist(func)
#FUNC2 <- vector("list",FUNC[length(FUNC)])
#for(xx in 1:FUNC[length(FUNC)]){
#	FUNC2[[xx]] <- which(xx==FUNC)
#}

return(list(LISTE=bayes_IN,FUNC=FUNC))

}

})

 # --------------------------------------------------------------------
 # Get FST: FSTH FSTN
 # --------------------------------------------------------------------
 
 setGeneric("get.F_ST", function(object,mode=FALSE,pairwise=FALSE) standardGeneric("get.F_ST"))
 setMethod("get.F_ST", "GENOME",

 function(object,mode,pairwise){

if(!pairwise){

if(mode[1]==FALSE){

  if(!object@Pop_FSTH$calculated){stop("Statistics have to be calculated first !")}

  res           <- cbind(object@haplotype.F_ST,object@nucleotide.F_ST,object@Nei.G_ST,object@Hudson.G_ST,object@Hudson.H_ST,object@Hudson.K_ST)
  colnames(res) <- c("haplotype.F_ST","nucleotide.F_ST","Nei.G_ST","Hudson.G_ST","Hudson.H_ST","Hudson.K_ST")
  return(res)


}
 
  if(mode[1]=="haplotype"){
   res           <- cbind(object@haplotype.F_ST,object@Nei.G_ST,object@Hudson.G_ST,object@Hudson.H_ST,object@Hudson.K_ST)
   colnames(res) <- c("haplotype.F_ST","Nei.G_ST","Hudson.G_ST","Hudson.H_ST","Hudson.K_ST")
   return(res)
  }   

  if(mode[1]=="nucleotide"){
   res           <- cbind(object@nucleotide.F_ST)
   colnames(res) <- c("nucleotide.F_ST")
   return(res)
  }   

}else{# End if !pairwise

 res       <- vector("list",3)
 res[[1]]  <- t(object@nuc.F_ST.pairwise)
 res[[2]]  <- t(object@hap.F_ST.pairwise)
 res[[3]]  <- t(object@Nei.G_ST.pairwise)
 res       <- as.matrix(res)
 rownames(res) <- c("nuc.F_ST.pairwise","hap.F_ST.pairwise","Nei.G_ST.pairwise")
 colnames(res) <- "pairwise F_ST"
 return(res)


}
 })


## GETMK
setGeneric("get.MKT", function(object) standardGeneric("get.MKT"))
setMethod("get.MKT", "GENOME",

 function(object){

 return(object@MKT)

})

# ---------------------------------------------------------------------
# Helper function : get a certain biallelic matrix from a GENOME object
# ----------------------------------------------------------------------
popGetBial <- function( XX , bialmatNr )
{

# ff or not
#if(XX@big.data){

# if no biallelic sites return NULL
if(length(XX@region.data@biallelic.sites[[bialmatNr]])==0){return(NULL)}


    if(length(XX@BIG.BIAL)==0){
     if(XX@big.data){

         open(XX@region.data@biallelic.matrix[[bialmatNr]]) # open ff file

     } # da zuviele ff files open... muss nich !

     bial <- XX@region.data@biallelic.matrix[[bialmatNr]][,,drop=FALSE]

     if(XX@big.data){

        close(XX@region.data@biallelic.matrix[[bialmatNr]]) # close ff file

     }# da zuviele ff files open... muss nich !

    }else{ 
     if(length(XX@SLIDE.POS[[bialmatNr]])==0){return(NULL)} # muss nur wegen BIGMEMORY package !
     # open(XX@BIG.BIAL[[1]])
     if(length(XX@jack.knife)==0){
     bial <- XX@BIG.BIAL[[1]][,XX@SLIDE.POS[[bialmatNr]],drop=FALSE]
     }else{
     jack.positions   <- unique(unlist(XX@SLIDE.POS[-bialmatNr]))
     bial <- XX@BIG.BIAL[[1]][,jack.positions,drop=FALSE] 
     }
    }

#}else{

#  bial <- XX@region.data@biallelic.matrix[[bialmatNr]]}
###########

  if(length(bial)>1){ # sollte das UNTERE loesen :)

            #if(!is.na(bial[1])){ # is na wegen sliding window mode letztes NULL verschwindet :( FIXME muessen diese Zeilen ? include.unknown problem
	 
             return(bial)
 
            #}else{
	    #return(NULL)
            #}

  }else{return(NULL)}

}
 
# ---------------------------------------------------------------------
# Calc Neutrality (calc_freqstats)
# ----------------------------------------------------------------------
 setGeneric("neutrality.stats", function(object,new.populations=FALSE,new.outgroup=FALSE,subsites=FALSE,detail=FALSE,FAST=FALSE,do.R2=FALSE) standardGeneric("neutrality.stats"))
 setMethod("neutrality.stats", "GENOME",
 
 function(object,new.populations,new.outgroup,subsites,detail,FAST,do.R2){
 
 region.names                         <- object@region.names
 n.region.names                       <- length(region.names)

 if(object@big.data){region.names <- NULL} # because of memory space
 
 
 object@Pop_Neutrality$sites       <- "ALL"
 object@Pop_Neutrality$calculated  <- TRUE
 
 # Populations 
 if(missing(new.populations)){
 npops                             <- length(object@populations)
 object@Pop_Neutrality$Populations <- object@populations
 }else{
 npops                             <- length(new.populations)
 object@Pop_Neutrality$Populations <- new.populations
 }

 # Outgroup
 if(missing(new.outgroup)){
 object@Pop_Neutrality$Outgroup <- object@populations
 }else{
 object@Pop_Neutrality$Outgroup <- new.outgroup
 }
 
 
## bial        <- object@biallelics



 # Init
 init        <- matrix(,n.region.names,npops)
 TajD        <- init
 thetaT      <- init
 thetaS      <- init
 thetaFL     <- init
 thetaSA     <- init
 thetaTA     <- init
 thetaFW     <- init
 thetaL      <- init
 S           <- init
 FuLi_F      <- init
 FuLi_D      <- init
 R2          <- init
 HnFw        <- init
 Ez          <- init
 FS          <- init
 Strobeck    <- init
 change      <- object@region.stats

 #--- for Coalescent simulation
 Pop_Neutrality <- vector("list",n.region.names)
 #-----------------------------

 # Names
 nam <- paste("pop",1:npops)
 #------------------------------

 rownames(TajD)     <- region.names
 colnames(TajD)     <- nam
 rownames(thetaT)   <- region.names
 colnames(thetaT)   <- nam
 rownames(thetaS)   <- region.names
 colnames(thetaS)   <- nam
 rownames(thetaFL)  <- region.names
 colnames(thetaFL)  <- nam
 rownames(thetaSA)  <- region.names
 colnames(thetaSA)  <- nam 
 rownames(thetaTA)  <- region.names
 colnames(thetaTA)  <- nam
 rownames(thetaFW)  <- region.names
 colnames(thetaFW)  <- nam
 rownames(thetaL)   <- region.names
 colnames(thetaL)   <- nam
 rownames(S)        <- region.names
 colnames(S)        <- nam
 rownames(FuLi_F)   <- region.names
 colnames(FuLi_F)   <- nam
 rownames(FuLi_D)   <- region.names
 colnames(FuLi_D)   <- nam
 rownames(R2)       <- region.names
 colnames(R2)       <- nam
 rownames(HnFw)     <- region.names
 colnames(HnFw)     <- nam
 rownames(Ez)       <- region.names
 colnames(Ez)       <- nam
 rownames(FS)       <- region.names
 colnames(FS)       <- nam
 rownames(Strobeck) <- region.names
 colnames(Strobeck) <- nam
 
  # Populations
  if(!missing(new.populations)){
   NEWPOP <- TRUE
   populations <- vector("list",npops)
  }else{
   NEWPOP <- FALSE
  } 
  
  # Outgroup
  if(!missing(new.outgroup)){
   NEWOUT <- TRUE
  }else{
   NEWOUT <- FALSE
  } 
  


## PROGRESS #########################
 progr <- progressBar()
#####################################

for(xx in 1:n.region.names){


### if Subsites ----------------------------------
bial <- popGetBial(object,xx)

if(length(bial)==0){next} # when no biallelic positions in the region


if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
   #object@Pop_Neutrality$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_Neutrality$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
   #object@Pop_Neutrality$sites <- "synonymous"
}

if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
   #object@Pop_Neutrality$sites <- "nonsynonymous"
}

if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   bial          <- bial[,intron,drop=FALSE]
  # object@Pop_Neutrality$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
   #object@Pop_Neutrality$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
   #object@Pop_Neutrality$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]]))
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
   # object@Pop_Neutrality$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]])
   bial             <- bial[,gene,drop=FALSE]
  # object@Pop_Neutrality$sites <- "gene"
}

} # End of if subsites

############### ---------------------------------

    
  if(length(bial)!=0){ # if a biallelic position exists  
    
 
   ## Populations
   if(NEWPOP){ # wenn neu Populationen definiert
    
       for(yy in 1:npops){
           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]
           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }      
       }
           
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   
       if(length(populations)==0){next} # Keine Population vorhanden

    }else{
     populations  <- object@region.data@populations[[xx]] # if there is no new population
    }

  ## Outgroup
  if(NEWOUT){

	if(is.character(new.outgroup)){
           outgroup <- match(new.outgroup,rownames(bial)) 
           naids    <- which(!is.na(outgroup))
           outgroup <- outgroup[naids]  
        }else{
           outgroup <- new.outgroup
           ids      <- which(outgroup>dim(bial)[1])
           if(length(ids)>0){outgroup <- outgroup[-ids]}   
        }
        
        if(length(outgroup)==0){outgroup <- FALSE}

  }else{
    outgroup <- object@region.data@outgroup[[xx]]
  }
  # --------------------------------------
      
    
    # important for Coalescent Simulation
    # change@Pop_Neutrality[[xx]] <- list(Populations=populations,Outgroup=outgroup)
      Pop_Neutrality[[xx]] <- list(Populations=populations,Outgroup=outgroup)
 
    # important for Coalescent Simulation
    #change@Pop_Neutrality[[xx]]$Outgroup    <- outgroup

   # NON FAST C
   if(!FAST){ 
    if(object@Pop_Slide$calculated | subsites[1]!=FALSE | object@snp.data){
    res          <- calc_freqstats(bial,populations=populations,outgroup=outgroup)
    }else{
    res          <- calc_freqstats(bial,populations=populations,outgroup=outgroup,
                    data=list(n.nucleotides=object@region.data@n.nucleotides[[xx]],n.valid.sites=object@n.valid.sites[xx]))
    }
   } 
   # FAST C
   if(FAST){ 
    if(object@Pop_Slide$calculated | subsites[1]!=FALSE | object@snp.data){
    res          <- calc_freqstats_FAST(bial,populations=populations,outgroup=outgroup)
    }else{
    res          <- calc_freqstats_FAST(bial,populations=populations,outgroup=outgroup,
                    data=list(n.nucleotides=object@region.data@n.nucleotides[[xx]],n.valid.sites=object@n.valid.sites[xx]))
    }
   } 

    if(NEWPOP)  {if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}} # nur die Populationen, die existieren
    if(!NEWPOP) {if(length(object@region.data@popmissing[[xx]])!=0){popmissing <- object@region.data@popmissing[[xx]];respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}
    
    # -----------------------
    #   fill detailed Slots
    # -----------------------
    
    if(detail){

     
     res2                              <- jointfreqdist(bial,populations=populations,outgroup=outgroup)
     change@minor.allele.freqs[[xx]]   <- res2$jfd
     res2                              <- calc_FS(bial,populations,res$THETA["thetaT",])
     FS[xx,respop]                     <- res2$FS
     Strobeck[xx,respop]               <- res2$Strobeck	


     # misdis                  <- mismatch(bial,populations=populations,res$THETA["thetaT",])    
     # change[[xx]]@MISDIS     <- misdis
     #;cat("2");
     # object@GEN[[xx]]@thetaS  <- res$THETA["thetaS",]	;cat("1");
     # object@GEN[[xx]]@thetaT  <- res$THETA["thetaT",]	;cat("1");
     # object@GEN[[xx]]@thetaFL <- res$THETA["thetaFL",]	;cat("3");
     # object@GEN[[xx]]@thetaSA <- res$THETA["thetaSA",]	;cat("4");
     # object@GEN[[xx]]@thetaTA <- res$THETA["thetaTA",]	;cat("5");
     # object@GEN[[xx]]@thetaFW <- res$THETA["thetaFW",]	;cat("6");
     # object@GEN[[xx]]@thetaL  <- res$THETA["thetaL",]	;cat("7");
     # change[[xx]]@THETA       <- res$THETA
     
     }
    # -----------------------#
   
    #;cat("1");
    #change@sfreq[[xx]]  <- res$FREQ

    TajD[xx,respop]     <- res$taj_D 
    thetaT[xx,respop]   <- res$THETA[3,]
    thetaS[xx,respop]   <- res$THETA[2,]
    thetaFL[xx,respop]  <- res$THETA[4,]
    thetaSA[xx,respop]  <- res$THETA[5,]
    thetaTA[xx,respop]  <- res$THETA[6,]
    thetaFW[xx,respop]  <- res$THETA[7,]
    thetaL[xx,respop]   <- res$THETA[8,]
    
    S   [xx,respop]     <- res$THETA[1,]
    FuLi_F[xx,respop]   <- res$FuLi_F
    FuLi_D[xx,respop]   <- res$FuLi_D
    HnFw[xx,respop]     <- res$HnFw
    Ez[xx,respop]       <- res$Ez
    
    if(do.R2){
     r2                  <- calcR2(bial,populations,res$THETA["thetaT",],res$THETA["S",])
     R2[xx,respop]       <- r2
    }   

   
  # PROGRESS #######################################################
     progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
 


}
 
 change@Pop_Neutrality        <- Pop_Neutrality # important for Coalescent Simulation
 object@region.stats          <- change
 object@Tajima.D              <- TajD
 object@n.segregating.sites   <- S
 object@Fu.Li.F               <- FuLi_F
 object@Fu.Li.D               <- FuLi_D
 object@Rozas.R_2             <- R2
 object@theta_Tajima          <- thetaT
 object@theta_Watterson       <- thetaS
 object@theta_Fu.Li           <- thetaFL
 object@theta_Achaz.Watterson <- thetaSA
 object@theta_Achaz.Tajima    <- thetaTA
 object@theta_Fay.Wu          <- thetaFW
 object@theta_Zeng            <- thetaL
 object@Fay.Wu.H              <- HnFw
 object@Zeng.E                <- Ez
 object@Fu.F_S                <- FS
 object@Strobeck.S            <- Strobeck

  return(object)
  
 })
 
# ------------------------------------------------------------
# popLinkage  (linkdisequ)
# ------------------------------------------------------------

setGeneric("linkage.stats", function(object,new.populations=FALSE,subsites=FALSE,detail=FALSE, do.ZnS=TRUE, do.WALL=TRUE) standardGeneric("linkage.stats"))
 setMethod("linkage.stats", "GENOME",
 function(object,new.populations,subsites,detail,do.ZnS,do.WALL){

# if(include.unknown){detail <-TRUE} #FIXME

 region.names                       <- object@region.names
 n.region.names                     <- length(region.names)
 if(object@big.data){region.names <- NULL} # because of memory space
 
 
 object@Pop_Linkage$sites        <- "ALL"
 object@Pop_Linkage$calculated   <- TRUE
 
 if(missing(new.populations)){
 npops                           <- length(object@populations)
 object@Pop_Linkage$Populations  <- object@populations
 }else{
 npops                          <- length(new.populations)
 object@Pop_Linkage$Populations <- new.populations
 
 }
 
# bial        <- object@biallelics
 
# Init
 init        <- matrix(,n.region.names,npops)
 Zns         <- init
 ZA          <- init
 ZZ          <- init
 WALLB       <- init
 WALLQ       <- init
 linkage.disequilibrium <- vector("list",n.region.names)
 change      <- object@region.stats
 Pop_Linkage <- vector("list",n.region.names) # important for Coalescent simulation


# Names
 nam               <- paste("pop",1:npops)
 rownames(Zns)     <- region.names
 colnames(Zns)     <- nam
 rownames(ZA)      <- region.names
 colnames(ZA)      <- nam
 rownames(ZZ)      <- region.names
 colnames(ZZ)      <- nam
 rownames(WALLB)   <- region.names
 colnames(WALLB)   <- nam
 rownames(WALLQ)   <- region.names
 colnames(WALLQ)   <- nam

  if(!missing(new.populations)){
   NEWPOP <- TRUE
   populations <- vector("list",npops)
  }else{
   NEWPOP <- FALSE
  } 

## PROGRESS #########################
 progr <- progressBar()
#####################################


for(xx in 1:n.region.names){


### if Subsites ----------------------------------
bial <- popGetBial(object,xx)


if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
   #object@Pop_Linkage$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
   #object@Pop_Linkage$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
   #object@Pop_Linkage$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
   #object@Pop_Linkage$sites <- "nonsynonymous"
}


if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
   #object@Pop_Linkage$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
   #object@Pop_Linkage$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
   #object@Pop_Linkage$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]]))
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
   #object@Pop_Linkage$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]])
   bial             <- bial[,gene,drop=FALSE]
   #object@Pop_Linkage$sites <- "gene"
}
}# End if subsites
############### ---------------------------------



  if(length(bial)!=0){ # if a biallelic position exists  
       
    if(NEWPOP){ # wenn neu Populationen definiert
       for(yy in 1:npops){
           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]
           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }   
       }
       
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   
       if(length(populations)==0){next} # Keine Population vorhanden
       
    }else{
     populations <- object@region.data@populations[[xx]] # if there is no new population
    }

    
    if(NEWPOP)  {if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}} # nur die Populationen, die existieren
    if(!NEWPOP) {if(length(object@region.data@popmissing[[xx]])!=0){popmissing <- object@region.data@popmissing[[xx]];respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}
    

     # important for Coalescent Simulation
     # change@Pop_Linkage[[xx]] <- list(Populations=populations,Outgroup=NULL)
     Pop_Linkage[[xx]]        <- list(Populations=populations,Outgroup=NULL)
     # ------------------- fill detail slots


  if(do.ZnS){

    if(!detail){
    res                          <- linkdisequ_FAST(bial,populations)
    Zns[xx,respop]               <- res$Zns
    ZA [xx,respop]               <- res$ZA
    ZZ [xx,respop]               <- res$ZZ
    }   

    if(detail){
     res                         <- linkdisequ(bial,populations)
     linkage.disequilibrium[[xx]]   <- res$res
     Zns[xx,respop]              <- res$Zns
     ZA [xx,respop]              <- res$ZA
     ZZ [xx,respop]              <- res$ZZ
    }
    # -------------------
  }  
       
    if(do.WALL){
    res              <- wall99bq(bial,populations)
    WALLB[xx,respop] <- res$B 
    WALLQ[xx,respop] <- res$Q
    }

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}
 change@linkage.disequilibrium <- linkage.disequilibrium
 change@Pop_Linkage  <- Pop_Linkage # important for Coalescent Simulation
 object@region.stats <- change 
 object@Kelly.Z_nS   <- Zns
 object@Rozas.ZA     <- ZA
 object@Rozas.ZZ     <- ZZ
 object@Wall.B       <- WALLB
 object@Wall.Q       <- WALLQ

  
  return(object)
 })
 
 # --------------------------------------------------------------------
 # popFSTN: FST (Nucleotide)  --->  fstcalc < ---
 # --------------------------------------------------------------------
 
 setGeneric("popFSTN", function(object,new.populations="list",subsites=FALSE,detail=TRUE,mode="nucleotide") standardGeneric("popFSTN"))
 setMethod("popFSTN", "GENOME",

 function(object,new.populations,subsites,detail,mode){
  
  region.names     <- object@region.names
  n.region.names  <- length(region.names)
  if(object@big.data){region.names <- NULL} # because of memory space
 
  
  #bial        <- object@biallelics
  
  object@Pop_FSTN$sites        <- "ALL"
  object@Pop_FSTN$calculated   <- TRUE
 
 # Populations 
  if(missing(new.populations)){
   npops       <- length(object@populations)
   object@Pop_FSTN$Populations <- object@populations
  }else{
   npops           <- length(new.populations)
   object@Pop_FSTN$Populations <- new.populations
  }
  
# Outgroup
# if(missing(new.outgroup)){
# object@Pop_FSTN$Outgroup <- object@populations
# }else{
# object@Pop_FSTN$Outgroup <- new.outgroup
# }
 

  
  # Get the names ----------for pairwaise comparison --------------------------------------------------
  #############################################################################
  
if(npops>1){
 #if(outgroup[1]!=F){
 # poppairs <- choose(npops+1,2) # Outgroup is included !!
 # pairs    <- combn(1:(npops+1),2)
 #}else{
  poppairs <- choose(npops,2)   # Outgroup is not included !!
  pairs    <- combn(1:(npops),2)
 #} 
 
#### --- Names of population pairs --- #### 
 nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
 if(dim(pairs)[2]>1){ # more than 2 Populations
  for(xx in 2:dim(pairs)[2]){
    m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
  } 
 }#END if
}# End npops > 1
else{poppairs <- 1;nn <- "pop1"} 
##### ------------------------------ ####------------------------------------------------ 
#########################################################################################  
  
   
  # INIT
  init  <- matrix(0,n.region.names,npops)
  init1 <- matrix(0,length(nn),n.region.names)
  init2 <- matrix(0,n.region.names,1)
  
  nuc.F_ST.vs.all   <- init
  PIW               <- matrix(0,n.region.names,npops)
  PIA               <- init1
  nuc.F_ST.pairwise <- init1
  FSTN              <- init2
  nucleotide.diversity <- vector("list",n.region.names) # region.stats
  change            <- object@region.stats
 
  popnames          <- paste("pop",1:npops) 
  # Names
  rownames(FSTN)    <- region.names
  colnames(FSTN)    <- "FST (Nucleotide)"
  
  rownames(PIA)     <- nn
  colnames(PIA)     <- region.names
  rownames(nuc.F_ST.pairwise) <- nn
  colnames(nuc.F_ST.pairwise) <- region.names
  rownames(PIW)      <- region.names
  colnames(PIW)      <- popnames 
  rownames(nuc.F_ST.vs.all)      <- region.names
  colnames(nuc.F_ST.vs.all)      <- popnames

  # Populations
  if(!missing(new.populations)){
   NEWPOP <- TRUE
   populations <- vector("list",npops)
  }else{
   NEWPOP <- FALSE
  } 

  # Outgroup
 # if(!missing(new.outgroup)){
 #  NEWOUT <- TRUE
 # }else{
 #  NEWOUT <- FALSE
 # } 

Pout                    <- TRUE
if(mode[1]=="ALL"){Pout <- TRUE}

if(Pout){
## PROGRESS #########################
 progr <- progressBar()
#####################################
}


for(xx in 1:n.region.names){

### if Subsites ----------------------------------

bial <- popGetBial(object,xx)

if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
   #object@Pop_FSTN$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_FSTN$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
  # object@Pop_FSTN$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
  # object@Pop_FSTN$sites <- "nonsynonymous"
}

if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
  # object@Pop_Linkage$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
  # object@Pop_FSTN$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
  # object@Pop_FSTN$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]])==TRUE)
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
  # object@Pop_FSTN$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]]==TRUE)
   bial             <- bial[,gene,drop=FALSE]
  # object@Pop_FSTN$sites <- "gene"
}

if(subsites=="intergenic"){
  intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   if(length(intron)==0){
     intron <- !object@region.data@ExonSNPS[[xx]]	  
   }
  utr            <- object@region.data@UTRSNPS[[xx]]
  exon           <- object@region.data@ExonSNPS[[xx]]
  gene           <- object@region.data@GeneSNPS[[xx]]
  coding         <- !is.na(object@region.data@synonymous[[xx]])  

  inter          <- !(intron|utr|exon|gene|coding)
  bial           <- bial[,inter,drop=FALSE]
 # object@Pop_FSTN$sites <- "intergenic"
}
}# end if subsites

############### ---------------------------------

 if(length(bial)!=0){ # if a biallelic position exists
    
    # Populations
    if(NEWPOP){ # wenn neu Populationen definiert
    
       for(yy in 1:npops){
           if(is.character(new.populations[[yy]])){
             # populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))

              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]
           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }   
       }
       
       #----------------------#
       temp         <- delNULLpop(populations) # L\F6sche nicht vorhandene Populationen
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   
       if(length(populations)==0){next} # Keine Population vorhanden

       

    }else{
     populations <- object@region.data@populations[[xx]]

    }
     
  # Outgroup
  #if(NEWOUT){

#	if(is.character(new.outgroup)){
#           outgroup <- match(new.outgroup,rownames(object@biallelics[[xx]])) 
#           naids    <- which(!is.na(outgroup))
#           outgroup <- outgroup[naids]  
#        }else{
#           outgroup <- new.outgroup
#           ids      <- which(outgroup>dim(bial)[1])
#           if(length(ids)>0){outgroup <- outgroup[-ids]}   
#        }
#        
#        if(length(outgroup)==0){outgroup <- FALSE}
#
#  }else{
#    outgroup <- object@DATA[[xx]]@outgroup
#  }
  # --------------------------------------


    # important for Coalescent Simulation wird zur Zeit eh nich benutzt !
    # change@Pop_FSTN[[xx]]$Populations <- populations
    # change[[xx]]@Pop_FSTN$Outgroup    <- outgroup


    if(!object@Pop_Slide$calculated){
    data.list              <- list(n.nucleotides=object@region.data@n.nucleotides[[xx]],n.valid.sites=object@n.valid.sites[xx],
                              transitions=object@region.data@transitions[[xx]],biallelic.compositions=object@region.data@biallelic.compositions[[xx]])
    }else{
    data.list              <- list(n.nucleotides=NULL,n.valid.sites=NULL,
                              transitions=object@region.data@transitions[[xx]],biallelic.compositions=NULL)
    }
    
    res                    <- fstcalc(bial,populations,data=data.list,outgroup=FALSE)
    if(NEWPOP) {temp       <- checkpoppairs(npops,popmissing,pairs,nn)} # welche populationen wurden \FCberhaupt berechnet
    if(!NEWPOP){temp       <- checkpoppairs(npops,object@region.data@popmissing[[xx]],pairs,nn)} 
   
    respop     <- temp$respop
    respairpop <- temp$respairpop
   
   # Some detailed statistic --------#
   if(detail){ 
      # change[[xx]]@FSTPAIR                <- res$FSTPAIR 
      # change[[xx]]@FST1ALL                <- res$FST1ALL
      nucleotide.diversity[[xx]]            <- res$PIA   
     # change[[xx]]@SV1         <- res$SV1
     # change[[xx]]@SV2         <- res$SV2
    
   }
   # --------------------------------#  

              
    FSTN[xx]                        <- res$FSTALL
    PIA[respairpop,xx]              <- res$PIA2
    nuc.F_ST.pairwise[respairpop,xx]<- res$FSTPAIR2
    PIW[xx,respop]                  <- res$PIW
    nuc.F_ST.vs.all[xx,respop]      <- res$FST1ALL

if(Pout){
    # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
    ###################################################################
}

 }
}  

 change@nucleotide.diversity  <- nucleotide.diversity
 object@region.stats          <- change  
 object@nucleotide.F_ST       <- FSTN
 object@nuc.diversity.between <- PIA
 object@nuc.diversity.within  <- PIW
 object@nuc.F_ST.vs.all       <- nuc.F_ST.vs.all
 object@nuc.F_ST.pairwise     <- nuc.F_ST.pairwise

 return(object)
 })


 # --------------------------------------------------------------------
 # popFSTH: FST (Haplotype) -----> calc_hwhafsth <------------------
 # --------------------------------------------------------------------

 #setGeneric("getFSTH", function(object) standardGeneric("getFSTH"))
 #setMethod("getFSTH", "GENOME", function(object) return(object@FSTH))

 setGeneric("F_ST.stats", function(object,new.populations=FALSE,subsites=FALSE,detail=TRUE,mode="ALL",only.haplotype.counts=FALSE,FAST=FALSE) standardGeneric("F_ST.stats"))
 setMethod("F_ST.stats","GENOME",function(object,new.populations,subsites,detail,mode,only.haplotype.counts,FAST){
  
 # mode nur wegen Progress Balken !

 # hier die FST abgespeckte FAST Version !

 if(FAST){

    if(!missing(new.populations)){
         
         object  <- SNPFST(object,new.populations,subsites=subsites,detail=detail,mode=mode,only.haplotype.counts=only.haplotype.counts)

    }else{
        
         object  <- SNPFST(object,subsites=subsites,detail=detail,mode=mode,only.haplotype.counts=only.haplotype.counts)
    }

    
    
    return(object)
  
 }
# ----------------------------------------



 # also calculate popFSTN
 if(mode[1]=="ALL"){

   cat("nucleotide \n")
   if(!missing(new.populations)){
         object <- popFSTN(object,new.populations,subsites=subsites,detail=detail,mode=mode)
   }else{
         object <- popFSTN(object,subsites=subsites,detail=detail,mode=mode)
   }
   cat("\n")
   cat("haplotype \n")

 }


 if(mode[1]=="nucleotide"){
   
    if(!missing(new.populations)){
         object <- popFSTN(object,new.populations,subsites=subsites,detail=detail,mode=mode)
    }else{object <- popFSTN(object,subsites=subsites,detail=detail,mode=mode)
    }
 
 return(object)
 }
 #
 
  region.names                  <- object@region.names
   n.region.names               <- length(region.names)
  if(object@big.data){region.names <- NULL} # because of memory space
 
  object@Pop_FSTH$sites        <- "ALL"
  object@Pop_FSTH$calculated   <- TRUE
   
  if(!missing(new.populations)){
    NEWPOP <- TRUE
    populations <- vector("list",length(new.populations))
    npops       <- length(populations)            # Wenn mehr Pops definiert werden
    object@Pop_FSTH$Populations  <- new.populations
  }else{
    NEWPOP <- FALSE
    npops                       <- length(object@populations)     # alte Anzahl der Populationen
    object@Pop_FSTH$Populations <- object@populations
  }
 
 #########################################
 # INIT
 #########################################
  
   # Get the names ----------for pairwaise comparison --------------------------------------------------
  #############################################################################
  
if(npops>1){
 #if(outgroup[1]!=F){
 # poppairs <- choose(npops+1,2) # Outgroup is included !!
 # pairs    <- combn(1:(npops+1),2)
 #}else{
  poppairs <- choose(npops,2)   # Outgroup is not included !!
  pairs    <- combn(1:(npops),2)
 #} 
 
#### --- Names of population pairs --- #### 
 nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
 if(dim(pairs)[2]>1){ # more than 2 Populations
  for(xx in 2:dim(pairs)[2]){
    m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
  } 
 }#END if
}# End npops > 1
else{poppairs <- 1;nn <- "pop1"} 
##### ------------------------------ ####------------------------------------------------ 
#########################################################################################  
  
  init3  <- matrix(0,length(nn),n.region.names)
  nam    <- paste("pop",1:npops)
  init1  <- matrix(0,n.region.names,npops)
  init2  <- matrix(0,n.region.names,1)
  
  FSTN2  <- init2
  FSTH   <- init2
  GST    <- init2 
  KST    <- init2
  GSTH   <- init2
  HST    <- init2
  Pi     <- init1
  hapw   <- init1
  hap.F_ST.vs.all       <- init1
  hap.F_ST.pairwise     <- init3
  Nei.G_ST.pairwise     <- init3
  hap.diversity.between <- init3
  haplotype.diversity   <- vector("list",n.region.names) # region stats
  haplotype.counts      <- vector("list",n.region.names) # region stats
 
  change    <- object@region.stats
  Pop_FSTH  <- vector("list",n.region.names)
  
#--------------------------------------------------

  # Names ----------------------------------------
  rownames(FSTH)     <- region.names
  colnames(FSTH)     <- "FST (Haplotype)"
  rownames(GST)      <- region.names
  colnames(GST)      <- "GST (Nei)"
  rownames(GSTH)     <- region.names
  colnames(GSTH)     <- "GST (Hudson)"
  rownames(HST)      <- region.names
  colnames(HST)      <- "HST (Hudson)"
  rownames(KST)      <- region.names
  colnames(KST)      <- "KST (Hudson)"
  rownames(Pi)       <- region.names
  colnames(Pi)       <- nam
  rownames(hapw)     <- region.names


  colnames(hapw)     <- nam
  rownames(hap.F_ST.vs.all)     <- region.names
  colnames(hap.F_ST.vs.all)     <- nam
  rownames(FSTN2)    <- region.names
  colnames(FSTN2)    <- "FST (Nucleotide)"
  rownames(hap.F_ST.pairwise) <- nn
  colnames(hap.F_ST.pairwise) <- region.names
  rownames(Nei.G_ST.pairwise) <- nn
  colnames(Nei.G_ST.pairwise) <- region.names
  rownames(hap.diversity.between) <- nn
  colnames(hap.diversity.between) <- region.names
  # ----------------------------------------------


## PROGRESS #########################
 progr <- progressBar()
#####################################



   
for(xx in 1:n.region.names){

### if Subsites ----------------------------------

bial <- popGetBial(object,xx)

if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
  # object@Pop_FSTH$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_FSTH$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
  # object@Pop_FSTH$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
  # object@Pop_FSTH$sites <- "nonsynonymous"
}

if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
  # object@Pop_Linkage$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
  # object@Pop_FSTH$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
  # object@Pop_FSTH$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]])==TRUE)
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
  # object@Pop_FSTH$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]]==TRUE)
   bial             <- bial[,gene,drop=FALSE]
  # object@Pop_FSTH$sites <- "gene"
}

if(subsites=="intergenic"){
  intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   if(length(intron)==0){
     intron <- !object@region.data@ExonSNPS[[xx]]	  
   }

  utr            <- object@region.data@UTRSNPS[[xx]]
  exon           <- object@region.data@ExonSNPS[[xx]]
  gene           <- object@region.data@GeneSNPS[[xx]]
  coding         <- !is.na(object@region.data@synonymous[[xx]])  

  inter          <- !(intron|utr|exon|gene|coding)
  bial           <- bial[,inter,drop=FALSE]
  #object@Pop_FSTH$sites <- "intergenic"
}
}# End if subsites
############### ---------------------------------


 if(length(bial)!=0){ # if a biallelic position exists
     
     if(NEWPOP){ # Wenn eine andere Population definiert !
          
       for(yy in 1:npops){       
           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]             
           }else{ # numeric values
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}
           }                      
       }     
       
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   
       if(length(populations)==0){next} # Keine Population vorhanden
       
        
    }else{populations <- object@region.data@populations[[xx]]} # Wenn keine neue Pop definiert !
    
    # Calculations
    # mat         <- bial[[xx]][unique(unlist(populations)),,drop=FALSE]
    res           <- calc_hwhafsth(bial,populations,only.haplotype.counts=only.haplotype.counts)




    #if(NEWPOP)  {if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}} # nur die Populationen, die existieren
    #if(!NEWPOP) {if(length(object@DATA[[xx]]@popmissing)!=0){popmissing<- object@DATA[[xx]]@popmissing;respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}
    
    if(NEWPOP) {temp       <- checkpoppairs(npops,popmissing,pairs,nn)} # welche populationen wurden \FCberhaupt berechnet
    if(!NEWPOP){temp       <- checkpoppairs(npops,object@region.data@popmissing[[xx]],pairs,nn)} 
   
    # important for Coalescent Simulation
    # change@Pop_FSTH[[xx]] <- list(Populations=populations,Outgroup=NULL)
      Pop_FSTH[[xx]]        <-  list(Populations=populations,Outgroup=NULL)   

    respop     <- temp$respop
    respairpop <- temp$respairpop
    
    
   # fill detailed Slots --------------------------------#
    if(detail){
    
     haplotype.counts[[xx]]              <- res$sfreqh
     # change[[xx]]@FSTHMATRIX           <- res$fsthmatrix
     # change[[xx]]@Gst_Hudson_pair      <- res$Gst_Hudson_pair
     # change[[xx]]@HSTpair              <- res$HSTpair
     # change[[xx]]@GstMATRIX            <- res$Gstmatrix                
     haplotype.diversity[[xx]]           <- res$hapamatrix
    # change[[xx]]@KST                   <- res$KST
   
   } 
  # ----------------------------------------------------# 
    
    hap.F_ST.pairwise[respairpop,xx]     <- res$fsth
    Nei.G_ST.pairwise[respairpop,xx]     <- res$Gst
    hap.diversity.between[respairpop,xx] <- res$hapa
    FSTH[xx]          <- res$fsthALL
    FSTN2[xx]         <- res$fstnALL
    Pi[xx,respop]     <- res$PIW_nei
    hapw[xx,respop]   <- res$hapw
    GST[xx]           <- res$GstAll
    KST[xx]           <- res$KST
    # sfreqh[[xx]]      <- res$sfreqh
    GSTH[xx]          <- res$Gst_Hudson
    HST[xx]           <- res$HST
    hap.F_ST.vs.all[xx,respop] <- res$fsthALL
    

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}
 
 change@haplotype.diversity   <- haplotype.diversity
 change@haplotype.counts      <- haplotype.counts
 change@Pop_FSTH              <- Pop_FSTH
 object@region.stats          <- change
 object@hap.F_ST.pairwise     <- hap.F_ST.pairwise
 object@Nei.G_ST.pairwise     <- Nei.G_ST.pairwise
 object@hap.diversity.between <- hap.diversity.between
 object@hap.diversity.within  <- hapw
 object@haplotype.F_ST        <- FSTH
 object@nucleotide.F_ST2      <- FSTN2 
 object@Pi                    <- Pi
 object@Nei.G_ST              <- GST
 object@Hudson.K_ST           <- KST
 #object@sfreqh  <- sfreqh
 object@Hudson.G_ST           <- GSTH
 object@Hudson.H_ST           <- HST
 object@hap.F_ST.vs.all       <- hap.F_ST.vs.all
 

 return(object)
 })




###################### SNN #######################################################################

setGeneric("F_ST.stats.2", function(object,new.populations="list",subsites=FALSE,snn=TRUE,Phi_ST=FALSE) standardGeneric("F_ST.stats.2"))
 setMethod("F_ST.stats.2", "GENOME",

 function(object,new.populations,subsites,snn,Phi_ST){

  region.names                     <- object@region.names  
  n.region.names                   <- length(region.names)
  if(object@big.data){region.names <- NULL} # because of memory space

  if(!missing(new.populations)){
    NEWPOP <- TRUE
    populations <- vector("list",length(new.populations))
    npops       <- length(populations)                            # Wenn mehr Pops definiert werden
    
  }else{
    NEWPOP <- FALSE
    npops                       <- length(object@populations)     # alte Anzahl der Populationen
   
  }
 
 #########################################
 # INIT
 #########################################
  
   # Get the names ----------for pairwaise comparison --------------------------------------------------
  #############################################################################
  
if(npops>1){
 #if(outgroup[1]!=F){
 # poppairs <- choose(npops+1,2) # Outgroup is included !!
 # pairs    <- combn(1:(npops+1),2)
 #}else{
  poppairs <- choose(npops,2)   # Outgroup is not included !!
  pairs    <- combn(1:(npops),2)
 #} 
 
#### --- Names of population pairs --- #### 
 nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
 if(dim(pairs)[2]>1){ # more than 2 Populations
  for(xx in 2:dim(pairs)[2]){
    m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
  } 
 }#END if
}# End npops > 1
else{poppairs <- 1;nn <- "pop1"} 
##### ------------------------------ ####------------------------------------------------ 
#########################################################################################  
  
  init   <- matrix(,n.region.names,1)
  Snn    <- init
  PhiST  <- init
  Pi     <- init 
#--------------------------------------------------

  # Names ----------------------------------------
  rownames(Snn)     <- region.names
  colnames(Snn)     <- "Snn"
  rownames(PhiST)   <- region.names
  colnames(PhiST)   <- "Phi_ST"
  rownames(Pi)      <- region.names
  colnames(Pi)      <- "Pi"
  # ----------------------------------------------


## PROGRESS #########################
 progr <- progressBar()
#####################################
   
for(xx in 1:n.region.names){

bial <- popGetBial(object,xx)

# if subsites
if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
   #object@Pop_FSTN$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_FSTN$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
  # object@Pop_FSTN$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
  # object@Pop_FSTN$sites <- "nonsynonymous"
}

if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   #if(length(intron)==0){
   #       intron <- object@region.data@GeneSNPS[[xx]] & !object@region.data@ExonSNPS[[xx]]	  
   #}
   bial          <- bial[,intron,drop=FALSE]
  # object@Pop_Linkage$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
  # object@Pop_FSTN$sites <- "utr"
}

if(subsites=="exon" & length(bial!=0)){
   exon           <- which(object@region.data@ExonSNPS[[xx]]==TRUE)
   bial           <- bial[,exon,drop=FALSE]
  # object@Pop_FSTN$sites <- "exon"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]])==TRUE)
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
  # object@Pop_FSTN$sites <- "coding"
}

if(subsites=="gene" & length(bial!=0)){
   gene             <- which(object@region.data@GeneSNPS[[xx]]==TRUE)
   bial             <- bial[,gene,drop=FALSE]
  # object@Pop_FSTN$sites <- "gene"
}

if(subsites=="intergenic"){
  intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   if(length(intron)==0){
     intron <- !object@region.data@ExonSNPS[[xx]]	  
   }
  utr            <- object@region.data@UTRSNPS[[xx]]
  exon           <- object@region.data@ExonSNPS[[xx]]
  gene           <- object@region.data@GeneSNPS[[xx]]
  coding         <- !is.na(object@region.data@synonymous[[xx]])  

  inter          <- !(intron|utr|exon|gene|coding)
  bial           <- bial[,inter,drop=FALSE]
 # object@Pop_FSTN$sites <- "intergenic"
}
}# end if subsites


 if(length(bial)!=0){ # if a biallelic position exists
     
     if(NEWPOP){ # Wenn eine andere Population definiert !       
       for(yy in 1:npops){       
           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]             
           }else{ # numeric values
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}
           }                      
       }     
       
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   
       if(length(populations)==0){next} # Keine Population vorhanden
       
        
    }else{populations <- object@region.data@populations[[xx]]} # Wenn keine neue Pop definiert !
    
    # Calculations
    if(snn){
     res                    <- snn(bial,populations)  
     Snn[xx]                <- res
    }
    if(Phi_ST){
     res                    <- calc_phi_st(bial, populations)
     PhiST[xx]              <- res$phi_ST
     Pi   [xx]              <- res$whole.pi
    }
    

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}
 
 object@Pi                   <- Pi
 object@Phi_ST               <- PhiST
 object@Hudson.Snn           <- Snn
 

 return(object)

 })


##################################################################################################




 #########################################################################
 # MKT ####################################################################
 setGeneric("MKT", function(object,new.populations=FALSE, do.fisher.test=FALSE) standardGeneric("MKT"))
 setMethod("MKT", "GENOME",

 function(object,new.populations, do.fisher.test){
  
  region.names   <- object@region.names
  n.region.names <- length(object@region.names)
  if(object@big.data){region.names <- NULL} # because of memory space

  object@Pop_MK$sites        <- "ALL"
  object@Pop_MK$calculated   <- TRUE
  
  if(missing(new.populations)){
   npops                     <- length(object@populations)
   object@Pop_MK$Populations <- object@populations
  }else{
   npops                     <- length(new.populations)
   object@Pop_MK$Populations <- new.populations
  }
  
  # Get the names ----------for pairwaise comparison --------------------------------------------------
  #############################################################################
  
if(npops>1){
 #if(outgroup[1]!=F){
 # poppairs <- choose(npops+1,2) # Outgroup is included !!
 # pairs    <- combn(1:(npops+1),2)
 #}else{
  poppairs <- choose(npops,2)   # Outgroup is not included !!
  pairs    <- combn(1:(npops),2)
 #} 
 
#### --- Names of population pairs --- ################################################# 
 nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
 if(dim(pairs)[2]>1){ # more than 2 Populations
  for(xx in 2:dim(pairs)[2]){
    m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
    nn <- c(nn,m)
  } 
 }#END if
}# End npops > 1
else{stop("You have to define more than one population to calculate the MK Test")} 
##### ------------------------------ ####------------------------------------------------ 
#########################################################################################  
  
   
  # INIT
  MKmulti  <- vector("list",n.region.names)
  MKTWO    <- matrix(,n.region.names,7)
  # Pn/Ps/Dn/Ds
 
  # Names
  MKnames          <- c("P_nonsyn","P_syn","D_nonsyn","D_syn","neutrality.index","alpha","fisher.P.value")
  
  if(!missing(new.populations)){
   NEWPOP <- TRUE
   populations <- vector("list",npops)
  }else{
   NEWPOP <- FALSE
  } 

## PROGRESS #########################
 progr <- progressBar()
#####################################



for(xx in 1:n.region.names){

 bial <- popGetBial(object,xx)

 if(length(bial)!=0){ # if a biallelic position exists
    
    ### neue Populationen
    if(NEWPOP){                                              # wenn neu Populationen definiert
       for(yy in 1:npops){
           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]
           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }   
       }
    }else{
     populations <- object@region.data@populations[[xx]]
    }
     
    # synonym = 0
    ### Statistic
    synonymous  <- object@region.data@synonymous[[xx]]
    bialcount   <- 1:length(synonymous) # Anzahl der biallelic sites
    
    res <- apply(pairs,2,function(pop){
           
           #segsites  <- get_segsites(bial,populations)
	   segsites   <- get_segsites(bial,list(populations[[pop[1]]],populations[[pop[2]]]))

           if( length(segsites[[1]])!=0 ){
           pop1mono  <- bialcount[-segsites[[1]]]
           }else{
           pop1mono  <- bialcount
           }

	   if( length(segsites[[2]])!=0 ){
           pop2mono  <- bialcount[-segsites[[2]]]
	   }else{
	   pop2mono  <- bialcount
           }    

           pop12mono <- intersect(pop1mono,pop2mono)
           
           Pnpop1   <- sum(synonymous[segsites[[1]]]==FALSE,na.rm=TRUE)
           Pspop1   <- sum(synonymous[segsites[[1]]]==TRUE,na.rm=TRUE)
           Pnpop2   <- sum(synonymous[segsites[[2]]]==FALSE,na.rm=TRUE)
           Pspop2   <- sum(synonymous[segsites[[2]]]==TRUE,na.rm=TRUE)
           
           pop1pop2        <- vector("list",1)
           pop1pop2[[1]]   <- unique(c(populations[[pop[1]]],populations[[pop[2]]]))
           
           segsites    <- get_segsites(bial,pop1pop2)
           # check fixed polymorphisms
           fixed       <- intersect(pop12mono,segsites[[1]])
           
           Dn          <- sum(synonymous[fixed]==FALSE,na.rm=TRUE)
           Ds          <- sum(synonymous[fixed]==TRUE,na.rm=TRUE)
     
           Pn    <- Pnpop1 + Pnpop2
           Ps    <- Pspop1 + Pspop2
        
           NI    <- (Pn/Ps)/(Dn/Ds)
           alpha <- 1-NI

           fisher.P <- NaN
           if(do.fisher.test){
           fisher.P <- fisher.test(rbind(c(Ps,Ds),c(Pn,Dn)))$p.value
           }           

     return(c(Pn,Ps,Dn,Ds,NI,alpha,fisher.P))
            
    })
    
    if(length(populations)>2){
      res           <- t(res)
      colnames(res) <- MKnames
      rownames(res) <- nn
      MKmulti[[xx]] <- res
    }else{
      MKTWO[xx,]    <- res
    }

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################
    
 }
}  
  
 if(length(populations)>2){
  MKmulti            <- as.matrix(MKmulti)
  rownames(MKmulti)  <- region.names
  colnames(MKmulti)  <- "McDonald & Kreitman Test"
  object@MKT         <- MKmulti
 }else{
  rownames(MKTWO)    <- region.names
  colnames(MKTWO)    <- c("P_nonsyn","P_syn","D_nonsyn","D_syn","neutrality.index","alpha","fisher.P.value")
  object@MKT         <- MKTWO
 }
 
 return(object)
 })

 # --------------------------------------------------------------------
 # popAll
 # --------------------------------------------------------------------
 
# setGeneric("popAll", function(object,new.populations="list") standardGeneric("popAll"))
# setMethod("popAll", "GENOME",
# function(object,new.populations){
 
#   if(!missing(new.populations)){
#   res    <- popNeutrality(object,new.populations)
#   res    <- popFSTN(res,new.populations)
#   res    <- popFSTH(res,new.populations)
#   }else{
#   res    <- popNeutrality(object)
#   res    <- popFSTN(res)
#   res    <- popFSTH(res)
#   }
  
#  return(res)
# })

# ---------------------------------------------------------------------
#
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Calc Neutrality (calc_freqstats)
# ----------------------------------------------------------------------
 setGeneric("detail.stats", function(object,new.populations=FALSE, new.outgroup=FALSE, subsites=FALSE,biallelic.structure=FALSE, mismatch.distribution=FALSE, site.spectrum=TRUE, site.FST=FALSE) standardGeneric("detail.stats"))
 setMethod("detail.stats", "GENOME",
 
 function(object,new.populations,new.outgroup,subsites,biallelic.structure,mismatch.distribution,site.spectrum, site.FST){
 
 region.names   <- object@region.names 
 n.region.names <- length(region.names)

 if(object@big.data){region.names <- NULL} # because of memory space

 
 
 object@Pop_Detail$sites       <- "ALL"
 object@Pop_Detail$calculated  <- TRUE
  
 # new populations

 if(missing(new.populations)){
 npops                         <- length(object@populations)
 object@Pop_Detail$Populations <- object@populations

 }else{
 npops                             <- length(new.populations)
 object@Pop_Detail$Populations     <- new.populations
 }

# Outgroup
 if(missing(new.outgroup)){
 object@Pop_Detail$Outgroup <- object@populations
 }else{
 object@Pop_Detail$Outgroup <- new.outgroup
 }
 
 
 # Outgroup
  if(!missing(new.outgroup)){
   NEWOUT <- TRUE
  }else{
   NEWOUT <- FALSE
  } 

 change      <- object@region.stats
 # Names
 # nam <- paste("pop",1:npops)
 #------------------------------

 
  if(!missing(new.populations)){
   NEWPOP <- TRUE
   populations <- vector("list",npops)
  }else{
   NEWPOP <- FALSE
  } 

# for change 
Xminor.allele.freqs     <- vector("list",n.region.names) 
Xbiallelic.structure    <- vector("list",n.region.names) 
Xsite.FST               <- vector("list",n.region.names)

# Init
 init        <- matrix(,n.region.names,npops)
 MDSD        <- init
 MDG1        <- init
 MDG2        <- init

 # Names
 nam <- paste("pop",1:npops)
 #------------------------------
 rownames(MDSD)     <- region.names
 colnames(MDSD)     <- nam
 rownames(MDG1)     <- region.names
 colnames(MDG1)     <- nam
 rownames(MDG2)     <- region.names
 colnames(MDG2)     <- nam


## PROGRESS #########################
 progr <- progressBar()
#####################################

for(xx in 1:n.region.names){

### if Subsites ----------------------------------

bial <- popGetBial(object,xx)

if(subsites[1]!=FALSE){

if(subsites=="transitions" & length(bial!=0)){
   tran       <- which(object@region.data@transitions[[xx]]==TRUE)
   bial       <- bial[,tran,drop=FALSE]
  # object@Pop_Detail$sites <- "transitions"   
}

if(subsites=="transversions" & length(bial!=0)){
   transv       <- which(object@region.data@transitions[[xx]]==FALSE)
   bial         <- bial[,transv,drop=FALSE]
  # object@Pop_Detail$sites <- "transversions"   
} 

if(subsites=="syn" & length(bial!=0)){
   syn        <- which(object@region.data@synonymous[[xx]]==TRUE)
   bial       <- bial[,syn,drop=FALSE]
   #object@Pop_Detail$sites <- "synonymous"
}
if(subsites=="nonsyn" & length(bial!=0)){
   nonsyn        <- which(object@region.data@synonymous[[xx]]==FALSE)
   bial          <- bial[,nonsyn,drop=FALSE]
   #object@Pop_Detail$sites <- "nonsynonymous"
}

if(subsites=="intron" & length(bial!=0)){
   intron        <- which(object@region.data@IntronSNPS[[xx]]==TRUE)
   bial          <- bial[,intron,drop=FALSE]
   #object@Pop_Detail$sites <- "introns"
}

if(subsites=="utr" & length(bial!=0)){
   utr           <- which(object@region.data@UTRSNPS[[xx]]==TRUE)
   bial          <- bial[,utr,drop=FALSE]
  # object@Pop_Detail$sites <- "utr"
}

if(subsites=="coding" & length(bial!=0)){
   #coding           <- which(!is.na(object@region.data@synonymous[[xx]])==TRUE)
   coding           <- which(object@region.data@CodingSNPS[[xx]]==TRUE)
   bial             <- bial[,coding,drop=FALSE]
  # object@Pop_Detail$sites <- "coding"
}


}# End if subsites

############### ---------------------------------


  if(length(bial)!=0){ # if a biallelic position exists  
    
    if(NEWPOP){ # wenn neu Populationen definiert
    
       for(yy in 1:npops){
           if(is.character(new.populations[[yy]])){
              #populations[[yy]] <- match(new.populations[[yy]],rownames(object@DATA[[xx]]@matrix_pol))
              populations[[yy]] <- match(new.populations[[yy]],rownames(bial))
              naids             <- which(!is.na(populations[[yy]]))
              populations[[yy]] <- populations[[yy]][naids]
           }else{
              populations[[yy]] <- new.populations[[yy]]
              ids               <- which(populations[[yy]]>dim(bial)[1])
              if(length(ids)>0){ populations[[yy]] <- populations[[yy]][-ids]}   
           }      
       }
           
       #----------------------#
       temp         <- delNULLpop(populations)
       populations  <- temp$Populations
       popmissing   <- temp$popmissing
       #----------------------#   
       if(length(populations)==0){next} # Keine Population vorhanden


    }else{
     populations <- object@region.data@populations[[xx]] # if there is no new population
    }
    
      if(NEWPOP)  {if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}} # nur die Populationen, die existieren
      if(!NEWPOP) {if(length(object@region.data@popmissing[[xx]])!=0){popmissing <- object@region.data@popmissing[[xx]];respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}}

    
    ## Outgroup
  if(NEWOUT){

	if(is.character(new.outgroup)){
           outgroup <- match(new.outgroup,rownames(bial)) 
           naids    <- which(!is.na(outgroup))
           outgroup <- outgroup[naids]  
        }else{
           outgroup <- new.outgroup
           ids      <- which(outgroup>dim(bial)[1])
           if(length(ids)>0){outgroup <- outgroup[-ids]}   
        }
        
        if(length(outgroup)==0){outgroup <- FALSE}

  }else{
    outgroup <- object@region.data@outgroup[[xx]]
  }
  # --------------------------------------

# Statistics

    # biallelic.structure
     if(biallelic.structure){
       data.list                            <- list(biallelic.sites=object@region.data@biallelic.sites[[xx]])
       sxx                                  <- calc_sxsfss(bial,populations=populations,data=data.list)    
       Xbiallelic.structure[[xx]]           <- sxx
     }

     # mismatch distribution
     if(mismatch.distribution){
      if(!object@Pop_Slide$calculated){
      data.list              <- list(n.nucleotides=object@region.data@n.nucleotides[[xx]],n.valid.sites=object@n.valid.sites[xx],
                                transitions=object@region.data@transitions[[xx]],biallelic.compositions=object@region.data@biallelic.compositions[[xx]])
      }else{
      data.list              <- list(n.nucleotides=NULL,n.valid.sites=NULL,
                                transitions=object@region.data@transitions[[xx]],biallelic.compositions=NULL)
      }
      res                     <- fstcalc(bial,populations=populations,data=data.list)
      thetaT                  <- res$PIW  
      misdis                  <- mismatch(bial,populations=populations,thetaT) 
      MDSD[xx,respop]         <- misdis[1,]
      MDG1[xx,respop]         <- misdis[2,]
      MDG2[xx,respop]         <- misdis[3,]
     }
     # site.spectrum
    if(site.spectrum){ 
     res2                              <- jointfreqdist(bial,populations=populations,outgroup=outgroup)
     Xminor.allele.freqs[[xx]]         <- res2$jfd
    }

    if(site.FST){
     #res3 <- apply(bial,2,function(x){return(fstcalc(as.matrix(x),populations,outgroup,data=list())$FSTALL)})
     #Xsite.FST[[xx]] <- res3
     Xsite.FST[[xx]] <- site_FST(bial, populations)
    }

    # -----------------------
    #   fill detailed Slots
    # -----------------------
   

  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################
 
 }
  
}
   change@minor.allele.freqs  <- Xminor.allele.freqs
   change@biallelic.structure <- Xbiallelic.structure
   change@site.FST            <- Xsite.FST 
   object@MDSD <- MDSD
   object@MDG1 <- MDG1
   object@MDG2 <- MDG2
   object@region.stats  <- change
  
  return(object)
 })

setGeneric("get.detail", function(object, biallelic.structure=FALSE) standardGeneric("get.detail"))
 setMethod("get.detail", "GENOME",
 
 function(object,biallelic.structure){


if(!object@Pop_Detail$calculated){stop("Statistics have to be calculated first !")}
 

# biallelic.structure
if(biallelic.structure){

res   <- vector("list",length(object@region.names)) 
npops <- length(object@Pop_Detail$Populations)

 for(yy in 1:length(object@region.names)){
   
   bial.sites <- object@region.data@biallelic.sites[[yy]]
   if(length(bial.sites)==0){next}

   popmat     <- matrix(0,npops,length(bial.sites))
   colnames(popmat) <- bial.sites
   rownames(popmat) <- paste("pop",1:npops)

   for(xx in 1:npops){
    
    data <- object@region.stats@biallelic.structure[[yy]]$POP[,1][[xx]] # SX "I am polymorph rest is monomorph "
    ids  <- match(data,bial.sites)
    popmat[xx,ids] <- 1

    data <- object@region.stats@biallelic.structure[[yy]]$POP[,2][[xx]] # SXF "I am monomorph rest is polymorph "
    ids  <- match(data,bial.sites)
    popmat[xx,ids] <- 2

    data <- object@region.stats@biallelic.structure[[yy]]$POP[,3][[xx]] # SF  "I am mono rest is mono with same mono value"
    ids  <- match(data,bial.sites)
    popmat[xx,ids] <- 3
   
    data <- object@region.stats@biallelic.structure[[yy]]$POP[,4][[xx]] # SS "I am mono rest is mono with different mono value"
    ids  <- match(data,bial.sites)
    popmat[xx,ids] <- 4


   }

 res[[yy]] <- popmat

 }

return(res)

}# end of if biallelic.structure


  res   <- matrix(,length(object@region.names),3)
  pops1 <- vector("list",length(object@Pop_Detail$Populations))
  
  for(xx in 1:length(object@Pop_Detail$Populations)){
     res[,1]       <- object@MDSD[,xx]
     res[,2]       <- object@MDG1[,xx]
     res[,3]       <- object@MDG2[,xx]
     
     colnames(res) <- c("MDSD","MDG1","MDG2")
     rownames(res) <- object@region.names
     pops1[[xx]]   <- res 
  }
 
 pops1            <- as.matrix(pops1)
 rownames(pops1)  <- paste("pop",1:length(object@Pop_Detail$Populations))
 colnames(pops1)   <- "Mismatch Distribution"
return(pops1)

})

##################################################################################
## Sliding Window Transform 
##################################################################################


setGeneric("sliding.window.transform", function(object,width=7,jump=5,type=1,start.pos=FALSE,end.pos=FALSE,whole.data=TRUE) standardGeneric("sliding.window.transform"))
 setMethod("sliding.window.transform", "GENOME",

function(object,width,jump,type,start.pos,end.pos,whole.data){

n.region.names  <- length(object@region.names)

if(whole.data){

  if(n.region.names > 1 && length(object@BIG.BIAL)==0){ # check BIG.BIAL when you split the data before

    object         <- concatenate_to_whole_genome(object,n.region.names)
    n.region.names <- 1
  
  }


 SPLIT <- sliding.window.transform.fast(object,width,jump,type,start.pos,end.pos)

 return(SPLIT)


}else{
 
 if(object@snp.data){
  RETURN <- slide.snp.sep(object,width,jump,type) 
 }else{
  RETURN <- sliding.window.transform.new(object,width,jump,type,start.pos,end.pos)
 }
 return(RETURN)

}
  
  

 # Check if there is a new population defined ! 
 # Init 

## PROGRESS #########################
 progr <- progressBar()
#####################################

genomeobj               <-  new("GENOME") 
ddatt                   <-  new("region.data")
NAM                     <-  NULL
count                   <-  1

for(xx in 1:n.region.names){

 if(object@big.data){
  bial <- object@region.data@biallelic.matrix[[xx]] 
 }else{
  bial <- popGetBial(object,xx)
 }

 if(length(bial)==0){
     
     # genomeobj@DATA[[count]] <- NULL
     # genomeobj@GEN[[count]]  <- NULL
     ddatt@biallelic.matrix[[count]] <- as.matrix(NaN)
     NAM                             <- c(NAM,object@region.names[xx])
     count <- count + 1
     next
     
 }

 if(length(bial)!=0){ # if a biallelic position exists
     bial.sites            <- as.numeric(colnames(bial))                 
  # Calculate type 3 = reference Biallelic Matrix. type 5 = reference GEN   
  if(type==1){repeatlength <- ceiling( (dim(bial)[2]-width+1)/jump)}
  if(type==2){repeatlength <- ceiling( (object@n.sites[xx]-width+1)/jump)} 
      
  if(repeatlength<=0){

     #genomeobj@DATA[[count]] <- NULL
     #genomeobj@GEN[[count]]  <- NULL
     ddatt@biallelic.matrix[[count]]                     <- as.matrix(NaN)
     NAM                                                 <- c(NAM,object@region.names[xx])
     count <- count + 1

     next
  }
   
   for(zz in 1:repeatlength){
 
        
        start      <- ((zz-1) * jump + 1)
        end 	   <- ((zz-1) * jump + width) 
       

       if(type==1){
        
       
	 window	   <- start:end 
        #genomeobj@DATA[[count]] <- new("DATA")
        #ddatt@biallelic.matrix[[count]] <- object@region.data@biallelic.matrix[[xx]][,window,drop=FALSE]

         if(object@big.data){
               ddatt@biallelic.matrix[[count]]  <- ff(bial[,window,drop=FALSE],dim=dim(bial[,window,drop=FALSE]))               
         }else{
               ddatt@biallelic.matrix[[count]]  <- bial[,window,drop=FALSE]
         }

        # ddatt@biallelic.matrix[[count]] <- bial[,window,drop=FALSE]
        ddatt@outgroup[[count]]         <- object@region.data@outgroup[[xx]]
        ddatt@populations[[count]]      <- object@region.data@populations[[xx]]
        ddatt@popmissing[[count]]       <- object@region.data@popmissing[[xx]]
        ddatt@synonymous[[count]]       <- object@region.data@synonymous[[xx]][window]
        ddatt@transitions[[count]]      <- object@region.data@transitions[[xx]][window]
        ddatt@biallelic.sites[[count]]  <- object@region.data@biallelic.sites[[xx]][window]
        
        # genomeobj@GEN [[count]] <- new("GEN")
        count                   <- count + 1     
        NAM                     <- c(NAM,object@region.names[xx])

       }
        
       
        
        #   value                 <- calc_freqstats(bial[,window,drop=FALSE],populations)
        #   TajD[zz,respop]       <- value$taj_D
        #   S   [zz,respop]       <- value$THETA["S",]
        #}

        if(type==2){
      
            # bialpos        <- as.numeric(colnames(bial))
            ids            <- (bial.sites >= start) & (bial.sites<=end) 
            # ids            <- is.element(bialpos,window)
            bialpos        <- which(ids)

           if(length(bialpos)!=0){
               
              #genomeobj@DATA[[count]] <- new("DATA")
              #genomeobj@DATA[[count]]@biallelic.matrix <- object@DATA[[xx]]@biallelic.matrix[,bialpos,drop=FALSE]
              #genomeobj@DATA[[count]]@outgroup         <- object@DATA[[xx]]@outgroup
              #genomeobj@GEN [[count]] <- new("GEN")

              if(object@big.data){
               ddatt@biallelic.matrix[[count]]  <- ff(bial[,bialpos,drop=FALSE],dim=dim(bial[,bialpos,drop=FALSE]))               
              }else{
               ddatt@biallelic.matrix[[count]]  <- bial[,bialpos,drop=FALSE]
              }

              ddatt@outgroup[[count]]          <- object@region.data@outgroup[[xx]]
              ddatt@populations[[count]]       <- object@region.data@populations[[xx]]
              ddatt@popmissing[[count]]        <- object@region.data@popmissing[[xx]]
              ddatt@synonymous[[count]]        <- object@region.data@synonymous[[xx]] [bialpos]
              ddatt@transitions[[count]]       <- object@region.data@transitions[[xx]][bialpos]
              ddatt@biallelic.sites[[count]]   <- object@region.data@biallelic.sites[[xx]][bialpos]
              count                            <- count + 1     
              NAM                              <- c(NAM,object@region.names[xx])

         #     value                 <- calc_freqstats(bial[,bialpos,drop=F],populations)
         #     TajD[zz,respop]       <- value$taj_D
         #     S   [zz,respop]       <- value$THETA["S",]

           }else{

             #genomeobj@DATA[[count]] <- NULL
             #genomeobj@GEN [[count]] <- NULL
             ddatt@biallelic.matrix[[count]] <- as.matrix(NaN)
             NAM                             <- c(NAM,object@region.names[xx])
             count                           <- count + 1

           }
       } 
 
# slidenames[zz]  <- paste(start,"-",end,sep="")
 

} # end for Sliding

 
#colnames(TajD) <- paste("pop",1:npops)
#rownames(TajD) <- slidenames

#colnames(TajD) <- paste("pop",1:npops)
#rownames(TajD) <- slidenames

#TAJSLIDE[[xx]] <- TajD
#SSLIDE[[xx]]   <- S


  # PROGRESS #######################################################
    progr <- progressBar(xx,n.region.names, progr)
  ###################################################################

 }
}# End over all genes

#TAJSLIDE           <- as.matrix(TAJSLIDE)
#SSLIDE             <- as.matrix(SSLIDE) 
#rownames(TAJSLIDE) <- object@region.names
#colnames(TAJSLIDE) <- c("Tajimas.D")
#rownames(SSLIDE)   <- object@region.names
#colnames(SSLIDE)   <- c("n.segregating.sites")
#object@SLIDE       <- cbind(TAJSLIDE,SSLIDE)

genomeobj@populations               <- object@populations
genomeobj@region.names              <- paste(1:length(NAM),":",NAM,sep="")
genomeobj@genelength                <- length(NAM)

 if(length(ddatt@popmissing)==0){ddatt@popmissing <- vector("list",length(genomeobj@region.names))}
        # weil liste fuellen mit NULL nicht funktioniert, eh besser die Sachen schon vorher zu definieren

genomeobj@region.data               <- ddatt
genomeobj@Pop_Neutrality$calculated <- FALSE
genomeobj@Pop_FSTN$calculated       <- FALSE
genomeobj@Pop_FSTH$calculated       <- FALSE
genomeobj@Pop_MK$calculated         <- FALSE
genomeobj@Pop_Recomb$calculated     <- FALSE
genomeobj@Pop_Linkage$calculated    <- FALSE
genomeobj@Pop_Slide$calculated      <- TRUE
genomeobj@Pop_Detail$calculated     <- FALSE
genomeobj@big.data                  <- object@big.data

return(genomeobj)

 })
 
 






