ALLELES <- NLOCI <- NALLELES <- ALLELEFREQS <- adults <- loci_number <- offspring <- adults <- MOMS <- DADS <- OFFSPRING <- NULL
if (getRversion() >= '2.15.1') globalVariables("solomon")
solomon.env <- new.env()

solomon <- function() {

fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "Welcome to SOLOMON"                                             # Name the window
heading <- tklabel(tt, text="SOLOMON: Parentage Analysis Program",font=fontHeading) # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
fontSUB <- tkfont.create(family="times",size=15,weight="bold")
heading <- tklabel(tt, text="Data Set Creation and Power Analysis",font=fontSUB)      # add a heading
tkgrid(heading, columnspan=5)
#Allele_Frequency_Module###############################################################################################################################################################################################
PressedFREQ <- function()
{
#Create and name toplevel window ===============================================#
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Parentage Analysis"                                    # Name the window
heading <- tklabel(tt, text="Allele Frequencies",font=fontHeading)              #add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set working directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function() {
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
  }
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
tkgrid(tklabel(tt,text="     "))
#Create allele frequency file===================================================#
fontHeading <- tkfont.create(family="times",size=15,weight="bold")
heading <- tklabel(tt, text="Create Allele Frequency File with a Data Set",font=fontHeading)   #add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))
#Load Genotype File=============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("ALLELES", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Genotype File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Create button to run allele calc===============================================#
PressedOK <- function()      {
gdata<-ALLELES
population=gdata[,-1]                                                           #Here we are accessing only the genotypes
L=ncol(population)                                                              #how many columns are there?
locus_positions=(2*(unique(round((1:(L-2))/2)))+1)                              #find the starting column number for each locus
#lnamedata=population[1,locus_positions]                                        #create a dummy dataset, in case user names each column (instead of locus)differently
#lnames=colnames(lnamedata)
lnames=colnames(population)                                                     #locus names, from the header
OUT=NULL                                                                        #create a null dataset to append allele freqs to
for (x in locus_positions) {                                                    #begin for loop, to calculate frequencies for each locus
  alleles=c(population[,x],population[,x+1])                                    #For example, combine columns 1 and 2 for locus 1 (two columns because they are diploid)
  alleles2=as.data.frame(table(alleles))                                        #count each allele at locus x
  missing=alleles2[which(alleles2[,1]==0),2]                                    #count missing data at locus x, entered as '0' in this dataset (not used further for simplicity)
  if (length(which(alleles2[,1]==0))>0) {
  alleles2=alleles2[-which(alleles2[,1]==0),]}                                  #remove missing data (otherwise 0 would be counted in total number of alleles)
  alleles4=cbind(alleles2,alleles2[,2]/sum(alleles2[,2]))                       #calculate frequencies
  output=cbind(x,lnames[x],alleles4)                                            #combine x, locusname, and frequencies
  OUT <- rbind(OUT,output)
}
colnames(OUT) <- c("Number","Locus","allele","count","frequency")               #add column headers
Allelefreqs=OUT[,-1]
write.table(Allelefreqs,file="AlleleFrequencies_UserFriendly.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
Allelefreqs2=Allelefreqs[,-c(2,3)]
Allelefreqs2=cbind(as.numeric(Allelefreqs2[,1]),Allelefreqs2[,2])
write.table(Allelefreqs2,file="AlleleFrequencies_SOLOMON.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
}
label.bayes <- tklabel(tt, text="Press 'Run' to calculate allele frequencies:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label.bayes,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.bayes, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
#Create allele frequency file===================================================#
fontHeading <- tkfont.create(family="times",size=15,weight="bold")
heading <- tklabel(tt, text="Create an Allele Frequency File",font=fontHeading)   #add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set Number of Loci=============================================================#
label_sim5 <- tklabel(tt, text="Number of wanted loci:",font=fontTextLabel)
Name.sim5 <- tclVar()
entry.label5 <- tkentry(tt, width="10",textvariable=Name.sim5)                  #create entry fields
OnOK25 <- function()  {
	NameVal <- tclvalue(Name.sim5)
	msg <- paste("You have now set the number of loci to",NameVal)
	tkmessageBox(message=msg)
	assign("NLOCI", NameVal, envir = solomon.env)
}
OK.but25 <-tkbutton(tt,text="   OK   ",command=OnOK25)
tkbind(entry.label5, "<Return>",OnOK25)
tkgrid(tklabel(tt,text="     "), label_sim5, tklabel(tt,text="     "), entry.label5, tklabel(tt,text="     "), OK.but25, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim5, sticky="w")
tkgrid.configure(entry.label5, sticky="w")
#Set Number of alleles per locus================================================#
label_sim.number6 <- tklabel(tt, text="Number of alleles per locus:",font=fontTextLabel)
Name.sim6 <- tclVar()
entry.label.sim6 <- tkentry(tt, width="10",textvariable=Name.sim6)              #create entry fields
OnOK26 <- function()  {
	NameVal <- tclvalue(Name.sim6)
	msg <- paste("You have now set the number of alleles per locus to",NameVal)
	tkmessageBox(message=msg)
	assign("NALLELES", NameVal, envir = solomon.env)
}
OK.but26 <-tkbutton(tt,text="   OK   ",command=OnOK26)
tkbind(entry.label.sim6, "<Return>",OnOK26)
tkgrid(tklabel(tt,text="     "), label_sim.number6, tklabel(tt,text="     "), entry.label.sim6, tklabel(tt,text="     "), OK.but26, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number6, sticky="w")
tkgrid.configure(entry.label.sim6, sticky="w")
PressedOK7 <- function()      {
Nloci=as.numeric(NLOCI)
wantednumberalleles=as.numeric(NALLELES)
a=c(1:wantednumberalleles)
b=rev(a)
c=b/a
d=sum(c)
freqs=c/d
lowestallele=140
alleles2=cbind(seq(lowestallele,lowestallele+wantednumberalleles-1,1),freqs)
afreqs = rep(alleles2[,2],Nloci)
afreqs=cbind(sort(rep(1:Nloci,wantednumberalleles)),afreqs)
afreqs=cbind(afreqs[,1],afreqs)
write.table(afreqs[,-1],file="AlleleFrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
}
label.bayes2 <- tklabel(tt, text="Press 'Run' to create allele frequencies:",font=fontTextLabel)
run.button2 <- tkbutton(tt, text = "Run", command = PressedOK7)
tkgrid(tklabel(tt,text="     "),label.bayes2,tklabel(tt,text="     "),run.button2,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.bayes2, sticky="w")
tkgrid.configure(run.button2, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
label.FREQ <- tklabel(tt, text="Create allele frequency file:",font=fontTextLabel)
OK.but.FREQ <- tkbutton(tt, text = "FREQUENCIES", command = PressedFREQ)
tkgrid(tklabel(tt,text="     "),label.FREQ,tklabel(tt,text="     "),OK.but.FREQ,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.FREQ, sticky="w")
tkgrid.configure(OK.but.FREQ, sticky="w")
tkfocus(tt)
#Module2_test_data_sets#############################################################################################################################################################################################
#Module2_test_data_sets_11.26.12
PressedSIMS <- function()
{
#Create and name toplevel window ==============================================#
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Parentage Analysis"                                    # Name the window
heading <- tklabel(tt, text="SOLOMON: Create Simulated Data Sets",font=fontHeading)         # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set working directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function() {
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
  }
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Set Number of Parents==========================================================#
label_sim.number <- tklabel(tt, text="Number of wanted parents:",font=fontTextLabel)
Name.sim <- tclVar()
entry.label.sim <- tkentry(tt, width="10",textvariable=Name.sim)                #create entry fields
OnOK2 <- function() {
	NameVal <- tclvalue(Name.sim)
	msg <- paste("You have now set the number of wanted parents to",NameVal)
	tkmessageBox(message=msg)
	assign("Nparents", NameVal, envir = solomon.env)
}
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label1, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_sim.number, tklabel(tt,text="     "), entry.label.sim, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number, sticky="w")
tkgrid.configure(entry.label.sim, sticky="w")
#Set Number of offspring per pair===============================================#
label_sim.number2 <- tklabel(tt, text="Number of offspring per parent:",font=fontTextLabel)
Name.sim2 <- tclVar()
entry.label.sim2 <- tkentry(tt, width="10",textvariable=Name.sim2)              #create entry fields
OnOK22 <- function()  {
	NameVal <- tclvalue(Name.sim2)
	msg <- paste("You have now set the number of offspring per parent to",NameVal)
	tkmessageBox(message=msg)
	assign("Noffs_perpair", NameVal, envir = solomon.env)
}
OK.but22 <-tkbutton(tt,text="   OK   ",command=OnOK22)
tkbind(entry.label.sim2, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_sim.number2, tklabel(tt,text="     "), entry.label.sim2, tklabel(tt,text="     "), OK.but22, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number2, sticky="w")
tkgrid.configure(entry.label.sim2, sticky="w")
#Set genotyping error rate======================================================#
label_sim.number3 <- tklabel(tt, text="Genotyping error rate (eg 0.01):",font=fontTextLabel)
Name.sim3 <- tclVar()
entry.label.sim3 <- tkentry(tt, width="10",textvariable=Name.sim3)              #create entry fields
OnOK23 <- function()  {
	NameVal <- tclvalue(Name.sim3)
	msg <- paste("You have now set the genotyping error rate to",NameVal)
	tkmessageBox(message=msg)
	assign("error", NameVal, envir = solomon.env)
  }
OK.but23 <-tkbutton(tt,text="   OK   ",command=OnOK23)
tkbind(entry.label.sim3, "<Return>",OnOK23)
tkgrid(tklabel(tt,text="     "), label_sim.number3, tklabel(tt,text="     "), entry.label.sim3, tklabel(tt,text="     "), OK.but23, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number3, sticky="w")
tkgrid.configure(entry.label.sim3, sticky="w")
#Set Number of unrelated individuals============================================#
label_sim.number4 <- tklabel(tt, text="Number of unrelated individuals:",font=fontTextLabel)
Name.sim4 <- tclVar()
entry.label.sim4 <- tkentry(tt, width="10",textvariable=Name.sim4)              #create entry fields
OnOK24 <- function() {
	NameVal <- tclvalue(Name.sim4)
	msg <- paste("You have now set the number of unrelated individuals to",NameVal)
	tkmessageBox(message=msg)
	assign("Nunrelated", NameVal, envir = solomon.env)
}
OK.but24 <-tkbutton(tt,text="   OK   ",command=OnOK24)
tkbind(entry.label.sim4, "<Return>",OnOK24)
tkgrid(tklabel(tt,text="     "), label_sim.number4, tklabel(tt,text="     "), entry.label.sim4, tklabel(tt,text="     "), OK.but24, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number4, sticky="w")
tkgrid.configure(entry.label.sim4, sticky="w")
#Set Number of siblings============================================#
label_sib.number4 <- tklabel(tt, text="Number of full-siblings:",font=fontTextLabel)
Name.sib4 <- tclVar()
entry.label.sib4 <- tkentry(tt, width="10",textvariable=Name.sib4)              #create entry fields
OnOK24 <- function() {
	NameVal <- tclvalue(Name.sib4)
	msg <- paste("You have now set the number of full-siblings to",NameVal)
	tkmessageBox(message=msg)
	assign("Nsibs", NameVal, envir = solomon.env)
}
OK.but24 <-tkbutton(tt,text="   OK   ",command=OnOK24)
tkbind(entry.label.sib4, "<Return>",OnOK24)
tkgrid(tklabel(tt,text="     "), label_sib.number4, tklabel(tt,text="     "), entry.label.sib4, tklabel(tt,text="     "), OK.but24, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sib.number4, sticky="w")
tkgrid.configure(entry.label.sib4, sticky="w")
#Load Allele frequency==========================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=TRUE, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("ALLELEFREQS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Allele Frequency File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing allele frequencies:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Create button to run parentage script==========================================#
PressedOK <- function()      {
Nparents=as.numeric(Nparents)
Noffs_perpair=as.numeric(Noffs_perpair)
error=as.numeric(error)
Nunrelated=as.numeric(Nunrelated)
Nadults=Nparents*2                                                              #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
afreqs <- ALLELEFREQS
OUT=NULL
sims=function(sims)
{
alleles2=afreqs[which(afreqs[,1]==z),]
alleles3=cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   #table allele frequencies
homos=(alleles3[,2])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
hets=t(combn(alleles3[,2],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles3[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
n=1000000                                                                       #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
gfreqs1=rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            #create genotypes(by coloumn, 1 for each allele)
gfreqs2=rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
gtypes=cbind(gfreqs1,gfreqs2)
gtypes=gtypes[sample(1:length(gtypes[,1]),replace=F),]
sg1=gtypes[sample(1:length(gtypes[,1]),Nadults),]
OUT<<-cbind(OUT,sg1)
}
z=length(unique(afreqs[,1]))
C1=for(z in 1:z) {lapply(z,sims)}

parents=OUT
c=c(1:(ncol(OUT)))
odd=2*(unique(round(((c-2))/2)))+1
l=length(odd) * 1000
codes=seq(from=1,to=l,by=1000)
cols=sort(rep(codes,2))-1
Anumbs=matrix(cols,Nadults,ncol(OUT),byrow=T)
parents=as.numeric(parents)+Anumbs
#create full sib families (go down the list in pair)============================#
OUT2=NULL
sims=function(sims)  {
p1=parents[z,]
p2=parents[z+1,]
als=rep(1:2,length(p1)/2)
Noffs=Noffs_perpair                                                             #number of offspring per pair
OUT2=NULL
for (b in 1:Noffs){
pos1=sample(als,length(p1)/2,replace=TRUE)                                      #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
pos2=sample(als,length(p1)/2,replace=TRUE)
pos11=pos1+(seq(0,(length(p1)-1),2))
pos22=pos2+(seq(0,(length(p2)-1),2))
o1=p1[pos11]
o2=p2[pos22]
o3=sort(c(o1,o2))
o3=t(c(z,o3))
write.table(o3,file="SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
}
z=length(parents[,1])
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,sims)
Dads=parents[seq(from=1,to=length(parents[,1]),by=2),]
Moms=parents[seq(from=2,to=length(parents[,1]),by=2),]
Anumbs=matrix(cols,Nadults,ncol(OUT),byrow=T)                                   #see code before functions for adding 1000s (here am removing 1000s)
parents=as.numeric(parents)-Anumbs
Dads=parents[seq(from=1,to=length(parents[,1]),by=2),]
Moms=parents[seq(from=2,to=length(parents[,1]),by=2),]
Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Offs = as.matrix(Offs2[,-1])
Anumbs=matrix(cols,length(Offs[,1]),ncol(OUT),byrow=T)
Offs=as.numeric(Offs)-Anumbs
if (Noffs_perpair>1){
Offnames=ceiling(Offs2[,1]/2)                                                   #naming of offpsring
Offnames2=paste(Offnames,".",1:length(Offnames))
Offs=cbind(paste("Offspring",Offnames2),Offs)} else {
Offnames=ceiling(Offs2[,1]/2)
Offs=cbind(paste("Offspring",Offnames),Offs)}
Moms=cbind(paste("Mom",1:length(Moms[,1])),Moms)
Dads=cbind(paste("Dad",1:length(Dads[,1])),Dads)
#add error======================================================================#
Dadsg=Dads[,-1]
ldad=length(Dadsg)*error
pdad1=sample(1:length(Dadsg),ldad,replace=FALSE)
pdad2=Dadsg[sample(1:length(Dadsg),ldad,replace=FALSE)]
Dadsg2=replace(Dadsg,pdad1,pdad2)
Dads=cbind(Dads[,1],Dadsg2)

Offsg=Offs[,-1]
loff=length(Offsg)*error
poff1=sample(1:length(Offsg),loff,replace=FALSE)
poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
Offsg2=replace(Offsg,poff1,poff2)
Offs=cbind(Offs[,1],Offsg2)

#Momsg=Moms[,-1]
#ldad=length(Momsg)*error
#pdad1=sample(1:length(Momsg),ldad,replace=FALSE)
#pdad2=Momsg[sample(1:length(Momsg),ldad,replace=FALSE)]
#Momsg2=replace(Momsg,pdad1,pdad2)
#Moms=cbind(Moms[,1],Momsg2)
#===============================================================================#
colsnmz=rep("Locus",ncol(Offs)-1)
colsnmz2=sort(rep(1:(length(colsnmz)/2),2))
colsnmz3=paste(colsnmz,colsnmz2)
colsnmz3=c("IDs",colsnmz3)
colnames(Offs)<- colsnmz3
colnames(Dads)<- colsnmz3
colnames(Moms)<- colsnmz3
#Add unrealated individuals=====================================================#
if (Nunrelated>0){
Nadults=Nunrelated*3                                                            #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
afreqs <- ALLELEFREQS
OUT=NULL
sims=function(sims) {
alleles2=afreqs[which(afreqs[,1]==z),]
alleles3=cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   #table allele frequencies
homos=(alleles3[,2])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
hets=t(combn(alleles3[,2],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles3[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
n=1000000                                                                       #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
gfreqs1=rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            #create genotypes(by coloumn, 1 for each allele)
gfreqs2=rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
gtypes=cbind(gfreqs1,gfreqs2)
gtypes=gtypes[sample(1:length(gtypes[,1]),replace=F),]
sg1=gtypes[sample(1:length(gtypes[,1]),Nadults),]
OUT<<-cbind(OUT,sg1)
}
z=length(unique(afreqs[,1]))
C1=for(z in 1:z) {lapply(z,sims)}

parents=OUT
c=c(1:(ncol(OUT)))
odd=2*(unique(round(((c-2))/2)))+1
l=length(odd) * 1000
codes=seq(from=1,to=l,by=1000)
cols=sort(rep(codes,2))-1
Anumbs=matrix(cols,Nadults,ncol(OUT),byrow=T)
parents=as.numeric(parents)+Anumbs
parents=as.numeric(parents)-Anumbs
unrelated=cbind(paste("Individual",1:length(parents[,1])),parents)
Lid=length(unrelated[,1])/3
m1=unrelated[1:Lid,]
d1=unrelated[(Lid+1):(Lid*2),]
j1=unrelated[((Lid*2)+1):(Lid*3),]
colsnmz=rep("Locus",ncol(Offs)-1)
colsnmz2=sort(rep(1:(length(colsnmz)/2),2))
colsnmz3=paste(colsnmz,colsnmz2)
colsnmz3=c("IDs",colsnmz3)
colnames(Offs)<- colsnmz3
colnames(Dads)<- colsnmz3
colnames(Moms)<- colsnmz3
Moms=rbind(Moms,m1)
Dads=rbind(Dads,d1)
Offs=rbind(Offs,j1)
}
write.table(Dads,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
write.table(Moms,file="Moms_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
write.table(Offs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
unlink("SimOffs.txt")
#Begin add siblings=============================================================#
#This scripts creates a desired number of siblings and splits them between adults and offspring files
#setwd("C:/POPS/")
#ALLELEFREQS<- read.table("Steelhead_freqs.txt", header=TRUE, sep="\t", na.strings="?",dec=".", strip.white=TRUE)
Nsibs=as.numeric(Nsibs)
#error=0.01
if (Nsibs > 0) {
Noffs_per_pair=2                                                                #could change this, but as stands only evaluating 1 pair of siblings per parent-pair
Nadults=Nsibs*2                                                                 #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
afreqs <- ALLELEFREQS
OUT=NULL
sims=function(sims)
{
alleles2=afreqs[which(afreqs[,1]==z),]
alleles3=cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   #table allele frequencies
homos=(alleles3[,2])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)
hets=t(combn(alleles3[,2],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles3[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
n=1000000                                                                       #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
gfreqs1=rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            #create genotypes(by coloumn, 1 for each allele)
gfreqs2=rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
gtypes=cbind(gfreqs1,gfreqs2)
gtypes=gtypes[sample(1:length(gtypes[,1]),replace=F),]
sg1=gtypes[sample(1:length(gtypes[,1]),Nadults),]
OUT<<-cbind(OUT,sg1)
}
z=length(unique(afreqs[,1]))
C1=for(z in 1:z) {lapply(z,sims)}
parents=OUT
c=c(1:(ncol(OUT)))
odd=2*(unique(round(((c-2))/2)))+1
l=length(odd) * 1000
codes=seq(from=1,to=l,by=1000)
cols=sort(rep(codes,2))-1
Anumbs=matrix(cols,Nadults,ncol(OUT),byrow=T)
parents=as.numeric(parents)+Anumbs
#create full sib families (go down the list in pair)============================#
OUT2=NULL
sims=function(sims)  {
N=1:length(parents[,1])
u=sample(N,1)
u2=sample(N,1)
p1=parents[u,]
p2=parents[u2,]
als=rep(1:2,length(p1)/2)
Noffs=Noffs_per_pair                                                            #number of offspring per pair
OUT2=NULL
for (b in 1:Noffs){
pos1=sample(als,length(p1)/2,replace=TRUE)                                      #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
pos2=sample(als,length(p1)/2,replace=TRUE)
pos11=pos1+(seq(0,(length(p1)-1),2))
pos22=pos2+(seq(0,(length(p2)-1),2))
o1=p1[pos11]
o2=p2[pos22]
o3=sort(c(o1,o2))
o3=t(c(z,o3))
write.table(o3,file="SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
}
z=length(parents[,1])
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,sims)                 #used to move down list, now sample randomly so is just used to produce wanted number of offspring
Dads=parents[seq(from=1,to=length(parents[,1]),by=2),]
Moms=parents[seq(from=2,to=length(parents[,1]),by=2),]
Anumbs=matrix(cols,Nadults,ncol(OUT),byrow=T)                                   #see code before functions for adding 1000s (here am removing 1000s)
parents=as.numeric(parents)-Anumbs
Dads=parents[seq(from=1,to=length(parents[,1]),by=2),]
Moms=parents[seq(from=2,to=length(parents[,1]),by=2),]
Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Offs = as.matrix(Offs2[,-1])
Anumbs=matrix(cols,length(Offs[,1]),ncol(OUT),byrow=T)
Offs=as.numeric(Offs)-Anumbs
if (Noffs_per_pair>1){
Offnames=ceiling(Offs2[,1]/2)                                                   #naming of offpsring
Offnames2=paste(Offnames,".",1:length(Offnames))
Offs=cbind(paste("Sibling",Offnames2),Offs)} else {
Offnames=ceiling(Offs2[,1]/2)
Offs=cbind(paste("Sibling",Offnames),Offs)}
Moms=cbind(paste("Mom",1:length(Moms[,1])),Moms)
Dads=cbind(paste("Dad",1:length(Dads[,1])),Dads)
#add error======================================================================#
Offsg=Offs[,-1]
loff=length(Offsg)*error
poff1=sample(1:length(Offsg),loff,replace=FALSE)
poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
Offsg2=replace(Offsg,poff1,poff2)
Offs=cbind(Offs[,1],Offsg2)
#===============================================================================#
colsnmz=rep("Locus",ncol(Offs)-1)
colsnmz2=sort(rep(1:(length(colsnmz)/2),2))
colsnmz3=paste(colsnmz,colsnmz2)
colsnmz3=c("IDs",colsnmz3)
colnames(Offs)<- colsnmz3
#calculate shared alleles among pairs of siblings===============================#
sib1=Offs[seq(from=1,to=length(Offs[,1]),by=2),]
sib2=Offs[seq(from=2,to=length(Offs[,1]),by=2),]
write.table(sib1,file="Sib1.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
write.table(sib2,file="Sib2.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
dadsibs <- read.table("Dads_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
sib1 <- read.table("Sib1.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
dadsibs=rbind(dadsibs,sib1)
write.table(dadsibs,file="Dads_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
juvsibs <- read.table("Juveniles_sim.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
sib2 <- read.table("Sib2.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
juvsibs=rbind(juvsibs,sib2)
write.table(juvsibs,file="Juveniles_sim.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
unlink("Sib1.txt")
unlink("Sib2.txt")
unlink("SimOffs.txt")
}
}
label.bayes <- tklabel(tt, text="Press 'Run' to create test data sets:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label.bayes,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label.bayes, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
label.SIMS <- tklabel(tt, text="Create test data sets:",font=fontTextLabel)
OK.but.SIMS <- tkbutton(tt, text = "SIMS", command = PressedSIMS)
tkgrid(tklabel(tt,text="     "),label.SIMS,tklabel(tt,text="     "),OK.but.SIMS,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.SIMS, sticky="w")
tkgrid.configure(OK.but.SIMS, sticky="w")
tkfocus(tt)
#Modue3_power_analysis########################################################################################################################################################################################################################
PressedOK_Power <- function() {
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Power Analysis "                                   # Name the window
heading <- tklabel(tt, text="Power Analysis",font=fontHeading)              # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set genotyping error rate======================================================#
label_mismatches <- tklabel(tt, text="Choose minimum number of loci:",font=fontTextLabel)
Name1 <- tclVar("")
entry.label2 <- tkentry(tt, width="4",textvariable=Name1)                       #create entry fields
OnOK2 <- function()  {
	NameVal <- tclvalue(Name1)
 	msg <- paste("You have now set the number of loci to:",NameVal)
	tkmessageBox(message=msg)
	assign("loci_number", NameVal, envir = solomon.env)
}
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label2, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_mismatches, tklabel(tt,text="     "), entry.label2, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_mismatches, sticky="w")
tkgrid.configure(entry.label2,sticky="w")
#Set working Directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()   {
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Load Adults File===============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("adults", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Adults File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing adult genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Offspring File============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("offspring", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
#Create button to run parentage script==========================================#
PressedOK <- function()     {
Adults1=adults
a=ncol(Adults1)
lociset=seq(from=2,to=a,by=2)
loci=as.numeric(loci_number)
s1=sample(lociset,loci,replace=FALSE)
s2=sort(c(s1,s1+1))
Adults=Adults1[,s2]
Offspring1<- offspring
Offspring=Offspring1[,s2]
#Begin calcualtion of Pr(Z) and Pr(delta)=======================================#
massive=function(massive)
{
AAT40=c(Adults[,a],Adults[,(a+1)])
locus=as.data.frame(table(AAT40))
Frequency=locus[,2]/sum(locus[,2])
Data=cbind(locus, Frequency)
if (Data[1,1]==0) {Data=Data[-1,]} else Data=Data
n1=((length(AAT40))/2)
A1=(Data[,3]*(2*n1))
A2=(Data[,3]^2)*n1
AA=A1-A2
AAA=cbind(Data,AA)
AAT402=c(Offspring[,a],Offspring[,(a+1)])
locus2=as.data.frame(table(AAT402))
Frequency2=locus2[,2]/sum(locus2[,2])
Data2=cbind(locus2, Frequency2)
if (Data2[1,1]==0) {Data2=Data2[-1,]} else Data2=Data2
n2=((length(AAT402))/2)
B1=(Data2[,3]*(2*n2))
B2=(Data2[,3]^2)*n2
BB=B1-B2
BBB=cbind(Data2,BB)
AAA1=match(BBB[,1],AAA[,1])
g=which(AAA1>0)
g1=AAA1[g]
AAA1=AAA[g1,]
BBB1=match(AAA[,1],BBB[,1])
j=which(BBB1>0)
j1=BBB1[j]
BBB1=BBB[j1,]
AB=AAA1[,4]*BBB1[,4]
AB=sum(AB)
y=AAA1[,1]
Agenotypes=expand.grid(y,y)
remove=which(Agenotypes[,1]==Agenotypes[,2])
Agenotypes=Agenotypes[-remove,]
z=BBB1[,1]
Bgenotypes=expand.grid(z,z)
remove=which(Bgenotypes[,1]==Bgenotypes[,2])
Bgenotypes=Bgenotypes[-remove,]
one=match(Agenotypes[,1],AAA1[,1])
two=match(Agenotypes[,2],AAA1[,1])
one=AAA1[one,3]
two=AAA1[two,3]
Agfreq=one*two*n1*2
oneb=match(Bgenotypes[,1],BBB1[,1])
twob=match(Bgenotypes[,2],BBB1[,1])
oneb=BBB1[oneb,3]
twob=BBB1[twob,3]
Ogfreq=oneb*twob*n2*2
Gfreqs=Agfreq*Ogfreq
Gfreqs=floor(Gfreqs)
Gfreq=sum(Gfreqs)/2
PrB=(AB-Gfreq)/(n1*n2)
write.table(PrB,file="Output_Per_Locus_Exclusion_Probabilities.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
a=ncol(Adults)
C1=for(a in (2*(unique(round((1:(a-2))/2)))+1)) lapply(a,massive)
PrBs <- read.table("Output_Per_Locus_Exclusion_Probabilities.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
nn1=length(Adults[,1])
nn2=length(Offspring[,1])
Pr_delta=prod(PrBs[,1])
Expected_Number_of_False_Pairs=Pr_delta*(nn1*nn2)
pvalue1=cbind(loci_number,Expected_Number_of_False_Pairs,Pr_delta)
write.table(pvalue1,file="Output_Expected_Number_of_False_Pairs.txt",row.names=FALSE,col.names=T,sep="\t",append=T)
unlink("Output_Per_Locus_Exclusion_Probabilities.txt")
}
label_exclusion <- tklabel(tt, text="Press 'Run' to calculate expected number of false pairs:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label_exclusion,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label_exclusion, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
label.POWER <- tklabel(tt, text="Power analysis:",font=fontTextLabel)
OK.but.POWER <- tkbutton(tt, text = "POWER", command = PressedOK_Power)
tkgrid(tklabel(tt,text="     "),label.POWER,tklabel(tt,text="     "),OK.but.POWER,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.POWER, sticky="w")
tkgrid.configure(OK.but.POWER, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkfocus(tt)
#Module4 _Exclusion_Noparents#############################################################################################################################################################################################################################
fontSUB <- tkfont.create(family="times",size=15,weight="bold")
heading <- tklabel(tt, text="No Known Parents",font=fontSUB)                   # add a heading
tkgrid(heading, columnspan=5)
PressedOK <- function() {
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Parentage Analysis "                                   # Name the window
heading <- tklabel(tt, text="Exclusion",font=fontHeading)                       # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set number of mismatching loci=================================================#
label_mismatches <- tklabel(tt, text="Set number of loci to let mismatch:",font=fontTextLabel)
Name1 <- tclVar("0")
entry.label2 <- tkentry(tt, width="4",textvariable=Name1)                       #create entry fields
OnOK2 <- function()  {
	NameVal <- tclvalue(Name1)
 	msg <- paste("You have now set the number of mismatches to",NameVal)
	tkmessageBox(message=msg)
	assign("mismatch", NameVal, envir = solomon.env)
}
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label2, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_mismatches, tklabel(tt,text="     "), entry.label2, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_mismatches, sticky="w")
tkgrid.configure(entry.label2,sticky="w")
#Set working Directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()  {
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Load Adults File===============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings=-900, dec=".", strip.white=TRUE)
  assign("adults", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Adults File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing adult genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Offspring File============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings=-900, dec=".", strip.white=TRUE)
  assign("offspring", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
#Create button to run parentage script==========================================#
PressedOK <- function()     {
mismatch=as.numeric(mismatch)                                                   #used because it is entered as a character through GUI
Adults1=adults
Offspring1=offspring
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Offspring=Offspring1[,c(2:a)]
Anames=Adults1[,1]
Onames=Offspring1[,1]
Adults=Adults1[,-1]
Offspring=Offspring1[,-1]
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
A=1:Aindivids
O=1:Oindivids
G=expand.grid(A,O)
AG=G[,1]
AO=G[,2]
Ads=Adults[AG,]
Offs=Offspring[AO,]
IdnamesA=Anames[AG]
IdnamesO=Onames[AO]
write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Names=cbind(IdnamesA,IdnamesO)
matches=function(matches)
{
A=Ads[,z]-Offs[,z]
B=Ads[,(z+1)]-Offs[,(z+1)]
C=Ads[,z]-Offs[,(z+1)]
D=Ads[,(z+1)]-Offs[,z]
f=A*B*C*D
f=(f^2)*10
ss=which(is.na(f)==TRUE)
f=replace(f,ss,0)
identify=which(f>0)
f=replace(f,identify,1)
f=cbind(z,f)
write.table(f,file="Sort.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=T)
}
z=ncol(Ads)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
a=length(U[,1])
stuff=rowSums(U)
Sorted=cbind(Names,stuff)
matches=which(Sorted[,3]<(mismatch+1))
Actual=sort(stuff)
IDS=which(stuff<(mismatch+1))
PAdults=Ads[IDS,]
POffspring=Offs[IDS,]
nput=length(matches)
Putativepairs=Sorted[matches,]
PAdults=cbind(Putativepairs[,c(1,3)],PAdults)
POffspring=cbind(Putativepairs[,c(2,3)],POffspring)
names(PAdults)[1]<-"ID"
names(POffspring)[1]<-"ID"
sorts=function(sorts) {
tell=rbind(PAdults[f,],POffspring[f,])
write.table(tell,file="Output_genotypes_tmp.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}
f=length(PAdults[,1])
C1= for(f in 1:f) lapply(f,sorts)
unlink("Sort.txt")
unlink("IdnamesA.txt")
unlink("IdnamesO.txt")
#Outputfile processing==========================================================#
putative<- read.table("Output_genotypes_tmp.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
write.table(putative,file="Output_genotypes.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=FALSE)
unlink("Output_genotypes_tmp.txt")                                              #strange 3 lines allow genotypes to be overwritten with multiple runs
parents=putative[seq(from=1,to=length(putative[,1]),by=2),c(1,2)]
offspring=putative[seq(from=2,to=length(putative[,1]),by=2),c(1,2)]
out1=cbind(parents,offspring)
out2=out1[order(out1[,1]),-2]
out2=cbind(out2,"")
out2[,1] <- as.character(out2[,1])
out2[,4] <- as.character(out2[,4])
ids=2:length(out2[,1])
for (x in 2:length(out2[,1])) {                                                 #using 2 here (because 1 minus nothing is invalid : add in first parent later
    if (out2[x,1]!=out2[x-1,1]) {out2[x,4] <- out2[x,1]}
}
out2[1,4] <- out2[1,1]
out3=out2[,c(4,2,3)]
out4=cbind(out3,"")
out4[,4] <- as.character(out4[,4])
count=data.frame(table(out2[,1]))
for (x in 1:length(count[,1])){
out4[match(count[x,1],out3[,1]),4]=count[x,2]
}
noffs=sort(unique(out4[,4]))
if (noffs[1]=="") {noffs=noffs[-1]}
noffs=sort(as.numeric(noffs),decreasing=TRUE)
noffs=as.character(noffs)                                                       #comment out this line if want highest reproductive success at bottom of file
OUT=NULL
for (x in 1:length(noffs)){
  sim1=which(out4[,4]==noffs[x])
      VALS=NULL
      for (y in 1:length(sim1)){
      vals=seq(sim1[y],sim1[y]+(as.numeric(noffs[x])-1),1)
      vals=cbind(noffs[x],vals)
      VALS=rbind(VALS,vals)
      }
  sim2=out4[as.numeric(VALS[,2]),]
  OUT=rbind(OUT,sim2)
}
parent.centric=OUT[,c(1,4,2,3)]
colnames(parent.centric)<-c("Parent","Number_Offspring","Offspring","Number_Loci_mismatching")
write.table(parent.centric,file="Output_by_parent.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
#Begin Ouputfiles_offspringcentric==============================================#
parents=putative[seq(from=1,to=length(putative[,1]),by=2),c(1,2)]
offspring=putative[seq(from=2,to=length(putative[,1]),by=2),c(1,2)]
out1=cbind(parents,offspring)
out2=out1[,-2]
out2=cbind(out2,"")
out2[,2] <- as.character(out2[,2])
out2[,4] <- as.character(out2[,4])
out2=out2[order(out2[,2]),]
ids=2:length(out2[,2])
for (x in 2:length(out2[,2])) {                                                 #using 2 here (because 1 minus nothing is invalid : add in first parent later
    if (out2[x,2]!=out2[x-1,2]) {out2[x,4] <- out2[x,2]}
}
out2[1,4] <- out2[1,2]
out3=out2[,c(4,1,3)]
out4=cbind(out3,"")
out4[,4] <- as.character(out4[,4])
count=data.frame(table(out2[,2]))
count=count[which(count[,2]>=1),]
for (x in 1:length(count[,1])){
out4[match(count[x,1],out3[,1]),4]=count[x,2]
}
noffs=sort(unique(out4[,4]))
if (noffs[1]=="") {noffs=noffs[-1]}
noffs=sort(as.numeric(noffs),decreasing=TRUE)
noffs=as.character(noffs)
OUT=NULL
for (x in 1:length(noffs)){
  sim1=which(out4[,4]==noffs[x])
      VALS=NULL
      for (y in 1:length(sim1)){
      vals=seq(sim1[y],sim1[y]+(as.numeric(noffs[x])-1),1)
      vals=cbind(noffs[x],vals)
      VALS=rbind(VALS,vals)
      }
  sim2=out4[as.numeric(VALS[,2]),]
  OUT=rbind(OUT,sim2)
}
parent.centric=OUT[,c(1,4,2,3)]
colnames(parent.centric)<-c("Offspring","Number_Parents","Parents","Number_Loci_mismatching")
write.table(parent.centric,file="Output_by_offspring.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
}
label_exclusion <- tklabel(tt, text="Press 'Run' to perform exlcusion:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label_exclusion,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label_exclusion, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
fontTextLabel <- tkfont.create(family="times",size=14)
label.exclusion1 <- tklabel(tt, text="Perform exlcusion with no known parents:",font=fontTextLabel)
OK.but <- tkbutton(tt, text = "EXCLUSION", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label.exclusion1,tklabel(tt,text="     "),OK.but,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label.exclusion1, sticky="w")
tkgrid.configure(OK.but, sticky="w")
tkfocus(tt)
#Module5_bayes_noparents#####################################################################################################################################################################################################################
#Bayesian parentage_noParents
PressedOK2 <- function()     {
#Create and name toplevel window================================================#
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Bayesian Parentage Analysis "                           # Name the window
heading <- tklabel(tt, text="SOLOMON: Bayesian Parentage Analysis",font=fontHeading)# add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set working Directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()  {
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Set Number of Sims=============================================================#
label_sim.number <- tklabel(tt, text="Number of simulated data sets:",font=fontTextLabel)
Name.sim <- tclVar(1000)
entry.label.sim <- tkentry(tt, width="10",textvariable=Name.sim)                #create entry fields
OnOK2 <- function()  {
	NameVal <- tclvalue(Name.sim)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("wanted_reps", NameVal, envir = solomon.env)
  }
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label1, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_sim.number, tklabel(tt,text="     "), entry.label.sim, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special1=tklabel(tt, text="Recommended: Microsatellites=1000 SNPs=100",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special1, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number, sticky="w")
tkgrid.configure(entry.label.sim,label_Special1, sticky="w")
#Set Number of "Geonotypes"=============================================================#
label_Ntotal <- tklabel(tt, text="Number of simulated genotypes:",font=fontTextLabel)
Name.Ntotal <- tclVar(50000000)
entry.label.Ntotal <- tkentry(tt, width="10",textvariable=Name.Ntotal)                #create entry fields
OnOKN <- function()  {
	NameVal <- tclvalue(Name.Ntotal)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("Ntotal", NameVal, envir = solomon.env)
  }
OK.butN <-tkbutton(tt,text="   OK   ",command=OnOKN)
tkbind(entry.label1, "<Return>",OnOKN)
tkgrid(tklabel(tt,text="     "), label_Ntotal, tklabel(tt,text="     "), entry.label.Ntotal, tklabel(tt,text="     "), OK.butN, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special=tklabel(tt, text="Recommended: Microsatellites=50,000,000 SNPs=500,000",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_Ntotal, sticky="w")
tkgrid.configure(entry.label.Ntotal,label_Special, sticky="w")
#Load Adults File===============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("adults", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Adults File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing adult genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Offspring File============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("offspring", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
#Create button to run parentage script==========================================#
PressedOK <- function()     {
  Adults1=adults
  Offspring1=offspring
  wanted_reps=as.numeric(wanted_reps)
  loci=ncol(Adults1)
  Adults=Adults1[,c(2:loci)]                                                    #assumes that there is one column of id names ; could modify this as needed.
  Offspring=Offspring1[,c(2:loci)]
  total <- ncol(Adults)/2                                                       #For Progress bar
  if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = total, width = 300)}
afreqs=function(afreqs)                                                         #Begin Master simulation function
{
locus_name=L
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, ceiling(locus_name/2), title=paste("Locus", ceiling(locus_name/2),"of ",total))}
vect=c(Adults[,L],Adults[,L+1])                                                 #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
alleles=data.frame(table(vect))
alleles=alleles[order(alleles[,1]),]
if (as.numeric(as.character(alleles[1,1]))==0) {alleles=alleles[-1,]}           #removes missing data
if(length(alleles[,1])==1) {                                                    #deals with monomorphic loci by adding 1 very strange allele (799)
alleles=cbind(vect[1],alleles[2])
alleles[2,]<-c(799,1)}
alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
homos=(alleles2[,3])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
hets=t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
csum=cumsum(as.numeric(gfreqs[,3]))
gfreqs1=cbind(gfreqs,csum)
Nadults=length(Adults[,1])
Noffs=length(Offspring[,1])
#===============================================================================#end locus-specific HWE genotype frequency calculations
alength=length(alleles2[,1])
for (Y in 1:wanted_reps) {
positions=1:length(gfreqs[,1])
sg3=sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                      #sample the repeated genotype positions, by the number of adults
sadults=gfreqs[sg3,1:2]                                                         #index gfreqs to create genotypes
og3=sample(positions,Noffs,replace=TRUE,prob=gfreqs[,3])                        #create juvenile genotyes
soffs=gfreqs[og3,1:2]
soffs=cbind(as.numeric(soffs[,1]),as.numeric(soffs[,2]))
asp=cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
asp=rbind(cbind(locus_name,0,0),asp)
for (i in 1:Nadults) {
parent1=as.numeric(sadults[i,1])                                                #first allele in parent
parent2=as.numeric(sadults[i,2])                                                #second allele in parent
p1=soffs-parent1
p2=soffs-parent2
pp1=which(p1[,1]==0)
pp2=which(p1[,2]==0)
allele1=unique(c(pp1,pp2))
p21=which(p2[,1]==0)
p22=which(p2[,2]==0)
allele2=unique(c(p21,p22))
Out51=cbind(parent1,length(allele1))
Out52=cbind(parent2,length(allele2))
Out53=cbind(0,Noffs-length(unique(c(allele1,allele2))))
Out5=rbind(Out51,Out52,Out53)
if(parent2==parent1) {Out5=Out5[-1,]}                                           #remove 1 of alleles count if homozygous
if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
  diffs=sum(Out5[,2])-Noffs                                                     #comment out to be more conservative!
  maxa=max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
  pos=which(Out5[,2]==maxa)
  Out5[pos,2]<-Out5[pos,2]-diffs}
m1=match(Out5[,1],asp[,2])
m2=asp[m1,3]+as.numeric(Out5[,2])
asp[m1,3]<-m2
asp<-asp
}
write.table(asp,file="out.sims",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}
}
L=ncol(Adults)
C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
if(.Platform$OS.type=="windows"){close(pb)}
#Bayes P-value==================================================================#
OUT<- read.table("out.sims", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
locname=unique(OUT[,1])                                                         #compile calculations for each locus
OUTALL=NULL
for (z in locname) {
  Loc1=OUT[which(OUT[,1]==z),]
  allfreqs=unique(Loc1[,2])
  OUT2=NULL
    for (x in allfreqs) {
    a1<-Loc1[which(Loc1[,2]==x),]
    a2=sum(a1[,3])
    a3=cbind(x,a2)
    OUT2 <- rbind(OUT2, a3)
  }
  OUT3=cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
  OUTALL <- rbind(OUTALL, cbind(z,OUT3))
}
#Create multilocus genotypes, 50000 at a time, calculate number of loci that mismatch, and calculate freqs of shared alleles
NL=length(unique(OUTALL[,1]))
ngtypes=1                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
Ntotal=as.numeric(Ntotal)
inreps=50000     #was tested as 10000 for SNPS
repnumber=round(Ntotal/inreps)                                                #this is the rep number to get 100,000 values.  can adjust accordingly
asp=cbind(0:NL,rep(0,length(0:NL)))
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = Ntotal , width = 300)}
for (n in 1:repnumber) {
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, n*inreps, title=paste("Genotype", n*inreps,"of ",Ntotal))}
    OUT=NULL
    for (r in unique(OUTALL[,1])) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*Loco[,4])
        findo=which(alleles3[,2]==0)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
  distm=apply(OUT, 1, function(x)sum(x == 1))
  distm2=data.frame(table(distm))
  m1=match(distm2[,1],asp[,1])
  m2=asp[m1,2]+distm2[,2]
  asp[m1,2]<-m2
  asp<-asp
}
if(.Platform$OS.type=="windows"){close(pb)}
#tabulate number of multilocus genotypes with x mismatches
d2=asp
d3=cbind(d2,d2[,2]/sum(d2[,2]))
#===============================================================================# Create plot of exclusionary power.  Also, some necessary data formatting
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "Posterior Processing", min = 0, max = 7 , width = 300)}
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 1, title=paste("Performing Exclusion"))}
Adults<- adults
Offs<- offspring
Nads=length(Adults[,1])
Noffs=length(Offs[,1])
asp=d3
asp=cbind(asp,NL-as.numeric(as.character(asp[,1])))  #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
asp=cbind(asp,cumsum(asp[,2]))
asp=cbind(asp,asp[,5]/max(asp[,5]))
#find minimum Nloci to mismatch (could modify this) and perform exclusion with decided-upon mismatches==#
distm=cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
#mismatch=min(which(round(distm[,6],1)==.9))                                    #deprecated
mismatch=which.min(abs(distm[,6] - .5))                                         #could cause issues with a value of .5 chosen here - could be too low.  Change to higher if all loci analyzed have phi < 1.
Adults1=Adults                                                                   #begin exclusion
Offspring1=Offs
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Offspring=Offspring1[,c(2:a)]
Anames=Adults1[,1]
Onames=Offspring1[,1]
Adults=Adults1[,-1]
Offspring=Offspring1[,-1]
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
A=1:Aindivids
O=1:Oindivids
G=expand.grid(A,O)
AG=G[,1]
AO=G[,2]
Ads=Adults[AG,]
Offs=Offspring[AO,]
IdnamesA=Anames[AG]
IdnamesO=Onames[AO]
write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Names=cbind(IdnamesA,IdnamesO)
matches=function(matches)
{
A=Ads[,z]-Offs[,z]
B=Ads[,(z+1)]-Offs[,(z+1)]
C=Ads[,z]-Offs[,(z+1)]
D=Ads[,(z+1)]-Offs[,z]
f=A*B*C*D
f=(f^2)*10
ss=which(is.na(f)==TRUE)
f=replace(f,ss,0)
identify=which(f>0)
f=replace(f,identify,1)
f=cbind(z,f)
write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Ads)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
a=length(U[,1])
stuff=rowSums(U)
Sorted=cbind(Names,stuff)
matches=which(Sorted[,3]<(mismatch+1))
Actual=sort(stuff)
IDS=which(stuff<(mismatch+1))
PAdults=Ads[IDS,]
POffspring=Offs[IDS,]
nput=length(matches)
Putativepairs=Sorted[matches,]
PAdults=cbind(Putativepairs[,c(1,3)],PAdults)
POffspring=cbind(Putativepairs[,c(2,3)],POffspring)
names(PAdults)[1]<-"ID"
names(POffspring)[1]<-"ID"
sorts=function(sorts)
{
tell=rbind(PAdults[f,],POffspring[f,])
write.table(tell,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
f=length(PAdults[,1])
C1= for(f in 1:f) lapply(f,sorts)
unlink("Sort.txt")
unlink("IdnamesA.txt")
unlink("IdnamesO.txt")
#Calculate phi for each number mismatching loci=================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 2, title=paste("Calculating phi")) }
Putative<- read.table("Output_genotypes.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
observed=data.frame(table(Putative[,2]))
observed=cbind(observed,observed[,2]/2)                                         #done becuase each number is written twice in file (once for parent and once for off)
zerom=0:mismatch                                                                #this chunk adds 0s for mismatches where there were no observed putative pairs
zerom2=which(is.na(match(zerom,observed[,1])))
if (length(zerom2>0)) {observed=observed[,3]
for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}  #not really 0, to prevent divide by 0
}   else {observed=observed[,3]}
expected=distm[1:(mismatch+1),7]                                                #using cumulatinve sum   (more conservative)
#expected=distm[1:(mismatch+1),3]*Nads*Noffs                                     #not using cumulative sum
phi=expected/observed
phi=replace(phi,which(phi>=1),1)
Offs<- offspring
actualTrue=length(grep("Off",Offs[,1]))
info=cbind(actualTrue,expected,observed,phi)
#calculate phi and index values ================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 3, title=paste("Creating distriubtion of falsely-shared alleles")) }
phibase=phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
observed=observed[min(which(phi[]==1))-1]                                       #do the same for observed
testob=which(observed==0.000000001)
phi2=cbind(1:length(phi),phi)
if (length(testob>0)) {phi4=phi2[-testob,]} else {phi4=phi2}
nmismatch=min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
index=phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
index=index[which(index>-1)]
if (length(index)>1) {
  if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
phi=phi[index]
index=index-1
#Create Plot ===================================================================#
pdf(file="Output_Dataset_Power.pdf")
x=0:(length(info[,1])-1)
y1=info[,3]
y2=info[,2]
p1=which(y2==0)
y2=replace(y2,p1,.000000001)
p1=which(y1<y2)
y3=replace(y1,p1,y2[p1])
y2=y2+1
y3=y3+1
par(mar=c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
ats=c(0,1,2,3,4,5,6,7,8,9)
labs=c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
axis(side=2,ats,labs)
lines(x,log10(y3),lwd=2)
lines(x,log10(y2),lwd=2)
points(x,log10(y3),pch=21,bg="green",cex=2)
points(x,log10(y2),pch=21,bg="blue",cex=2)
legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
yphi=y2/y3
par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
points(x,yphi,cex=2,pch=21,bg="gray")
mtext("Number of Mismatching Loci",side=1,line=1.94)
dev.off()
info2=cbind(x,info[,4])
colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
#===============================================================================#
ngtypes=20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
inreps=10000
repnumber=round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
#writes values at all loci, to be analyzed further below
for (n in 1:repnumber) {
    OUT=NULL
    for (r in unique(OUTALL[,1])) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*Loco[,4])
        findo=which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        alleles3=alleles3[-which(alleles3[,2]==0),]
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
for (i in index) {                                                              #loop over numbers of mismatched loci
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                    #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distp=cbind(i,DIST)
    write.table(distp,file="False_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
 }
}

#Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
OUT9=NULL
for (n in index){
  Putative2=Putative[which(Putative[,2]==n),]
  write.table(Putative2, file="Putative.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    #if (length(distp[,1])>1000) {                                                #need at least 1000 values , else assigned 0 (may be redundant now)
       OBS=NULL
       afreqs=function(afreqs) {
        PutL=Putative2[,c(L+2,L+3)]
        PutLadults=PutL[seq(from=1,to=length(PutL[,1]),by=2),]                  #Combine putative pairs alleles into a single row
        PutLoffs=PutL[seq(from=2,to=length(PutL[,1]),by=2),]
        Puts=cbind(PutLadults,PutLoffs)
        c1=Puts[,1]-Puts[,3]                                                    #find the matching alleles
        c2=Puts[,1]-Puts[,4]
        c3=Puts[,2]-Puts[,3]
        c4=Puts[,2]-Puts[,4]
        Puts2=cbind(Puts,c1,c2,c3,c4)
        P5=replace(Puts2[,5],which(Puts2[,5]!=0),-10)
        P6=replace(Puts2[,6],which(Puts2[,6]!=0),-10)
        P7=replace(Puts2[,7],which(Puts2[,7]!=0),-10)
        P8=replace(Puts2[,8],which(Puts2[,8]!=0),-10)
        P5=replace(P5,which(P5==0),1)
        P6=replace(P6,which(P6==0),1)
        P7=replace(P7,which(P7==0),1)
        P8=replace(P8,which(P8==0),1)
        Puts3=cbind(Puts,P5,P6,P7,P8)
        Puts4=cbind((Puts3[,1]*Puts3[,5]),(Puts3[,1]*Puts3[,6]),(Puts3[,2]*Puts3[,7]),(Puts3[,2]*Puts3[,8]))
        alleles2=OUTALL[which(OUTALL[,1]==L),]
        alfreq1=alleles2[match(Puts4[,1],alleles2[,2]),4]
        alfreq2=alleles2[match(Puts4[,2],alleles2[,2]),4]                       #find the actual allele values
        alfreq3=alleles2[match(Puts4[,3],alleles2[,2]),4]
        alfreq4=alleles2[match(Puts4[,4],alleles2[,2]),4]
        Puts5=cbind(alfreq1,alfreq2,alfreq3,alfreq4)                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
        R1=replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
        R2=replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
        R3=replace(Puts5[,3],which(is.na(Puts5[,3])==TRUE),1)
        R4=replace(Puts5[,4],which(is.na(Puts5[,4])==TRUE),1)
        Puts6=cbind(R1,R2,R3,R4)
        Put_share=apply(Puts6, 1, min)                                          #find row minimum
        Put_share2=apply(Puts4, 1, max)                                         #find shared allele name
        Put_share3=c(Put_share,Put_share2)
        OBS <<- cbind(OBS,Put_share3)
    }
    L=ncol(Putative2)-2
    C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    lengths=length(OBS[,1])/2
    if (lengths==1) {OBA3=t(OBS[(lengths+1):(2*lengths),])} else {OBA3=OBS[(lengths+1):(2*lengths),]}
    write.table(OBA3,file="True_shared_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)  #Actual shared alleles   #This file wouuld be useful as output for people to identify shared and mismatching loci
    if (lengths==1) {OBS=t(OBS[1:lengths,])} else {OBS=OBS[1:lengths,]}   #formatting for if there is only a single pair
    obsp=apply(OBS, 1, prod)
    #}  else obsp=rep(0,(length(Putative2[,1]))/2)
OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
}
#calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic)==#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 4, title=paste("Calculating Posterior Component 1"))}
OBA3<- read.table("True_shared_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
for (n in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
      OUT=NULL
      for (r in unique(OUTALL[,1])) {                                           #This first section calculates products of parental alleles (True distribution)
        vect=c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
        alleles=data.frame(table(vect))
        alleles=alleles[order(alleles[,1]),]
        if (as.numeric(as.character(alleles[1,1]))<=0) {alleles=alleles[-1,]}
        alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))
        gtrue=sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
        OUT <- cbind(OUT,gtrue)

        if(n==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele
                            mm=alleles2[match(OBA3[i,ceiling(r/2)],alleles2[,1]),3]
                            write.table(cbind(r,mm),file="Shared_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
                            }
                }
        }
    for (i in index) {
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                   #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distt<-cbind(i,DIST)
    write.table(distt,file="True_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
}
#Calculate lamdaphi=============================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 5, title=paste("Calculating Posterior Component 2"))}
Putative3<- read.table("Putative.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Putadults=Putative3[seq(from=1,to=length(Putative3[,1]),by=2),1]                #Combine putative pairs alleles into a single row
Putoffs=Putative3[seq(from=2,to=length(Putative3[,1]),by=2),1]
Names=cbind(as.character(Putadults),as.character(Putoffs))
empirical=cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
distp <- read.table("False_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
P2=NULL
for (i in index) {                                                              #loop over numbers of mismatched loci
  empirical2=empirical[which(as.numeric(as.character(empirical[,3]))==i),]
  if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
  if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
    a3=distp[which(distp[,1]==i),2]
    DIST<-as.data.frame(a3)
    P=NULL
    for (b in 1:length(empirical2[,1])) {
      p1=length(which(DIST[,1] <= as.numeric(empirical2[b,4]) ))
      if (p1==0) {p1=0.00001}
      p2=cbind(empirical2[b,1],empirical2[b,2],p1)
      p3=cbind(p2,as.numeric(p2[,3])/length(DIST[,1]))
      P <- rbind(P,p3)
      }
    }
  P2<-rbind(P2,cbind(i,P))
}
lamdaphi=as.numeric(as.character(P2[,5]))
#Calculate lamda|phic ==========================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 6, title=paste("Calculating Posterior Component 3"))}
lamdaphic_dist<- read.table("True_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Observed<- read.table("Shared_allele_freqs.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
mm1=which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
Observed[mm1,2] <- 1                                                            #replace NAs with 1
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
lamdaphic=apply(U, 1, prod)
l1=length(which(OUT9[,2]==0))
if (l1>0) lamdaphic=c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)    #double check values by hand
P3=cbind(P2,lamdaphic)
for (i in index) {                                                              #loop over numbers of mismatched loci
   e2=P3[which(as.numeric(as.character(P3[,1]))==i),]
   if(length(e2)==6) {e2=t(e2)}                                                 #deals with one indiviudal formatting
   if (e2[1,5]==0) {write.table(e2[,6], file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{                    #deals with not enough reps
    a3=lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
    DIST<-as.data.frame(a3)
  for (b in 1:length(e2[,1])) {                                                 #calculate p values
    p1=length(which(DIST[,1] <= e2[b,6]))
    p2=p1/length(DIST[,1])
    write.table(p2, file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    }
  }
}
lamdaphic<- read.table("lamdaphic.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
#Put it all together with Bayes theorem!========================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 7, title=paste("Calculating Posterior"))}
vals=cbind(P2,lamdaphic[,1])
philength=cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
phis=rep(philength[,2],philength[,3])
vals=cbind(vals,phis)
colnames(vals)<-c("Nmismatch","Parent","Off","ignore","lamdaphi","lamdaphic","phi")
phi=as.numeric(as.character(vals[,7]))
lamdaphi=as.numeric(as.character(vals[,5]))
lamdaphic=as.numeric(as.character(vals[,6]))
lamdaphi=replace(lamdaphi,which(lamdaphi==0),1)
lamdaphic=replace(lamdaphic,which(lamdaphic==0),1)
pval=(lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
pval=cbind(vals[,2],vals[,3],vals[,1],pval)
pval=pval[order(as.numeric(pval[,4])),]
colnames(pval) <- c("Adult","Offspring","NL_mismatch","Probability of pair being false given frequencies of shared alleles")
#write.table(vals, file="Posterior_Components_Bayes.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
if(.Platform$OS.type=="windows"){close(pb)}
unlink("False_allele_freqs.txt")                                                #clear all sims files
unlink("True_allele_freqs.txt")
unlink("Shared_allele_freqs.txt")
unlink("lamdaphic.txt")
unlink("Putative.txt")
unlink("True_shared_freqs.txt")
unlink("Output_genotypes.txt")
unlink("*.sims")
rm(list=ls())
}
label.bayes <- tklabel(tt, text="Press 'Run' to perform Bayesian parentage analysis:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label.bayes,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label.bayes, sticky="w")
tkgrid.configure(run.button, sticky="w")
#blank end spaces and configure lables#=========================================#
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
fontTextLabel <- tkfont.create(family="times",size=14)
label.exclusion1 <- tklabel(tt, text="Perform Bayesian parentage analysis with no known parents:",font=fontTextLabel)
OK.but <- tkbutton(tt, text = "BAYES", command = PressedOK2)
tkgrid(tklabel(tt,text="     "),label.exclusion1,tklabel(tt,text="     "),OK.but,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label.exclusion1, sticky="w")
tkgrid.configure(OK.but, sticky="w")
#Module6_Bayes_NoParents_withSiblings###########################################################################################################################################################################################################
#Bayesian parentage_noParents_withSiblings
PressedOKSIB <- function()     {
#Create and name toplevel window================================================#
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Bayesian Parentage Analysis"                           # Name the window
heading <- tklabel(tt, text="Bayesian Parentage Analysis with Siblings",font=fontHeading)# add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set working Directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()  {
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Set Number of Sims=============================================================#
label_sim.number <- tklabel(tt, text="Number of simulated data sets:",font=fontTextLabel)
Name.sim <- tclVar(1000)
entry.label.sim <- tkentry(tt, width="10",textvariable=Name.sim)                #create entry fields
OnOK2 <- function()  {
	NameVal <- tclvalue(Name.sim)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("wanted_reps", NameVal, envir = solomon.env)
  }
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label1, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_sim.number, tklabel(tt,text="     "), entry.label.sim, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))

SpecialFont <- tkfont.create(family="times",size=11)
label_Special1=tklabel(tt, text="Recommended: Microsatellites=1000 SNPs=100",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special1, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number, sticky="w")
tkgrid.configure(entry.label.sim,label_Special1, sticky="w")
#Set Number of "Geonotypes"=============================================================#
label_Ntotal <- tklabel(tt, text="Number of simulated genotypes:",font=fontTextLabel)
Name.Ntotal <- tclVar(50000000)
entry.label.Ntotal <- tkentry(tt, width="10",textvariable=Name.Ntotal)                #create entry fields
OnOKN <- function()  {
	NameVal <- tclvalue(Name.Ntotal)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("Ntotal", NameVal, envir = solomon.env)
  }
OK.butN <-tkbutton(tt,text="   OK   ",command=OnOKN)
tkbind(entry.label1, "<Return>",OnOKN)
tkgrid(tklabel(tt,text="     "), label_Ntotal, tklabel(tt,text="     "), entry.label.Ntotal, tklabel(tt,text="     "), OK.butN, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special=tklabel(tt, text="Recommended: Microsatellites=50,000,000 SNPs=500,000",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_Ntotal, sticky="w")
tkgrid.configure(entry.label.Ntotal,label_Special, sticky="w")
#Set genotyping error rate======================================================#
label_sim.number3 <- tklabel(tt, text="Genotyping error rate:",font=fontTextLabel)
Name.sim3 <- tclVar(0.01)
entry.label.sim3 <- tkentry(tt, width="10",textvariable=Name.sim3)              #create entry fields
OnOK23 <- function()  {
	NameVal <- tclvalue(Name.sim3)
	msg <- paste("You have now set the genotyping error rate to",NameVal)
	tkmessageBox(message=msg)
	assign("error", NameVal, envir = solomon.env)
  }
OK.but23 <-tkbutton(tt,text="   OK   ",command=OnOK23)
tkbind(entry.label.sim3, "<Return>",OnOK23)
tkgrid(tklabel(tt,text="     "), label_sim.number3, tklabel(tt,text="     "), entry.label.sim3, tklabel(tt,text="     "), OK.but23, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special=tklabel(tt, text="Typical rates: Microsatellites=0.01 SNPs=0.001",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_Ntotal, sticky="w")
tkgrid.configure(entry.label.Ntotal,label_Special, sticky="w")
tkgrid.configure(label_sim.number3, sticky="w")
tkgrid.configure(entry.label.sim3, sticky="w")
#Load Adults File===============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("adults", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Adults File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing adult genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Offspring File============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("offspring", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
##add in genotyping error and vestigal null values from ModuleX_sibships========#
Nparents=10000                                                                  #number of sibilings and POPs to generate (2 siblings per pair are also defaults)
Noffs_perpair=2
Nunrelated=0
#Create button to run parentage script==========================================#
PressedOK <- function()      {
Nparents=as.numeric(Nparents)
Noffs_perpair=as.numeric(Noffs_perpair)
error=as.numeric(error)
Nunrelated=as.numeric(Nunrelated)
Nadults=Nparents*2                                                              #here to be used as the number of breeders (2* the total number of pairs and number of offspring)
#creage ALLELEFREQS file
gdata<-adults
population=gdata[,-1]                                                           #Here we are accessing only the genotypes
L=ncol(population)                                                              #how many columns are there?
locus_positions=(2*(unique(round((1:(L-2))/2)))+1)                              #find the starting column number for each locus
#lnamedata=population[1,locus_positions]                                        #create a dummy dataset, in case user names each column (instead of locus)differently
#lnames=colnames(lnamedata)
lnames=colnames(population)                                                     #locus names, from the header
OUT=NULL                                                                        #create a null dataset to append allele freqs to
for (x in locus_positions) {                                                    #begin for loop, to calculate frequencies for each locus
  alleles=c(population[,x],population[,x+1])                                    #For example, combine columns 1 and 2 for locus 1 (two columns because they are diploid)
  alleles2=as.data.frame(table(alleles))                                        #count each allele at locus x
  missing=alleles2[which(alleles2[,1]==0),2]                                    #count missing data at locus x, entered as '0' in this dataset (not used further for simplicity)
  if (length(which(alleles2[,1]==0))>0) {
  alleles2=alleles2[-which(alleles2[,1]==0),]}                                  #remove missing data (otherwise 0 would be counted in total number of alleles)
  alleles4=cbind(alleles2,alleles2[,2]/sum(alleles2[,2]))                       #calculate frequencies
  output=cbind(x,lnames[x],alleles4)                                            #combine x, locusname, and frequencies
  OUT <- rbind(OUT,output)
}
colnames(OUT) <- c("Number","Locus","allele","count","frequency")               #add column headers
Allelefreqs=OUT[,-1]
#write.table(Allelefreqs,file="AlleleFrequencies_UserFriendly.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
Allelefreqs2=Allelefreqs[,-c(2,3)]
Allelefreqs2=cbind(as.numeric(Allelefreqs2[,1]),Allelefreqs2[,2])
write.table(Allelefreqs2,file="AlleleFrequencies_SOLOMON.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
#}
afreqs<- read.table("AlleleFrequencies_SOLOMON.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
#Begin creation of siblings and calculations of shared alleles between full sibs and Parent-Offspring pairs===============#
OUT=NULL
sims=function(sims)
{
alleles2=afreqs[which(afreqs[,1]==z),]
alleles3=cbind(100:(99+(length(alleles2[,1]))),alleles2[,-1])                   #table allele frequencies
homos=(alleles3[,2])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles3[,1]),as.character(alleles3[,1]),homos)

hets=t(combn(alleles3[,2],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles3[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
n=1000000                                                                       #sample size of all simulated genotypes (customized to indidvidual data sets) #plus 1000 is to make up for shorter simulated datsets
gfreqs1=rep(gfreqs[,1],(n*(as.numeric(gfreqs[,3]))))                            #create genotypes(by coloumn, 1 for each allele)
gfreqs2=rep(gfreqs[,2],(n*(as.numeric(gfreqs[,3]))))
gtypes=cbind(gfreqs1,gfreqs2)
gtypes=gtypes[sample(1:length(gtypes[,1]),replace=F),]
sg1=gtypes[sample(1:length(gtypes[,1]),Nadults),]
OUT<<-cbind(OUT,sg1)
}
z=length(unique(afreqs[,1]))
C1=for(z in 1:z) {lapply(z,sims)}
parents=OUT
c=c(1:(ncol(OUT)))
odd=2*(unique(round(((c-2))/2)))+1
l=length(odd) * 1000
codes=seq(from=1,to=l,by=1000)
cols=sort(rep(codes,2))-1
Anumbs=matrix(cols,Nadults,ncol(OUT),byrow=T)
parents=as.numeric(parents)+Anumbs
#create full sib families (go down the list in pair)============================#
OUT2=NULL
sims=function(sims)  {
p1=parents[z,]
p2=parents[z+1,]
als=rep(1:2,length(p1)/2)
Noffs=Noffs_perpair                                                             #number of offspring per pair
OUT2=NULL
for (b in 1:Noffs){
pos1=sample(als,length(p1)/2,replace=TRUE)                                      #note that this captures the variance, could just create the 4 genotypes in equal frequencies if you dont want that variance
pos2=sample(als,length(p1)/2,replace=TRUE)
pos11=pos1+(seq(0,(length(p1)-1),2))
pos22=pos2+(seq(0,(length(p2)-1),2))
o1=p1[pos11]
o2=p2[pos22]
o3=sort(c(o1,o2))
o3=t(c(z,o3))
write.table(o3,file="SimOffs.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
}
z=length(parents[,1])
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,sims)
Dads=parents[seq(from=1,to=length(parents[,1]),by=2),]
Moms=parents[seq(from=2,to=length(parents[,1]),by=2),]
Anumbs=matrix(cols,Nadults,ncol(OUT),byrow=T)                                   #see code before functions for adding 1000s (here am removing 1000s)
parents=as.numeric(parents)-Anumbs
Dads=parents[seq(from=1,to=length(parents[,1]),by=2),]
Moms=parents[seq(from=2,to=length(parents[,1]),by=2),]
Offs2 <- read.table("SimOffs.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Offs = as.matrix(Offs2[,-1])
Anumbs=matrix(cols,length(Offs[,1]),ncol(OUT),byrow=T)
Offs=as.numeric(Offs)-Anumbs
if (Noffs_perpair>1){
Offnames=ceiling(Offs2[,1]/2)                                                   #naming of offpsring
Offnames2=paste(Offnames,".",1:length(Offnames))
Offs=cbind(paste("Offspring",Offnames2),Offs)} else {
Offnames=ceiling(Offs2[,1]/2)
Offs=cbind(paste("Offspring",Offnames),Offs)}
Moms=cbind(paste("Mom",1:length(Moms[,1])),Moms)
Dads=cbind(paste("Dad",1:length(Dads[,1])),Dads)
#add error======================================================================#
Dadsg=Dads[,-1]
ldad=length(Dadsg)*error
pdad1=sample(1:length(Dadsg),ldad,replace=FALSE)
pdad2=Dadsg[sample(1:length(Dadsg),ldad,replace=FALSE)]
Dadsg2=replace(Dadsg,pdad1,pdad2)
Dads=cbind(Dads[,1],Dadsg2)
Offsg=Offs[,-1]
loff=length(Offsg)*error
poff1=sample(1:length(Offsg),loff,replace=FALSE)
poff2=Offsg[sample(1:length(Offsg),loff,replace=FALSE)]
Offsg2=replace(Offsg,poff1,poff2)
Offs=cbind(Offs[,1],Offsg2)
#Momsg=Moms[,-1]
#ldad=length(Momsg)*error
#pdad1=sample(1:length(Momsg),ldad,replace=FALSE)
#pdad2=Momsg[sample(1:length(Momsg),ldad,replace=FALSE)]
#Momsg2=replace(Momsg,pdad1,pdad2)
#Moms=cbind(Moms[,1],Momsg2)
#===============================================================================#
colsnmz=rep("Locus",ncol(Offs)-1)
colsnmz2=sort(rep(1:(length(colsnmz)/2),2))
colsnmz3=paste(colsnmz,colsnmz2)
colsnmz3=c("IDs",colsnmz3)
colnames(Offs)<- colsnmz3
colnames(Dads)<- colsnmz3
colnames(Moms)<- colsnmz3
#calculate shared alleles among pairs of siblings=====================================================#
sibs=(Offs[,-1])
OUT2=NULL
pairs=seq(from=1,to=length(sibs[,1]),by=2)
loci=seq(from=1,to=ncol(sibs),by=2)
for (i in pairs) {
  sib1=as.numeric(sibs[i,])
  sib2=as.numeric(sibs[i+1,])
    OUT=NULL
    for (x in loci) {
      sib1a=as.numeric(sib1[x])
      sib1b=as.numeric(sib1[x+1])
      sib2a=as.numeric(sib2[x])
      sib2b=as.numeric(sib2[x+1])
      test1=sib1a-sib2a
      test2=sib1a-sib2b
      test3=sib1b-sib2a
      test4=sib1b-sib2b
      out=c(test1,test2,test3,test4)
      out2=length(which(out==0))
      if(out2>0) {out3=1}  else {out3=0}
      OUT=rbind(OUT,out3)
      }
 out4=length(which(OUT[,1]==1))
 OUT2=rbind(OUT2,out4)
}
#calculate shared alleles among parents and offspring===========================#
sibs=Offs[seq(from=1,to=length(Offs[,1]),by=2),]
sibs=sibs[,-1]
dads=Dads[,-1]
OUT3=NULL
pairs=seq(from=1,to=length(sibs[,1]),by=1)
loci=seq(from=1,to=ncol(sibs),by=2)
for (i in pairs) {
  sib1=as.numeric(sibs[i,])
  sib2=as.numeric(dads[i,])
    OUT=NULL
    for (x in loci) {
      sib1a=sib1[x]
      sib1b=sib1[x+1]
      sib2a=sib2[x]
      sib2b=sib2[x+1]
      test1=sib1a-sib2a
      test2=sib1a-sib2b
      test3=sib1b-sib2a
      test4=sib1b-sib2b
      out=c(test1,test2,test3,test4)
      out2=length(which(out==0))
      if(out2>0) {out3=1}  else {out3=0}
      OUT=rbind(OUT,out3)
      }
 out4=length(which(OUT[,1]==1))
 OUT3=rbind(OUT3,out4)
}
#calculate shared alleles among 2 randomly chosen pairs=========================#
sibs=Offs[seq(from=1,to=length(Offs[,1]),by=2),]
sibs=sibs[,-1]
dads=Dads[,-1]
OUT4=NULL
pairs=seq(from=1,to=length(sibs[,1]),by=1)
loci=seq(from=1,to=ncol(sibs),by=2)
for (i in pairs) {
  sib1=as.numeric(sibs[sample(pairs,1),])
  sib2=as.numeric(dads[sample(pairs,1),])
    OUT=NULL
    for (x in loci) {
      sib1a=sib1[x]
      sib1b=sib1[x+1]
      sib2a=sib2[x]
      sib2b=sib2[x+1]
      test1=sib1a-sib2a
      test2=sib1a-sib2b
      test3=sib1b-sib2a
      test4=sib1b-sib2b
      out=c(test1,test2,test3,test4)
      out2=length(which(out==0))
      if(out2>0) {out3=1}  else {out3=0}
      OUT=rbind(OUT,out3)
      }
 out4=length(which(OUT[,1]==1))
 OUT4=rbind(OUT4,out4)
}
results=data.frame(table(OUT4))
x3=as.numeric(as.character(results[,1]))
y3=results[,2]/sum(results[,2])
results=data.frame(table(OUT3))            #parents and offspring
x2=as.numeric(as.character(results[,1]))
y2=results[,2]/sum(results[,2])
results=data.frame(table(OUT2))
x=as.numeric(as.character(results[,1]))
y=results[,2]/sum(results[,2])
maxx=max(c(x,x2,x3))
minx=min(c(x,x2,x3))
maxy=max(c(y,y2,y3))
miny=min(c(y,y2,y3))
pdf(file="Solomon_shared_alleles.pdf")
par(cex.lab=1.2,cex=1.05,lwd=2)
plot(x,y,xlim=c(minx,maxx),ylim=c(miny,maxy),pch=21,cex=2,bg="blue",xlab="Number of Loci that share an allele (IBD+IBS)",ylab="Frequency")
abline(0,0)
abline(0.05,0,lty=2)
lines(x,y,col="blue")
lines(x2,y2,col="green")
lines(x3,y3,col="orange")
points(x,y,pch=21,cex=2,bg="blue")
points(x2,y2,pch=21,cex=2,bg="green")
points(x3,y3,pch=21,cex=2,bg="orange")
legend("topleft",c("Unrelated pairs", "Full siblings", "Parent-offspring pairs"), pch = c(21,21,21),pt.bg=c("orange","blue","green"))
dev.off()
x2=sort(x2)
x=sort(x)
x.2=x[match(x2[1],x):length(x)]
y.2=y[match(x2[1],x):length(x)]
sumy=y.2+y2
proportion=y.2/sumy
write.table(proportion,file="proportion.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=FALSE)
unlink("SimOffs.txt")
#Begin modified Bayesian parentage analysis=====================================#
#(modified in last section to account for IBD alleles)

  Adults1=adults
  Offspring1=offspring
  wanted_reps=as.numeric(wanted_reps)
  loci=ncol(Adults1)
  Adults=Adults1[,c(2:loci)]                                                    #assumes that there is one column of id names ; could modify this as needed.
  Offspring=Offspring1[,c(2:loci)]
  total <- ncol(Adults)/2                                                       #For Progress bar
  if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = total, width = 300)}
afreqs=function(afreqs)                                                         #Begin Master simulation function
{
locus_name=L
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, ceiling(locus_name/2), title=paste("Locus", ceiling(locus_name/2),"of ",total))}
vect=c(Adults[,L],Adults[,L+1])                                                 #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
alleles=data.frame(table(vect))
alleles=alleles[order(alleles[,1]),]
if (as.numeric(as.character(alleles[1,1]))==0) {alleles=alleles[-1,]}           #removes missing data
alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
homos=(alleles2[,3])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
hets=t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
csum=cumsum(as.numeric(gfreqs[,3]))
gfreqs1=cbind(gfreqs,csum)
Nadults=length(Adults[,1])
Noffs=length(Offspring[,1])
#===============================================================================#end locus-specific HWE genotype frequency calculations
alength=length(alleles2[,1])
for (Y in 1:wanted_reps) {
positions=1:length(gfreqs[,1])
sg3=sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                      #sample the repeated genotype positions, by the number of adults
sadults=gfreqs[sg3,1:2]                                                         #index gfreqs to create genotypes
og3=sample(positions,Noffs,replace=TRUE,prob=gfreqs[,3])                        #create juvenile genotyes
soffs=gfreqs[og3,1:2]
soffs=cbind(as.numeric(soffs[,1]),as.numeric(soffs[,2]))
asp=cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
asp=rbind(cbind(locus_name,0,0),asp)
for (i in 1:Nadults) {
parent1=as.numeric(sadults[i,1])                                                #first allele in parent
parent2=as.numeric(sadults[i,2])                                                #second allele in parent
p1=soffs-parent1
p2=soffs-parent2
pp1=which(p1[,1]==0)
pp2=which(p1[,2]==0)
allele1=unique(c(pp1,pp2))
p21=which(p2[,1]==0)
p22=which(p2[,2]==0)
allele2=unique(c(p21,p22))
Out51=cbind(parent1,length(allele1))
Out52=cbind(parent2,length(allele2))
Out53=cbind(0,Noffs-length(unique(c(allele1,allele2))))
Out5=rbind(Out51,Out52,Out53)
if(parent2==parent1) {Out5=Out5[-1,]}                                           #remove 1 of alleles count if homozygous
if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
  diffs=sum(Out5[,2])-Noffs                                                     #comment out to be more conservative!
  maxa=max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
  pos=which(Out5[,2]==maxa)
  Out5[pos,2]<-Out5[pos,2]-diffs}
m1=match(Out5[,1],asp[,2])
m2=asp[m1,3]+as.numeric(Out5[,2])
asp[m1,3]<-m2
asp<-asp
}
write.table(asp,file="out.sims",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}
}
L=ncol(Adults)
C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
if(.Platform$OS.type=="windows"){close(pb)}
#Bayes P-value==================================================================#
OUT<- read.table("out.sims", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
locname=unique(OUT[,1])                                                         #compile calculations for each locus
OUTALL=NULL
for (z in locname) {
  Loc1=OUT[which(OUT[,1]==z),]
  allfreqs=unique(Loc1[,2])
  OUT2=NULL
    for (x in allfreqs) {
    a1<-Loc1[which(Loc1[,2]==x),]
    a2=sum(a1[,3])
    a3=cbind(x,a2)
    OUT2 <- rbind(OUT2, a3)
  }
  OUT3=cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
  OUTALL <- rbind(OUTALL, cbind(z,OUT3))
}
#Create multilocus genotypes, 50000 at a time, calculate number of loci that mismatch, and calculate freqs of shared alleles
NL=length(unique(OUTALL[,1]))
ngtypes=1                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
Ntotal=as.numeric(Ntotal)
inreps=50000     #was tested as 10000 for SNPS
repnumber=round(Ntotal/inreps)                                                #this is the rep number to get 100,000 values.  can adjust accordingly
asp=cbind(0:NL,rep(0,length(0:NL)))
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = Ntotal , width = 300)}
for (n in 1:repnumber) {
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, n*inreps, title=paste("Genotype", n*inreps,"of ",Ntotal))}
    OUT=NULL
    for (r in unique(OUTALL[,1])) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*Loco[,4])
        findo=which(alleles3[,2]==0)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
  distm=apply(OUT, 1, function(x)sum(x == 1))
  distm2=data.frame(table(distm))
  m1=match(distm2[,1],asp[,1])
  m2=asp[m1,2]+distm2[,2]
  asp[m1,2]<-m2
  asp<-asp
}
if(.Platform$OS.type=="windows"){close(pb)}
#tabulate number of multilocus genotypes with x mismatches
d2=asp
d3=cbind(d2,d2[,2]/sum(d2[,2]))
#===============================================================================# Create plot of exclusionary power.  Also, some necessary data formatting
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "Posterior Processing", min = 0, max = 7 , width = 300)}
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 1, title=paste("Performing Exclusion"))}
Adults<- adults
Offs<- offspring
Nads=length(Adults[,1])
Noffs=length(Offs[,1])
asp=d3
asp=cbind(asp,NL-as.numeric(as.character(asp[,1])))  #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
asp=cbind(asp,cumsum(asp[,2]))
asp=cbind(asp,asp[,5]/max(asp[,5]))
#find minimum Nloci to mismatch (could modify this) and perform exclusion with decided-upon mismatches==#
distm=cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
#mismatch=min(which(round(distm[,6],1)==.9))                                    #deprecated
mismatch=which.min(abs(distm[,6] - .5))                                         #could cause issues with a value of .5 chosen here - could be too low.  Change to higher if all loci analyzed have phi < 1.
Adults1=Adults                                                                   #begin exclusion
Offspring1=Offs
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Offspring=Offspring1[,c(2:a)]
Anames=Adults1[,1]
Onames=Offspring1[,1]
Adults=Adults1[,-1]
Offspring=Offspring1[,-1]
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
A=1:Aindivids
O=1:Oindivids
G=expand.grid(A,O)
AG=G[,1]
AO=G[,2]
Ads=Adults[AG,]
Offs=Offspring[AO,]
IdnamesA=Anames[AG]
IdnamesO=Onames[AO]
write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Names=cbind(IdnamesA,IdnamesO)
matches=function(matches)
{
A=Ads[,z]-Offs[,z]
B=Ads[,(z+1)]-Offs[,(z+1)]

C=Ads[,z]-Offs[,(z+1)]
D=Ads[,(z+1)]-Offs[,z]
f=A*B*C*D
f=(f^2)*10
ss=which(is.na(f)==TRUE)
f=replace(f,ss,0)

identify=which(f>0)
f=replace(f,identify,1)
f=cbind(z,f)
write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Ads)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
a=length(U[,1])
stuff=rowSums(U)
Sorted=cbind(Names,stuff)
matches=which(Sorted[,3]<(mismatch+1))
Actual=sort(stuff)
IDS=which(stuff<(mismatch+1))
PAdults=Ads[IDS,]
POffspring=Offs[IDS,]
nput=length(matches)
Putativepairs=Sorted[matches,]
PAdults=cbind(Putativepairs[,c(1,3)],PAdults)
POffspring=cbind(Putativepairs[,c(2,3)],POffspring)
names(PAdults)[1]<-"ID"
names(POffspring)[1]<-"ID"
sorts=function(sorts)
{
tell=rbind(PAdults[f,],POffspring[f,])
write.table(tell,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
f=length(PAdults[,1])
C1= for(f in 1:f) lapply(f,sorts)
unlink("Sort.txt")
unlink("IdnamesA.txt")
unlink("IdnamesO.txt")
#Calculate phi for each number mismatching loci=================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 2, title=paste("Calculating phi")) }
Putative<- read.table("Output_genotypes.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
observed=data.frame(table(Putative[,2]))
observed=cbind(observed,observed[,2]/2)                                         #done becuase each number is written twice in file (once for parent and once for off)
zerom=0:mismatch                                                                #this chunk adds 0s for mismatches where there were no observed putative pairs
zerom2=which(is.na(match(zerom,observed[,1])))
if (length(zerom2>0)) {observed=observed[,3]
for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}  #not really 0, to prevent divide by 0
}   else {observed=observed[,3]}
expected=distm[1:(mismatch+1),7]                                                #using cumulatinve sum   (more conservative)
#expected=distm[1:(mismatch+1),3]*Nads*Noffs                                     #not using cumulative sum
phi=expected/observed
phi=replace(phi,which(phi>=1),1)
Offs<- offspring
actualTrue=length(grep("Off",Offs[,1]))
info=cbind(actualTrue,expected,observed,phi)
#calculate phi and index values ================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 3, title=paste("Creating distriubtion of falsely-shared alleles")) }
phibase=phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
observed=observed[min(which(phi[]==1))-1]                                       #do the same for observed
testob=which(observed==0.000000001)
phi2=cbind(1:length(phi),phi)
if (length(testob>0)) {phi4=phi2[-testob,]} else {phi4=phi2}
nmismatch=min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
index=phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
index=index[which(index>-1)]
if (length(index)>1) {
  if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
phi=phi[index]
index=index-1
#Create Plot ===================================================================#
pdf(file="Output_Dataset_Power.pdf")
x=0:(length(info[,1])-1)
y1=info[,3]
y2=info[,2]
p1=which(y2==0)
y2=replace(y2,p1,.000000001)
p1=which(y1<y2)
y3=replace(y1,p1,y2[p1])
y2=y2+1
y3=y3+1
par(mar=c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
ats=c(0,1,2,3,4,5,6,7,8,9)
labs=c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
axis(side=2,ats,labs)
lines(x,log10(y3),lwd=2)
lines(x,log10(y2),lwd=2)
points(x,log10(y3),pch=21,bg="green",cex=2)
points(x,log10(y2),pch=21,bg="blue",cex=2)
legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
yphi=y2/y3
par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
points(x,yphi,cex=2,pch=21,bg="gray")
mtext("Number of Mismatching Loci",side=1,line=1.94)
dev.off()
info2=cbind(x,info[,4])
colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
#===============================================================================#
ngtypes=20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
inreps=10000
repnumber=round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
#writes values at all loci, to be analyzed further below
for (n in 1:repnumber) {
    OUT=NULL
    for (r in unique(OUTALL[,1])) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*Loco[,4])
        findo=which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        alleles3=alleles3[-which(alleles3[,2]==0),]
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
for (i in index) {                                                              #loop over numbers of mismatched loci
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                    #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distp=cbind(i,DIST)
    write.table(distp,file="False_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
 }
}

#Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
OUT9=NULL
for (n in index){
  Putative2=Putative[which(Putative[,2]==n),]
  write.table(Putative2, file="Putative.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    #if (length(distp[,1])>1000) {                                                #need at least 1000 values , else assigned 0 (may be redundant now)
       OBS=NULL
       afreqs=function(afreqs) {
        PutL=Putative2[,c(L+2,L+3)]
        PutLadults=PutL[seq(from=1,to=length(PutL[,1]),by=2),]                  #Combine putative pairs alleles into a single row
        PutLoffs=PutL[seq(from=2,to=length(PutL[,1]),by=2),]
        Puts=cbind(PutLadults,PutLoffs)
        c1=Puts[,1]-Puts[,3]                                                    #find the matching alleles
        c2=Puts[,1]-Puts[,4]
        c3=Puts[,2]-Puts[,3]
        c4=Puts[,2]-Puts[,4]
        Puts2=cbind(Puts,c1,c2,c3,c4)
        P5=replace(Puts2[,5],which(Puts2[,5]!=0),-10)
        P6=replace(Puts2[,6],which(Puts2[,6]!=0),-10)
        P7=replace(Puts2[,7],which(Puts2[,7]!=0),-10)
        P8=replace(Puts2[,8],which(Puts2[,8]!=0),-10)
        P5=replace(P5,which(P5==0),1)
        P6=replace(P6,which(P6==0),1)
        P7=replace(P7,which(P7==0),1)
        P8=replace(P8,which(P8==0),1)
        Puts3=cbind(Puts,P5,P6,P7,P8)
        Puts4=cbind((Puts3[,1]*Puts3[,5]),(Puts3[,1]*Puts3[,6]),(Puts3[,2]*Puts3[,7]),(Puts3[,2]*Puts3[,8]))
        alleles2=OUTALL[which(OUTALL[,1]==L),]
        alfreq1=alleles2[match(Puts4[,1],alleles2[,2]),4]
        alfreq2=alleles2[match(Puts4[,2],alleles2[,2]),4]                       #find the actual allele values
        alfreq3=alleles2[match(Puts4[,3],alleles2[,2]),4]
        alfreq4=alleles2[match(Puts4[,4],alleles2[,2]),4]
        Puts5=cbind(alfreq1,alfreq2,alfreq3,alfreq4)                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
        R1=replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
        R2=replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
        R3=replace(Puts5[,3],which(is.na(Puts5[,3])==TRUE),1)
        R4=replace(Puts5[,4],which(is.na(Puts5[,4])==TRUE),1)
        Puts6=cbind(R1,R2,R3,R4)
        Put_share=apply(Puts6, 1, min)                                          #find row minimum
        Put_share2=apply(Puts4, 1, max)                                         #find shared allele name
        Put_share3=c(Put_share,Put_share2)
        OBS <<- cbind(OBS,Put_share3)
    }
    L=ncol(Putative2)-2
    C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    lengths=length(OBS[,1])/2
    if (lengths==1) {OBA3=t(OBS[(lengths+1):(2*lengths),])} else {OBA3=OBS[(lengths+1):(2*lengths),]}
    write.table(OBA3,file="True_shared_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)  #Actual shared alleles   #This file wouuld be useful as output for people to identify shared and mismatching loci
    if (lengths==1) {OBS=t(OBS[1:lengths,])} else {OBS=OBS[1:lengths,]}   #formatting for if there is only a single pair
    obsp=apply(OBS, 1, prod)
    #}  else obsp=rep(0,(length(Putative2[,1]))/2)
OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
}
#calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic)==#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 4, title=paste("Calculating Posterior Component 1"))}
OBA3<- read.table("True_shared_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
for (n in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
      OUT=NULL
      for (r in unique(OUTALL[,1])) {                                           #This first section calculates products of parental alleles (True distribution)
        vect=c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
        alleles=data.frame(table(vect))
        alleles=alleles[order(alleles[,1]),]
        if (as.numeric(as.character(alleles[1,1]))<=0) {alleles=alleles[-1,]}
        alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))
        gtrue=sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
        OUT <- cbind(OUT,gtrue)

        if(n==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele
                            mm=alleles2[match(OBA3[i,ceiling(r/2)],alleles2[,1]),3]
                            write.table(cbind(r,mm),file="Shared_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
                            }
                }
        }
    for (i in index) {
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                   #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distt<-cbind(i,DIST)
    write.table(distt,file="True_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
}
#Calculate lamdaphi=============================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 5, title=paste("Calculating Posterior Component 2"))}
Putative3<- read.table("Putative.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Putadults=Putative3[seq(from=1,to=length(Putative3[,1]),by=2),1]                #Combine putative pairs alleles into a single row
Putoffs=Putative3[seq(from=2,to=length(Putative3[,1]),by=2),1]
Names=cbind(as.character(Putadults),as.character(Putoffs))
empirical=cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
distp <- read.table("False_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
P2=NULL
for (i in index) {                                                              #loop over numbers of mismatched loci
  empirical2=empirical[which(as.numeric(as.character(empirical[,3]))==i),]
  if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
  if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
    a3=distp[which(distp[,1]==i),2]
    DIST<-as.data.frame(a3)
    P=NULL
    for (b in 1:length(empirical2[,1])) {
      p1=length(which(DIST[,1] <= as.numeric(empirical2[b,4]) ))
      if (p1==0) {p1=0.00001}
      p2=cbind(empirical2[b,1],empirical2[b,2],p1)
      p3=cbind(p2,as.numeric(p2[,3])/length(DIST[,1]))
      P <- rbind(P,p3)
      }
    }
  P2<-rbind(P2,cbind(i,P))
}
lamdaphi=as.numeric(as.character(P2[,5]))
#Calculate lamda|phic ==========================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 6, title=paste("Calculating Posterior Component 3"))}
lamdaphic_dist<- read.table("True_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Observed<- read.table("Shared_allele_freqs.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
mm1=which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
Observed[mm1,2] <- 1                                                            #replace NAs with 1
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
lamdaphic=apply(U, 1, prod)
l1=length(which(OUT9[,2]==0))
if (l1>0) lamdaphic=c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)    #double check values by hand
P3=cbind(P2,lamdaphic)
for (i in index) {                                                              #loop over numbers of mismatched loci
   e2=P3[which(as.numeric(as.character(P3[,1]))==i),]
   if(length(e2)==6) {e2=t(e2)}                                                 #deals with one indiviudal formatting
   if (e2[1,5]==0) {write.table(e2[,6], file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{                    #deals with not enough reps
    a3=lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
    DIST<-as.data.frame(a3)
  for (b in 1:length(e2[,1])) {                                                 #calculate p values
    p1=length(which(DIST[,1] <= e2[b,6]))
    p2=p1/length(DIST[,1])
    write.table(p2, file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    }
  }
}
lamdaphic<- read.table("lamdaphic.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
#Put it all together with Bayes theorem!========================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 7, title=paste("Calculating Posterior"))}
vals=cbind(P2,lamdaphic[,1])
phires=phi
sib<- read.table("proportion.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
sib=rev(t(sib))
nphi=length(phi)
sib=sib[1:nphi]
#sib2=sib+1
#phi2=phi*sib2
phi2=phi+sib
phi=phi2                                                                        #add in alleles that are IBD
phi[which(phi>1)]<-1
philength=cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
phis=rep(philength[,2],philength[,3])
vals=cbind(vals,phis)
colnames(vals)<-c("Nmismatch","Parent","Off","ignore","lamdaphi","lamdaphic","phi")
phi=as.numeric(as.character(vals[,7]))
lamdaphi=as.numeric(as.character(vals[,5]))
lamdaphic=as.numeric(as.character(vals[,6]))
lamdaphi=replace(lamdaphi,which(lamdaphi==0),1)
lamdaphic=replace(lamdaphic,which(lamdaphic==0),1)
pval=(lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
pval=cbind(vals[,2],vals[,3],vals[,1],pval)
pval=pval[order(as.numeric(pval[,4])),]
colnames(pval) <- c("Adult","Offspring","NL_mismatch","Probability of pair being false given frequencies of shared alleles")
#write.table(vals, file="Posterior_Components_Bayes.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
if(.Platform$OS.type=="windows"){close(pb)}
unlink("False_allele_freqs.txt")                                                #clear all sims files
unlink("True_allele_freqs.txt")
unlink("Shared_allele_freqs.txt")
unlink("lamdaphic.txt")
unlink("Putative.txt")
unlink("True_shared_freqs.txt")
unlink("Output_genotypes.txt")
unlink("proportion.txt")
unlink("AlleleFrequencies_SOLOMON.txt")
unlink("*.sims")
rm(list=ls())
}
label.bayes <- tklabel(tt, text="Press 'Run' to perform Bayesian parentage analysis:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label.bayes,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label.bayes, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
label.SIMS <- tklabel(tt, text="Perform Bayesian parentage analysis accounting for siblings:",font=fontTextLabel)
OK.but.SIMS <- tkbutton(tt, text = "SIBS", command = PressedOKSIB)
tkgrid(tklabel(tt,text="     "),label.SIMS,tklabel(tt,text="     "),OK.but.SIMS,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.SIMS, sticky="w")
tkgrid.configure(OK.but.SIMS, sticky="w")
tkgrid(tklabel(tt,text="     "))
#Module7_exclusion_1parent######################################################################################################################################################################################################################
#Exclusion_OneKnownParent_Module
heading <- tklabel(tt, text="Parentage with 1 Known Parent",font=fontSUB)       # add a heading
tkgrid(heading, columnspan=5)
PressedEXCLUSION1 <- function()
{
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Parentage Analysis with 1 Known Parent"                          # Name the window
heading <- tklabel(tt, text="Exclusion",font=fontHeading)                       # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set number of mismatching loci ================================================#
label_mismatches <- tklabel(tt, text="Set number of loci to let mismatch:",font=fontTextLabel)
Name1 <- tclVar("0")
entry.label2 <- tkentry(tt, width="4",textvariable=Name1)                       #create entry fields
OnOK2 <- function()
{
	NameVal <- tclvalue(Name1)
 	msg <- paste("You have now set the number of mismatches to",NameVal)
	tkmessageBox(message=msg)
	assign("mismatch", NameVal, envir = solomon.env)
}
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label2, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_mismatches, tklabel(tt,text="     "), entry.label2, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_mismatches, sticky="w")
tkgrid.configure(entry.label2,sticky="w")
#Set working Directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()
{
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Load Known Adults File=========================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings=-900, dec=".", strip.white=TRUE)
  assign("MOMS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Known Parents File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing known-parent genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Putative Adults File======================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings=-900, dec=".", strip.white=TRUE)
  assign("DADS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Putative Parent File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing putative-parent genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")

tkgrid.configure(adults.button,sticky="w")
#Load Offspring File============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings=-900, dec=".", strip.white=TRUE)
  assign("OFFSPRING", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
#Create button to run parentage script==========================================#
PressedOK <- function()
    {
#begin exclusion to known mothers
moms<- MOMS
Adults1=moms
Offs<- OFFSPRING
mismatch=0
Offspring1=Offs
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Offspring=Offspring1[,c(2:a)]
Anames=Adults1[,1]
Onames=Offspring1[,1]
Adults=Adults1[,-1]
Offspring=Offspring1[,-1]
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
A=1:Aindivids
O=1:Oindivids
G=cbind(A,O)
AG=G[,1]
AO=G[,2]
Ads=Adults[AG,]
Offs=Offspring[AO,]
IdnamesA=Anames[AG]
IdnamesO=Onames[AO]
write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Names=cbind(IdnamesA,IdnamesO)
matches=function(matches)
{
A=Ads[,z]-Offs[,z]
B=Ads[,(z+1)]-Offs[,(z+1)]
C=Ads[,z]-Offs[,(z+1)]
D=Ads[,(z+1)]-Offs[,z]
f=A*B*C*D
f=(f^2)*10
ss=which(is.na(f)==TRUE)
f=replace(f,ss,0)
identify=which(f>0)
f=replace(f,identify,1)
f=cbind(z,f)
write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Ads)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
a=length(U[,1])
stuff=rowSums(U)
Sorted=cbind(Names,stuff)
matches=which(Sorted[,3]<(mismatch+1))
Actual=sort(stuff)
IDS=which(stuff<(mismatch+1))
PAdults=Ads[IDS,]
POffspring=Offs[IDS,]
nput=length(matches)
Putativepairs=Sorted[matches,]
PAdults=cbind(Putativepairs[,c(1,3)],PAdults)
POffspring=cbind(Putativepairs[,c(2,3)],POffspring)
names(PAdults)[1]<-"ID"
names(POffspring)[1]<-"ID"
sorts=function(sorts)
{
tell=rbind(PAdults[f,],POffspring[f,])
write.table(tell,file="Putative2.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
f=length(PAdults[,1])
C1= for(f in 1:f) lapply(f,sorts)
unlink("Sort.txt")
unlink("IdnamesA.txt")
unlink("IdnamesO.txt")
#===============================================================================#
#Begin exclude maternal allele
Putative<- read.table("Putative2.txt", header=F, sep="\t", na.strings=NA, dec=".", strip.white=TRUE) #make
Putative=Putative[,-2]
#Putative=as.matrix(Putative)  #new
#mpositions=which(Putative[]=="NA")  #new
#Putative[mpositions]<-900             #new
#Putative=as.data.frame(Putative)

categories=ncol(Putative[,-1])
Aindivids=length(Putative[,1])
c=c(1:categories)
odd=2*(unique(round(((c-2))/2)))+1
l=length(odd) * 1000
codes=seq(from=1,to=l,by=1000)
cols=sort(rep(codes,2))-1
Anumbs=matrix(cols,Aindivids,categories,byrow=T)
Putative2=Putative[,-1]+Anumbs
Putative=cbind(Putative[,1],Putative2)
allelefinder=function(allelefinder)
{
names=cbind(as.character(Putative[z,1]),as.character(Putative[z+1,1]))
put1=Putative[z,-1]+100000
put2=Putative[z+1,-1]+100000
asp=which(is.na(put1))
put1[asp]<-asp[1]*1000
put2[which(is.na(put2)==TRUE)]<-asp[1]*1000
matches=match(put1,put2)
matches2=match(put2,put1)
finder=which(matches>0)
finder2=which(matches2>0)
allelef=cbind(put2[,1],put1[finder])
reps=table(as.numeric(signif(allelef[-1],3)))
remov=which(reps>1)
seeit=(which(reps==1))
alleles=allelef[-1]
alength=length(reps)
rdata=seq(from=1,to=alength,by=1)
dstack=cbind(reps,rdata)
dstackf=rep(dstack[,2],dstack[,1])
remfin=which(match(dstackf,remov)>0)
sees=which(match(dstackf,seeit)>0)
allhet=alleles[remfin]
allssing=alleles[sees]
alldup=sort(as.integer(rep(allssing,2)))
allallele=sort(as.integer(c(allhet,alldup)))
allallele=unlist(allallele)
put3=put2
detr=cbind(as.data.frame(sort(allallele)),sort(t(put3)))
dets=detr[,1]-detr[,2]
detf=which(dets==0)
detd=detr[-detf,2]
alla=allallele
dets=dets*100
da=which(dets>0)
db=which(dets<0)
dc=c(da,db)
dc=ceiling(dc/2)
dd=c((dc*2),((dc*2)-1))
newa=rep(detd,2)
aba=allallele[-dd]
alleles=c(aba,newa)
alleles=sort(alleles)
locn=((ncol(put1)))/2
if(length(alleles)>0) {                                                           #if no alleles get excluded (e.g., all homozygous)
locd=seq(from=100000, to=floor(max(alleles)), by=1000)                          #!
removdiv=sort(rep(locd,2))
if(length(finder)==ncol(put1)) final=t(as.data.frame(rep(0,ncol(put1)))) else {
alla=as.data.frame(alleles)
aldiv=cbind(alla,removdiv)
final=aldiv[,1]-aldiv[,2]
final=t(final)}
final=cbind(names,final)
} else {final=cbind(names,Putative[z+1,-1]-Anumbs[1,])}
write.table(final,file="sharedallelestmp.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}
z=length(Putative[,1])
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,allelefinder)
#Begin exclusion to unknown parent==============================================#
mismatch=as.numeric(mismatch)                                                   #used because it is entered as a character through GUI
Adults1<- DADS
loci=ncol(Adults1)
Adults=Adults1[,c(2:loci)]                                                      #assumes that there is one column of id names ; could modify this as needed.
Offspring1<- read.table("sharedallelestmp.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
write.table(Offspring1,file="sharedalleles.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=FALSE)
Offspring=Offspring1[,c(3:(loci+1))]
Adults<- DADS
Offs<- read.table("sharedallelestmp.txt", header=FALSE, sep="\t", na.strings="0", dec=".", strip.white=TRUE)
unlink("sharedallelestmp.txt")
Nads=length(Adults[,1])
Noffs=length(Offs[,1])
Adults1=Adults
Offspring1=Offs
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Offspring=Offspring1[,c(3:(a+1))]                                               #special formatting for this script only!
colnames(Adults)<-""
colnames(Offspring)<-""
Anames=Adults1[,1]
Onames=Offspring1[,2]
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
A=1:Aindivids
O=1:Oindivids
G=expand.grid(A,O)
AG=G[,1]
AO=G[,2]
Ads=Adults[AG,]
Offs=Offspring[AO,]
IdnamesA=Anames[AG]
IdnamesO=Onames[AO]
write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Names=cbind(IdnamesA,IdnamesO)
matches=function(matches)
{
A=Ads[,z]-Offs[,z]
B=Ads[,(z+1)]-Offs[,(z+1)]
C=Ads[,z]-Offs[,(z+1)]
D=Ads[,(z+1)]-Offs[,z]
f=A*B*C*D
f=(f^2)*10
ss=which(is.na(f)==TRUE)
f=replace(f,ss,0)
identify=which(f>0)
f=replace(f,identify,1)
f=cbind(z,f)
write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Ads)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
a=length(U[,1])
stuff=rowSums(U)
Sorted=cbind(Names,stuff)
matches=which(Sorted[,3]<(mismatch+1))
Actual=sort(stuff)
IDS=which(stuff<(mismatch+1))
PAdults=Ads[IDS,]
POffspring=Offs[IDS,]
nput=length(matches)
Putativepairs=Sorted[matches,]
PAdults=cbind(Putativepairs[,c(1,3)],PAdults)
POffspring=cbind(Putativepairs[,c(2,3)],POffspring)
names(PAdults)[1]<-"ID"
names(POffspring)[1]<-"ID"
sorts=function(sorts)
{
tell=rbind(PAdults[f,],POffspring[f,])
write.table(tell,file="Output_genotypestmp.txt",row.names=FALSE,col.names=F,sep="\t",append=TRUE)
}
f=length(PAdults[,1])
C1= for(f in 1:f) lapply(f,sorts)
unlink("Sort.txt")
unlink("IdnamesA.txt")
unlink("IdnamesO.txt")
unlink("Putative2.txt")
#Outputfile processing==========================================================#
putative<- read.table("Output_genotypestmp.txt", header=FALSE, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
write.table(putative,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=FALSE)
unlink("Output_genotypestmp.txt")
parents=putative[seq(from=1,to=length(putative[,1]),by=2),c(1,2)]
offspring=putative[seq(from=2,to=length(putative[,1]),by=2),c(1,2)]
out1=cbind(parents,offspring)
out2=out1[order(out1[,1]),-2]
out2=cbind(out2,"")
out2[,1] <- as.character(out2[,1])
out2[,4] <- as.character(out2[,4])
ids=2:length(out2[,1])
for (x in 2:length(out2[,1])) {                                                 #using 2 here (because 1 minus nothing is invalid : add in first parent later
    if (out2[x,1]!=out2[x-1,1]) {out2[x,4] <- out2[x,1]}
}
out2[1,4] <- out2[1,1]
out3=out2[,c(4,2,3)]
out4=cbind(out3,"")
out4[,4] <- as.character(out4[,4])
count=data.frame(table(out2[,1]))
for (x in 1:length(count[,1])){
out4[match(count[x,1],out3[,1]),4]=count[x,2]
}
noffs=sort(unique(out4[,4]))
if (noffs[1]=="") {noffs=noffs[-1]}
noffs=sort(as.numeric(noffs),decreasing=TRUE)
noffs=as.character(noffs)                                                       #comment out this line if want highest reproductive success at bottom of file
OUT=NULL
for (x in 1:length(noffs)){
  sim1=which(out4[,4]==noffs[x])
      VALS=NULL
      for (y in 1:length(sim1)){                                                #possibly start at 2 not 1
      vals=seq(sim1[y],sim1[y]+(as.numeric(noffs[x])-1),1)
      vals=cbind(noffs[x],vals)
      VALS=rbind(VALS,vals)
      }
  sim2=out4[as.numeric(VALS[,2]),]
  OUT=rbind(OUT,sim2)
}

parent.centric=OUT[,c(1,4,2,3)]
colnames(parent.centric)<-c("Parent","Number_Offspring","Offspring","Number_Loci_mismatching")
write.table(parent.centric,file="Output_by_parent.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
#Begin Ouputfiles_offspringcentric==============================================#
parents=putative[seq(from=1,to=length(putative[,1]),by=2),c(1,2)]
offspring=putative[seq(from=2,to=length(putative[,1]),by=2),c(1,2)]
out1=cbind(parents,offspring)
out2=out1[,-2]
out2=cbind(out2,"")
out2[,2] <- as.character(out2[,2])
out2[,4] <- as.character(out2[,4])
out2=out2[order(out2[,2]),]
ids=2:length(out2[,2])
for (x in 2:length(out2[,2])) {                                                 #using 2 here (because 1 minus nothing is invalid : add in first parent later
    if (out2[x,2]!=out2[x-1,2]) {out2[x,4] <- out2[x,2]}
}
out2[1,4] <- out2[1,2]
out3=out2[,c(4,1,3)]
out4=cbind(out3,"")
out4[,4] <- as.character(out4[,4])
count=data.frame(table(out2[,2]))
count=count[which(count[,2]>=1),]
for (x in 1:length(count[,1])){
out4[match(count[x,1],out3[,1]),4]=count[x,2]
}
noffs=sort(unique(out4[,4]))
if (noffs[1]=="") {noffs=noffs[-1]}
noffs=sort(as.numeric(noffs),decreasing=TRUE)
noffs=as.character(noffs)                                                       #comment out this line if want highest reproductive success at bottom of file
OUT=NULL
for (x in 1:length(noffs)){
  sim1=which(out4[,4]==noffs[x])
      VALS=NULL
      for (y in 1:length(sim1)){
      vals=seq(sim1[y],sim1[y]+(as.numeric(noffs[x])-1),1)
      vals=cbind(noffs[x],vals)
      VALS=rbind(VALS,vals)
      }
  sim2=out4[as.numeric(VALS[,2]),]
  OUT=rbind(OUT,sim2)
}
parent.centric=OUT[,c(1,4,2,3)]
colnames(parent.centric)<-c("Offspring","Number_Parents","Parents","Number_Loci_mismatching")
write.table(parent.centric,file="Output_by_offspring.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
}
label_exclusion <- tklabel(tt, text="Press 'Run' to perform exlcusion:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label_exclusion,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label_exclusion, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
label.EXCLUSION1 <- tklabel(tt, text="Perform exclusion with 1 known parent:",font=fontTextLabel)
OK.but.EXCLUSION1 <- tkbutton(tt, text = "EXCLUSION", command = PressedEXCLUSION1)
tkgrid(tklabel(tt,text="     "),label.EXCLUSION1,tklabel(tt,text="     "),OK.but.EXCLUSION1,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.EXCLUSION1, sticky="w")
tkgrid.configure(OK.but.EXCLUSION1, sticky="w")
tkfocus(tt)
#Module8_Bayes_1parent#############################################################################################################################################################################################################
#Bayes_OneKnownParent_Module
PressedBAYES1 <- function()
{
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Bayesian Parentage Analysis with 1 Known Parent"                          # Name the window
heading <- tklabel(tt, text="SOLOMON: Bayesian Parentage Analysis",font=fontHeading)   # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set working Directory==========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()
{
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Set Number of Sims ============================================================#
label_sim.number <- tklabel(tt, text="Number of simulated data sets:",font=fontTextLabel)
Name.sim <- tclVar(1000)
entry.label.sim <- tkentry(tt, width="10",textvariable=Name.sim)                #create entry fields
OnOK2 <- function()  {
	NameVal <- tclvalue(Name.sim)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("wanted_reps", NameVal, envir = solomon.env)
  }
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label1, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_sim.number, tklabel(tt,text="     "), entry.label.sim, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special1=tklabel(tt, text="Recommended: Microsatellites=1000 SNPs=100",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special1, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number, sticky="w")
tkgrid.configure(entry.label.sim,label_Special1, sticky="w")
#Set Number of "Geonotypes"=============================================================#
label_Ntotal <- tklabel(tt, text="Number of simulated genotypes:",font=fontTextLabel)
Name.Ntotal <- tclVar(50000000)
entry.label.Ntotal <- tkentry(tt, width="10",textvariable=Name.Ntotal)                #create entry fields
OnOKN <- function()  {
	NameVal <- tclvalue(Name.Ntotal)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("Ntotal", NameVal, envir = solomon.env)
  }
OK.butN <-tkbutton(tt,text="   OK   ",command=OnOKN)
tkbind(entry.label1, "<Return>",OnOKN)
tkgrid(tklabel(tt,text="     "), label_Ntotal, tklabel(tt,text="     "), entry.label.Ntotal, tklabel(tt,text="     "), OK.butN, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special=tklabel(tt, text="Recommended: Microsatellites=50,000,000 SNPs=500,000",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_Ntotal, sticky="w")
tkgrid.configure(entry.label.Ntotal,label_Special, sticky="w")
#Load Known Adults File ========================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("MOMS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Known Parents File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing known-parent genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Putative Adults File =====================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("DADS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Putative Parent File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing putative-parent genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Offspring File============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("OFFSPRING", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
#Create button to run parentage script =========================================#
PressedOK <- function()
    {
#begin exclusion to known mothers
moms<- MOMS
Adults1=moms
Offs<- OFFSPRING
mismatch=0
Offspring1=Offs
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Offspring=Offspring1[,c(2:a)]
Anames=Adults1[,1]
Onames=Offspring1[,1]
Adults=Adults1[,-1]
Offspring=Offspring1[,-1]
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
A=1:Aindivids
O=1:Oindivids
G=cbind(A,O)
AG=G[,1]
AO=G[,2]
Ads=Adults[AG,]
Offs=Offspring[AO,]
IdnamesA=Anames[AG]
IdnamesO=Onames[AO]
write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Names=cbind(IdnamesA,IdnamesO)
matches=function(matches)
{
A=Ads[,z]-Offs[,z]
B=Ads[,(z+1)]-Offs[,(z+1)]
C=Ads[,z]-Offs[,(z+1)]
D=Ads[,(z+1)]-Offs[,z]
f=A*B*C*D
f=(f^2)*10
ss=which(is.na(f)==TRUE)
f=replace(f,ss,0)
identify=which(f>0)
f=replace(f,identify,1)
f=cbind(z,f)
write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Ads)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
a=length(U[,1])
stuff=rowSums(U)
Sorted=cbind(Names,stuff)
matches=which(Sorted[,3]<(mismatch+1))
Actual=sort(stuff)
IDS=which(stuff<(mismatch+1))
PAdults=Ads[IDS,]
POffspring=Offs[IDS,]
nput=length(matches)
Putativepairs=Sorted[matches,]
PAdults=cbind(Putativepairs[,c(1,3)],PAdults)
POffspring=cbind(Putativepairs[,c(2,3)],POffspring)
names(PAdults)[1]<-"ID"
names(POffspring)[1]<-"ID"
sorts=function(sorts)
{
tell=rbind(PAdults[f,],POffspring[f,])
write.table(tell,file="Putative2.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
f=length(PAdults[,1])
C1= for(f in 1:f) lapply(f,sorts)
unlink("Sort.txt")
unlink("IdnamesA.txt")
unlink("IdnamesO.txt")
#Begin exclude maternal allele==================================================#
Putative<- read.table("Putative2.txt", header=F, sep="\t", na.strings="?", dec=".", strip.white=TRUE)
Putative=Putative[,-2]
categories=ncol(Putative[,-1])
Aindivids=length(Putative[,1])
c=c(1:categories)
odd=2*(unique(round(((c-2))/2)))+1
l=length(odd) * 1000
codes=seq(from=1,to=l,by=1000)
cols=sort(rep(codes,2))-1
Anumbs=matrix(cols,Aindivids,categories,byrow=T)
Putative2=Putative[,-1]+Anumbs
Putative=cbind(Putative[,1],Putative2)
allelefinder=function(allelefinder)
{
names=cbind(as.character(Putative[z,1]),as.character(Putative[z+1,1]))
put1=Putative[z,-1]+100000                                             #if get error in this section (play with these numbers, and"reps" below
put2=Putative[z+1,-1]+100000
matches=match(put1,put2)
matches2=match(put2,put1)
finder=which(matches>0)
finder2=which(matches2>0)
allelef=cbind(put2[,1],put1[finder])
reps=table(as.numeric(signif(allelef[-1],3)))
remov=which(reps>1)
seeit=(which(reps==1))
alleles=allelef[-1]
alength=length(reps)
rdata=seq(from=1,to=alength,by=1)
dstack=cbind(reps,rdata)
dstackf=rep(dstack[,2],dstack[,1])
remfin=which(match(dstackf,remov)>0)
sees=which(match(dstackf,seeit)>0)
allhet=alleles[remfin]
allssing=alleles[sees]
alldup=sort(as.integer(rep(allssing,2)))
allallele=sort(as.integer(c(allhet,alldup)))
allallele=unlist(allallele)
put3=put2
detr=cbind(as.data.frame(sort(allallele)),sort(t(put3)))
dets=detr[,1]-detr[,2]
detf=which(dets==0)
detd=detr[-detf,2]
alla=allallele
dets=dets*100
da=which(dets>0)
db=which(dets<0)
dc=c(da,db)
dc=ceiling(dc/2)
dd=c((dc*2),((dc*2)-1))
newa=rep(detd,2)
aba=allallele[-dd]
alleles=c(aba,newa)
alleles=sort(alleles)
locn=((ncol(put1)))/2
locd=seq(from=100000, to=floor(max(alleles)), by=1000)
if(length(alleles)>0) {
removdiv=sort(rep(locd,2))
if(length(finder)==ncol(put1)) final=t(as.data.frame(rep(0,ncol(put1)))) else {
alla=as.data.frame(alleles)
aldiv=cbind(alla,removdiv)
final=aldiv[,1]-aldiv[,2]
final=t(final)}
final=cbind(names,final)
} else {final=cbind(names,Putative[z+1,-1]-Anumbs[1,])}
write.table(final,file="sharedalleles.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=length(Putative[,1])
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,allelefinder)
#===============================================================================#
#Start of Data Creation
Adults1<- DADS
loci=ncol(Adults1)
Adults=Adults1[,c(2:loci)]                                                      #assumes that there is one column of id names ; could modify this as needed.
Offspring1<- read.table("sharedalleles.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Offspring=Offspring1[,c(3:(loci+1))]
  total <- ncol(Adults)/2                                                          #For Progress bar
  if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = total, width = 300)}
afreqs=function(afreqs)                                                         #Begin Master simulation function
{
locus_name=L
wanted_reps=as.numeric(wanted_reps)
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, ceiling(locus_name/2), title=paste("Locus", ceiling(locus_name/2),"of ",total))}
vect=c(Adults[,L],Adults[,L+1])                                                 #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
alleles=data.frame(table(vect))
alleles=alleles[order(alleles[,1]),]
if (as.numeric(as.character(alleles[1,1]))<=0) {alleles=alleles[-1,]}
if(length(alleles[,1])==1) {                                                    #deals with monomorphic loci by adding 1 very strange allele (799)
alleles=cbind(vect[1],alleles[2])
alleles[2,]<-c(799,1)}
alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
homos=(alleles2[,3])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
hets=t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
csum=cumsum(as.numeric(gfreqs[,3]))
gfreqs1=cbind(gfreqs,csum)
Nadults=length(Adults[,1])
Noffs=length(Offspring[,1])
#end locus-specific HWE genotype frequency calculations  for adults
gfreqs2=cbind(Offspring[,L],Offspring[,L+1])
gf2=data.frame(table(paste(gfreqs2[,1],gfreqs2[,2])))                           #new to version 1.33 and up
gf2=cbind(gf2,gf2[,2]/sum(gf2[,2]))                                             #new to version 1.33 and up
gf3=(unlist(strsplit(as.character(gf2[,1]),"\\s")))                             #new to version 1.33 and up
gf4=cbind(gf3[seq(from=1,to=length(gf3),by=2)],gf3[seq(from=2,to=length(gf3),by=2)])
gfreqs2=cbind(gf4,gf2)                                                          #new to version 1.33 and up
gfreqs2=gfreqs2[,c(1,2,5)]                                                      #new to version 1.33 and up
alength=length(alleles2[,1])
for (Y in 1:wanted_reps) {
  positions=1:length(gfreqs[,1])
  sg3=sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                    #sample the repeated genotype positions, by the number of adults
  sadults=gfreqs[sg3,1:2]                                                       #index gfreqs to create genotypes
  positions=1:length(gfreqs2[,1])
  og3=sample(positions,Noffs,replace=TRUE,prob=gfreqs2[,3])                     #create juvenile genotyes
  soffs=gfreqs2[og3,1:2]
  soffs=cbind(as.numeric(as.character(soffs[,1])),as.numeric(as.character(soffs[,2])))
#Begin calcualtions of shared allele frequencies================================#
  #Create all pairwise compparisons between adults and offspring (Ads and Offs)
asp=cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
asp=rbind(cbind(locus_name,0,0),asp)
for (i in 1:Nadults) {
parent1=as.numeric(sadults[i,1])                                                #first allele in parent
parent2=as.numeric(sadults[i,2])                                                #second allele in parent
p1=soffs-parent1
p2=soffs-parent2
pp1=which(p1[,1]==0)
pp2=which(p1[,2]==0)
allele1=unique(c(pp1,pp2))
p21=which(p2[,1]==0)
p22=which(p2[,2]==0)
allele2=unique(c(p21,p22))
Out51=cbind(parent1,length(allele1))
Out52=cbind(parent2,length(allele2))
Out53=cbind(0,Noffs-length(unique(c(allele1,allele2))))
Out5=rbind(Out51,Out52,Out53)
if(parent2==parent1) {Out5=Out5[-1,]}                                           #remove 1 of alleles count if homozygous
if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
  diffs=sum(Out5[,2])-Noffs                                                     #remove to be more conservative!
  maxa=max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
  pos=which(Out5[,2]==maxa)
  Out5[pos,2]<-Out5[pos,2]-diffs}
m1=match(Out5[,1],asp[,2])
m2=asp[m1,3]+as.numeric(Out5[,2])
asp[m1,3]<-m2
asp<-asp
}
write.table(asp,file="out.sims",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
}
}
L=ncol(Adults)
C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
if(.Platform$OS.type=="windows"){close(pb)}
#Begin Sims_Read in=============================================================#
OUT<- read.table("out.sims", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
locname=unique(OUT[,1])                                                         #compile calculations for each locus
OUTALL=NULL
for (z in locname) {
  Loc1=OUT[which(OUT[,1]==z),]
  allfreqs=unique(Loc1[,2])
  OUT2=NULL
    for (x in allfreqs) {
    a1<-Loc1[which(Loc1[,2]==x),]
    a2=sum(a1[,3])
    a3=cbind(x,a2)
    OUT2 <- rbind(OUT2, a3)
  }
  OUT3=cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
  OUTALL <- rbind(OUTALL, cbind(z,OUT3))
}

NL=length(unique(OUTALL[,1]))
ngtypes=1                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
Ntotal=as.numeric(Ntotal)
inreps=50000     #was tested as 10000 for SNPS
repnumber=round(Ntotal/inreps)                                                #this is the rep number to get 100,000 values.  can adjust accordingly
asp=cbind(0:NL,rep(0,length(0:NL)))
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = Ntotal , width = 300)}
for (n in 1:repnumber) {
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, n*inreps, title=paste("Genotype", n*inreps,"of ",Ntotal))}
    OUT=NULL
    for (r in unique(OUTALL[,1])) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*Loco[,4])
        findo=which(alleles3[,2]==0)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
  distm=apply(OUT, 1, function(x)sum(x == 1))
  distm2=data.frame(table(distm))
  m1=match(distm2[,1],asp[,1])
  m2=asp[m1,2]+distm2[,2]
  asp[m1,2]<-m2

  asp<-asp
}
if(.Platform$OS.type=="windows"){close(pb)}
d2=asp
d3=cbind(d2,d2[,2]/sum(d2[,2]))
#Create plot of exclusionary power.  Also, some necessary data formatting========#
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "Posterior Processing", min = 0, max = 6, width = 300)}
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 1, title=paste("Performing Exclusion"))}
Adults<- DADS
Offs<- read.table("sharedalleles.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Nads=length(Adults[,1])
Noffs=length(Offs[,1])
asp=d3
asp=cbind(asp,NL-as.numeric(as.character(asp[,1])))  #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
asp=cbind(asp,cumsum(asp[,2]))
asp=cbind(asp,asp[,5]/max(asp[,5]))
distm=cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
mismatch=which.min(abs(distm[,6] - .5))                                         #could cause issues with a value of .5 chosen here - could be too low.  Change to higher if all loci analyzed have phi < 1.
Adults1=Adults                                                                  #begin exclusion
Offspring1=Offs
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Offspring=Offspring1[,c(3:(a+1))]                                               #special formatting for this script only!
colnames(Adults)<-""
colnames(Offspring)<-""
Anames=Adults1[,1]
Onames=Offspring1[,2]
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
A=1:Aindivids
O=1:Oindivids
G=expand.grid(A,O)
AG=G[,1]
AO=G[,2]
Ads=Adults[AG,]
Offs=Offspring[AO,]
IdnamesA=Anames[AG]
IdnamesO=Onames[AO]
write.table(IdnamesA,file="IdnamesA.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
write.table(IdnamesO,file="IdnamesO.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
IdnamesA<- read.table("IdnamesA.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
IdnamesO<- read.table("IdnamesO.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
Names=cbind(IdnamesA,IdnamesO)
matches=function(matches)
{
A=Ads[,z]-Offs[,z]
B=Ads[,(z+1)]-Offs[,(z+1)]
C=Ads[,z]-Offs[,(z+1)]
D=Ads[,(z+1)]-Offs[,z]
f=A*B*C*D
f=(f^2)*10
ss=which(is.na(f)==TRUE)
f=replace(f,ss,0)
identify=which(f>0)
f=replace(f,identify,1)
f=cbind(z,f)
write.table(f,file="Sort.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Ads)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("Sort.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
a=length(U[,1])
stuff=rowSums(U)
Sorted=cbind(Names,stuff)
matches=which(Sorted[,3]<(mismatch+1))
Actual=sort(stuff)
IDS=which(stuff<(mismatch+1))
PAdults=Ads[IDS,]
POffspring=Offs[IDS,]
nput=length(matches)
Putativepairs=Sorted[matches,]
PAdults=cbind(Putativepairs[,c(1,3)],PAdults)
POffspring=cbind(Putativepairs[,c(2,3)],POffspring)
names(PAdults)[1]<-"ID"
names(POffspring)[1]<-"ID"
sorts=function(sorts)
{
tell=rbind(PAdults[f,],POffspring[f,])
write.table(tell,file="Output_genotypes.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
f=length(PAdults[,1])
C1= for(f in 1:f) lapply(f,sorts)
unlink("Sort.txt")
unlink("IdnamesA.txt")
unlink("IdnamesO.txt")
#Calculate phi for each number mismatching loci=================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 2, title=paste("Calculating prior"))}
Putative<- read.table("Output_genotypes.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
observed=data.frame(table(Putative[,2]))
observed=cbind(observed,observed[,2]/2)                                         #done becuase each number is written twice in file (once for parent and once for off)
zerom=0:mismatch                                                                #this chunk adds 0s for mismatches where there were no observed putative pairs
zerom2=which(is.na(match(zerom,observed[,1])))
if (length(zerom2>0)) {observed=observed[,3]
for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}    #not really 0, to prevent divide by 0
}   else {observed=observed[,3]}
expected=distm[1:(mismatch+1),7]                                                #using cumulatinve sum   (more conservative)
#expected=distm[1:(mismatch+1),3]*Nads*Noffs                                    #not using cumulative sum
phi=expected/observed
phi=replace(phi,which(phi>=1),1)
Offs<- OFFSPRING
actualTrue=length(grep("Off",Offs[,1]))
info=cbind(actualTrue,expected,observed,phi)
#begin caluclation of distribution for lamda|phi ===============================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 3, title=paste("Creating distriubtion of falsely-shared alleles"))}
phibase=phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
observed=observed[min(which(phi[]==1))-1]                                       #do the same for observed
testob=which(observed==0.000000001)
phi2=cbind(1:length(phi),phi)
if (length(testob>0)) {phi4=phi2[-testob,]} else {phi4=phi2}
nmismatch=min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
index=phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
index=index[which(index>-1)]
if (length(index)>1) {
  if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
phi=phi[index]
index=index-1
#Create Plot ===================================================================#
pdf(file="Output_Dataset_Power.pdf")
x=0:(length(info[,1])-1)
y1=info[,3]
y2=info[,2]
p1=which(y2==0)
y2=replace(y2,p1,.000000001)
p1=which(y1<y2)
y3=replace(y1,p1,y2[p1])
y2=y2+1
y3=y3+1
par(mar=c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
ats=c(0,1,2,3,4,5,6,7,8,9)
labs=c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
axis(side=2,ats,labs)
lines(x,log10(y3),lwd=2)
lines(x,log10(y2),lwd=2)
points(x,log10(y3),pch=21,bg="green",cex=2)
points(x,log10(y2),pch=21,bg="blue",cex=2)
legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
yphi=y2/y3
par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
points(x,yphi,cex=2,pch=21,bg="gray")
mtext("Number of Mismatching Loci",side=1,line=1.94)
dev.off()
info2=cbind(x,info[,4])
colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
#Begin calculation of distribution for lamda|phi ===============================#
ngtypes=20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
inreps=10000
repnumber=round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
for (n in 1:repnumber) {
     OUT=NULL
    for (r in unique(OUTALL[,1])) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*Loco[,4])
        findo=which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        alleles3=alleles3[-which(alleles3[,2]==0),]
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
for (i in index) {                                                              #loop over numbers of mismatched loci
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                    #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distp=cbind(i,DIST)
    write.table(distp,file="False_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
 }
}
#Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
OUT9=NULL
for (n in index){
  Putative2=Putative[which(Putative[,2]==n),]
  write.table(Putative2, file="Putative.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    #if (length(distp[,1])>1000) {                                               #need at least 1000 values , else assigned 0 (may be redundant now)
      OBS=NULL
      afreqs=function(afreqs) {
        PutL=Putative2[,c(L+2,L+3)]
        PutLadults=PutL[seq(from=1,to=length(PutL[,1]),by=2),]                  #Combine putative pairs alleles into a single row
        PutLoffs=PutL[seq(from=2,to=length(PutL[,1]),by=2),]
        Puts=cbind(PutLadults,PutLoffs)
        c1=Puts[,1]-Puts[,3]                                                    #find the matching alleles
        c2=Puts[,1]-Puts[,4]
        c3=Puts[,2]-Puts[,3]
        c4=Puts[,2]-Puts[,4]
        Puts2=cbind(Puts,c1,c2,c3,c4)
        P5=replace(Puts2[,5],which(Puts2[,5]!=0),-10)
        P6=replace(Puts2[,6],which(Puts2[,6]!=0),-10)
        P7=replace(Puts2[,7],which(Puts2[,7]!=0),-10)
        P8=replace(Puts2[,8],which(Puts2[,8]!=0),-10)
        P5=replace(P5,which(P5==0),1)
        P6=replace(P6,which(P6==0),1)
        P7=replace(P7,which(P7==0),1)
        P8=replace(P8,which(P8==0),1)
        Puts3=cbind(Puts,P5,P6,P7,P8)
        Puts4=cbind((Puts3[,1]*Puts3[,5]),(Puts3[,1]*Puts3[,6]),(Puts3[,2]*Puts3[,7]),(Puts3[,2]*Puts3[,8]))
        alleles2=OUTALL[which(OUTALL[,1]==L),]
        alfreq1=alleles2[match(Puts4[,1],alleles2[,2]),4]
        alfreq2=alleles2[match(Puts4[,2],alleles2[,2]),4]                       #find the actual allele values
        alfreq3=alleles2[match(Puts4[,3],alleles2[,2]),4]
        alfreq4=alleles2[match(Puts4[,4],alleles2[,2]),4]
        Puts5=cbind(alfreq1,alfreq2,alfreq3,alfreq4)                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
        R1=replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
        R2=replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
        R3=replace(Puts5[,3],which(is.na(Puts5[,3])==TRUE),1)
        R4=replace(Puts5[,4],which(is.na(Puts5[,4])==TRUE),1)
        Puts6=cbind(R1,R2,R3,R4)
        Put_share=apply(Puts6, 1, min)                                          #find row minimum
        Put_share2=apply(Puts4, 1, max)                                         #find shared allele name
        Put_share3=c(Put_share,Put_share2)
        OBS <<- cbind(OBS,Put_share3)
    }
    L=ncol(Putative2)-2
    C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    lengths=length(OBS[,1])/2
    if (lengths==1) {OBA3=t(OBS[(lengths+1):(2*lengths),])} else {OBA3=OBS[(lengths+1):(2*lengths),]}
    write.table(OBA3,file="True_shared_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)  #Actual shared alleles   #This file wouuld be useful as output for people to identify shared and mismatching loci
    if (lengths==1) {OBS=t(OBS[1:lengths,])} else {OBS=OBS[1:lengths,]}         #formatting for if there is only a single pair
    obsp=apply(OBS, 1, prod)
    #}  else obsp=rep(0,(length(Putative2[,1]))/2)
OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
}
#calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic)==#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 4, title=paste("Calculating Posterior Component 1"))}
OBA3<- read.table("True_shared_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
for (n in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
      OUT=NULL
      for (r in unique(OUTALL[,1])) {                                           #This first section calculates products of parental alleles (True distribution)
        vect=c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
        alleles=data.frame(table(vect))
        alleles=alleles[order(alleles[,1]),]
        if (as.numeric(as.character(alleles[1,1]))<=0) {alleles=alleles[-1,]}
        alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))
        gtrue=sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
        OUT <- cbind(OUT,gtrue)

        if(n==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele
                            mm=alleles2[match(OBA3[i,ceiling(r/2)],alleles2[,1]),3]
                            write.table(cbind(r,mm),file="Shared_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
                            }
                }
        }
    for (i in index) {
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                   #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distt<-cbind(i,DIST)
    write.table(distt,file="True_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
}
#Calculate lamdaphi=============================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 5, title=paste("Calculating Posterior Component 2"))}
Putative3<- read.table("Putative.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Putadults=Putative3[seq(from=1,to=length(Putative3[,1]),by=2),1]                #Combine putative pairs alleles into a single row
Putoffs=Putative3[seq(from=2,to=length(Putative3[,1]),by=2),1]
Names=cbind(as.character(Putadults),as.character(Putoffs))
empirical=cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
distp <- read.table("False_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
P2=NULL
for (i in index) {                                                              #loop over numbers of mismatched loci
  empirical2=empirical[which(as.numeric(as.character(empirical[,3]))==i),]
  if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
  if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
    a3=distp[which(distp[,1]==i),2]
    DIST<-as.data.frame(a3)
    P=NULL
    for (b in 1:length(empirical2[,1])) {
      p1=length(which(DIST[,1] <= as.numeric(empirical2[b,4]) ))
      if (p1==0) {p1=0.00001}
      p2=cbind(empirical2[b,1],empirical2[b,2],p1)
      p3=cbind(p2,as.numeric(p2[,3])/length(DIST[,1]))
      P <- rbind(P,p3)
      }
    }
  P2<-rbind(P2,cbind(i,P))
}
lamdaphi=as.numeric(as.character(P2[,5]))
#Calculate lamda|phic ==========================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 6, title=paste("Calculating Posterior Component 3 and Finishing Calculations"))}
lamdaphic_dist<- read.table("True_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Observed<- read.table("Shared_allele_freqs.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
mm1=which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
Observed[mm1,2] <- 1                                                            #replace NAs with 1
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
lamdaphic=apply(U, 1, prod)
l1=length(which(OUT9[,2]==0))
if (l1>0) lamdaphic=c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)    #double check values by hand
P3=cbind(P2,lamdaphic)
for (i in index) {                                                              #loop over numbers of mismatched loci
   e2=P3[which(as.numeric(as.character(P3[,1]))==i),]
   if(length(e2)==6) {e2=t(e2)}                                                 #deals with one indiviudal formatting
   if (e2[1,5]==0) {write.table(e2[,6], file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{                    #deals with not enough reps
  a3=lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
  DIST<-as.data.frame(a3)
  for (b in 1:length(e2[,1])) {                                                 #calculate p values
    p1=length(which(DIST[,1] <= e2[b,6]))
    p2=p1/length(DIST[,1])
    write.table(p2, file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    }
  }
}
lamdaphic<- read.table("lamdaphic.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
#Put it all together with Bayes theorem! =======================================#
vals=cbind(P2,lamdaphic[,1])
philength=cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
phis=rep(philength[,2],philength[,3])
vals=cbind(vals,phis)
colnames(vals)<-c("Nmismatch","Parent","Off","ignore","lamdaphi","lamdaphic","phi")
phi=as.numeric(as.character(vals[,7]))
lamdaphi=as.numeric(as.character(vals[,5]))
lamdaphic=as.numeric(as.character(vals[,6]))
lamdaphi=replace(lamdaphi,which(lamdaphi==0),1)
lamdaphic=replace(lamdaphic,which(lamdaphic==0),1)
pval=(lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
pval=cbind(vals[,2],vals[,3],vals[,1],pval)
pval=pval[order(as.numeric(pval[,4])),]
colnames(pval) <- c("Adult","Offspring","NL_mismatch","Pvalue")
#write.table(vals, file="Bayesian_Output.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
if(.Platform$OS.type=="windows"){close(pb)}
unlink("False_allele_freqs.txt")                                                #clear all sims files
unlink("True_allele_freqs.txt")
unlink("Shared_allele_freqs.txt")
unlink("lamdaphic.txt")
unlink("Putative.txt")
unlink("True_shared_freqs.txt")
#unlink("Output_genotypes.txt")
unlink("*.sims")
unlink("sharedalleles.txt")
unlink("Putative2.txt")
rm(list=ls())
#gc()
}
label.bayes <- tklabel(tt, text="Press 'Run' to perform Bayesian parentage analysis:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label.bayes,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label.bayes, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
label.BAYES1 <- tklabel(tt, text="Perform Bayesian parentage analysis with 1 known parent:",font=fontTextLabel)
OK.but.BAYES1 <- tkbutton(tt, text = "BAYES", command = PressedBAYES1)
tkgrid(tklabel(tt,text="     "),label.BAYES1,tklabel(tt,text="     "),OK.but.BAYES1,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.BAYES1, sticky="w")
tkgrid.configure(OK.but.BAYES1, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkfocus(tt)
heading <- tklabel(tt, text="Parentage with Known Parent-Pairs",font=fontSUB)   # add a heading
tkgrid(heading, columnspan=5)
#Module9_exclusion_knownpairs############################################################################################################################################################################################################################
#Exclusion_KnownParentPair_Module_8.27
PressedEXCLUSION2 <- function()
{
#Create and name toplevel window
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Parentage Analysis with Known Parent-Pairs "                          # Name the window
heading <- tklabel(tt, text="Exclusion",font=fontHeading)                       # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set working Directory =========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()
{
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Load Known Adults File ========================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
 # if (missing == "yes") {missing1=as.integer(0)} else {missing1=" "}      #encode missing option
  Adults <- read.table(fileName, header=T, sep="\t", na.strings=-1, dec=".", strip.white=TRUE)
  assign("MOMS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Mother Genotype File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing mother genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Putative Adults File ======================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
#  if (missing == "yes") {missing1=as.integer(0)} else {missing1=" "}      #encode missing option
  Adults <- read.table(fileName, header=T, sep="\t", na.strings=-1, dec=".", strip.white=TRUE)
  assign("DADS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select Father Genotype File", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing father genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Offspring File============================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
 # if (missing == "yes") {missing1=as.integer(0)} else {missing1=" "}      #encode missing option
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings=-1, dec=".", strip.white=TRUE)
  assign("OFFSPRING", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
#Create button to run parentage script =========================================#
PressedOK <- function()
    {
Adults1<- MOMS
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Numberloci=ncol(Adults)
Offspring1<- DADS
Offspring=Offspring1[,c(2:a)]
Juvs<- OFFSPRING
Juv=Juvs[,c(2:a)]
#Begin combining of alleles between broodstock pairs (grandparents)
Anames=Adults1[,1]
Onames=Offspring1[,1]
names2=cbind(as.character(Anames),as.character(Onames))
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
matches=function(matches)
{
A=Adults[,z]
B=Adults[,(z+1)]
C=Offspring[,(z+1)]
D=Offspring[,z]
bale=cbind(z,A,B,C,D)
write.table(bale,file="pooled.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Adults)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("pooled.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=data.frame(rep(1,length(Adults1[,1])))
for (i in a) {
u=Observed[Observed[,1]==i,-1]
U=cbind(U,u)
}
gall=cbind(names2,U[,-1])
erase=""
write.table(erase,file="pooled.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
Jnames=Juvs[,1]
Juvs2=Juvs[,-1]
gnames=gall[,1:2]
gall2=gall[,-c(1,2)]
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = length(Juvs[,1]) , width = 300)}
OUT3=NULL
for  (i in 1:length(Juvs[,1])){
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, i, title=paste("Exclusion comparison", i,"of ",length(Juvs[,1])))}
j1=Juvs2[i,]
  OUT=data.frame(rep(Jnames[i],length(gall2[,1])))
  for (l in seq(from=1,to=ncol(j1),by=2)) {
    j2=j1[,c(l,l+1)]                                 #isolate gtypes by locus
    #gindex=seq(from=1,to=ncol(gall2),by=4)
    start1=(l*2)-1
    stop1=start1+3
    gall3=gall2[,start1:stop1]
    j3=as.numeric(j2[1])
    j4=as.numeric(j2[2])
    #first allele
    #match to mom?
    m11=abs(gall3[,1]-j3)
    m12=abs(gall3[,2]-j3)
    #match to dad?
    d11=abs(gall3[,3]-j3)
    d12=abs(gall3[,4]-j3)
    #second allele
    #match to mom?
    m21=abs(gall3[,1]-j4)
    m22=abs(gall3[,2]-j4)
    #match to dad?
    d21=abs(gall3[,3]-j4)
    d22=abs(gall3[,4]-j4)
    cbind(m11,m12,d11,d12,m21,m22,d21,d22)
    #tests
    t1=which((m11+d21)==0)
    t2=which((m11+d22)==0)
    t3=which((m12+d21)==0)
    t4=which((m12+d22)==0)
    t5=which((m21+d11)==0)
    t6=which((m22+d11)==0)
    t7=which((m21+d12)==0)
    t8=which((m22+d12)==0)
    results=sort(unique(c(t1,t2,t3,t4,t5,t6,t7,t8)))
    #data formatting
    dt1=1:length(gall3[,1])
    m1=match(dt1,results)
    r1=which(is.na(m1))
    m1[r1]<-0
    m1[-r1]<-1
    if((length(r1)>0)==FALSE) {m1=rep(1,length(m1))}   #deals with all match
    out=data.frame(m1)
    OUT=cbind(OUT,out)
    }
OUT2=cbind(gnames,OUT)
colnames(OUT2)<-c("Mothers","Fathers","Offspring")
OUT3=rbind(OUT3,OUT2)
}
write.table(OUT3,file="AllComparisons_Output.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=T)
if(.Platform$OS.type=="windows"){close(pb)}
unlink("pooled.txt")
}
label_exclusion <- tklabel(tt, text="Press 'Run' to perform exclusion:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label_exclusion,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label_exclusion, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
#Set number of mismatching loci ================================================#
label_mismatches <- tklabel(tt, text="Set number of loci to let mismatch:",font=fontTextLabel)
label_mismatches2 <- tklabel(tt, text="Tip: can use repeatedly after performing exclusion")
Name1 <- tclVar("0")
entry.label2 <- tkentry(tt, width="4",textvariable=Name1)                       #create entry fields
OnOK2 <- function()
{
	NameVal <- tclvalue(Name1)
	assign("mismatch", NameVal, envir = solomon.env)
	Alldata<- read.table("AllComparisons_Output.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) #make sure has header
  Sorts=cbind(Alldata,rowSums(Alldata[,4:ncol(Alldata)]))
  mismatch=as.numeric(mismatch)
  NL2=(ncol(OFFSPRING)-1)/2
  mismatch2=NL2-mismatch
  good=Sorts[(which(Sorts[,ncol(Sorts)]>=mismatch2)),]                               #set stringency of parent matches here
  good=good[order(good[,1]),]
  colsnmz=rep("Locus",ncol(good)-4)
  colsnmz2=sort(1:(length(colsnmz)))
  colsnmz3=paste(colsnmz,colsnmz2)
  colsnmz4=c("Mother","Father","Offspring",colsnmz3,"Loci matching")
  colnames(good)<-colsnmz4
  write.table(good,file="Exclusion_Output.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
  Females=data.frame(table(good[,2]))
  Females=Females[which(Females[,2]>0),]
  Females=Females[order(Females[,2]),]
  colnames(Females)<-c("ID","Frequency")
  males=data.frame(table(good[,1]))
  males=males[which(males[,2]>0),]
  males=males[order(males[,2]),]
  colnames(males)<-c("ID","Frequency")
  write.table(males,file="Mother_Assignments_Output.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
  write.table(Females,file="Father_Assignments_Output.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)
}
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label2, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_mismatches, tklabel(tt,text="     "), entry.label2, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "), label_mismatches2)
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_mismatches,label_mismatches2, sticky="w")
tkgrid.configure(entry.label2,sticky="w")
}
label.EXCLUSION2 <- tklabel(tt, text="Perform exclusion with known parent-pairs:",font=fontTextLabel)
OK.but.EXCLUSION2 <- tkbutton(tt, text = "EXCLUSION", command = PressedEXCLUSION2)
tkgrid(tklabel(tt,text="     "),label.EXCLUSION2,tklabel(tt,text="     "),OK.but.EXCLUSION2,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.EXCLUSION2, sticky="w")
tkgrid.configure(OK.but.EXCLUSION2, sticky="w")
tkfocus(tt)
#Module10########################################################################################################################################################################################################################
#BayesianParentage_KnownParentPair_Module_9.6.12
PressedBAYES2 <- function()
{

#Create and name toplevel window
fontHeading <- tkfont.create(family="times",size=18,weight="bold")
fontTextLabel <- tkfont.create(family="times",size=14)
tt <- tktoplevel()                                                              # Create a new toplevel window; Note this window is called tt (could create other windows with different names)
tktitle(tt) <- "SOLOMON: Bayesian Parentage Analysis with Known Parent Pairs  "                          # Name the window
heading <- tklabel(tt, text="SOLOMON: Bayesian Parentage Analysis",font=fontHeading)        # add a heading
tkgrid(heading, columnspan=5)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
#Set working Directory =========================================================#
label_working.directory <- tklabel(tt, text="Set working directory (use forward slash):",font=fontTextLabel)
Name <- tclVar("C:/SOLOMON")
entry.label1 <- tkentry(tt, width="20",textvariable=Name)                       #create entry fields
OnOK <- function()
{
	NameVal <- tclvalue(Name)
	msg <- paste("You have now set the working directory to",NameVal)
	tkmessageBox(message=msg)
	assign("directory", NameVal, envir = solomon.env)
	setwd(NameVal)
}
OK.but <-tkbutton(tt,text="   OK   ",command=OnOK)
tkbind(entry.label1, "<Return>",OnOK)
tkgrid(tklabel(tt,text="     "), label_working.directory, tklabel(tt,text="     "), entry.label1, tklabel(tt,text="     "), OK.but, tklabel(tt,text="     "))
tkfocus(tt)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_working.directory, sticky="w")
tkgrid.configure(entry.label1, sticky="w")
#Set Number of Sims=============================================================#
label_sim.number <- tklabel(tt, text="Number of simulated data sets:",font=fontTextLabel)
Name.sim <- tclVar(1000)
entry.label.sim <- tkentry(tt, width="10",textvariable=Name.sim)                #create entry fields
OnOK2 <- function()  {
	NameVal <- tclvalue(Name.sim)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("wanted_reps", NameVal, envir = solomon.env)
  }
OK.but2 <-tkbutton(tt,text="   OK   ",command=OnOK2)
tkbind(entry.label1, "<Return>",OnOK2)
tkgrid(tklabel(tt,text="     "), label_sim.number, tklabel(tt,text="     "), entry.label.sim, tklabel(tt,text="     "), OK.but2, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special1=tklabel(tt, text="Recommended: Microsatellites=1000 SNPs=100",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special1, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_sim.number, sticky="w")
tkgrid.configure(entry.label.sim,label_Special1, sticky="w")
#Set Number of "Geonotypes"=====================================================#
label_Ntotal <- tklabel(tt, text="Number of simulated genotypes:",font=fontTextLabel)
Name.Ntotal <- tclVar(50000000)
entry.label.Ntotal <- tkentry(tt, width="10",textvariable=Name.Ntotal)                #create entry fields
OnOKN <- function()  {
	NameVal <- tclvalue(Name.Ntotal)
	msg <- paste("You have now set the number of simulated data sets to",NameVal)
	tkmessageBox(message=msg)
	assign("Ntotal", NameVal, envir = solomon.env)
  }
OK.butN <-tkbutton(tt,text="   OK   ",command=OnOKN)
tkbind(entry.label1, "<Return>",OnOKN)
tkgrid(tklabel(tt,text="     "), label_Ntotal, tklabel(tt,text="     "), entry.label.Ntotal, tklabel(tt,text="     "), OK.butN, tklabel(tt,text="     "))
SpecialFont <- tkfont.create(family="times",size=11)
label_Special=tklabel(tt, text="Recommended: Microsatellites=50,000,000 SNPs=500,000",font=SpecialFont)
tkgrid(tklabel(tt,text="     "), label_Special, tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(label_Ntotal, sticky="w")
tkgrid.configure(entry.label.Ntotal,label_Special, sticky="w")
#Load Known Adults File  =======================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("MOMS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select File Containing Mothers", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing mother genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Putative Adults File =====================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
   Adults <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("DADS", Adults, envir = solomon.env)
}
adults.button <- tkbutton(tt, text = "Select File Containg Fathers", command = getfile)
adults_label <- tklabel(tt, text="Please select file containing father genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),adults_label,tklabel(tt,text="     "),adults.button)
tkgrid(tklabel(tt,text="     "))                                                #add a blank line
tkgrid.configure(adults_label, sticky="w")
tkgrid.configure(adults.button,sticky="w")
#Load Offspring File ===========================================================#
getfile <- function(){
  fileName <- tclvalue(tkgetOpenFile())
  if (!nchar(fileName)) {
  tkmessageBox(message = "No file was selected!")
  } else {
  tkmessageBox(message = paste("The file selected was", fileName))
  }
  Offspring <- read.table(fileName, header=T, sep="\t", na.strings="-1", dec=".", strip.white=TRUE)
  assign("OFFSPRING", Offspring, envir = solomon.env)
}
offspring.button <- tkbutton(tt, text = "Select Offspring File", command = getfile)
offspring_label <- tklabel(tt, text="Please select file containing offspring genotypes:",font=fontTextLabel)
tkgrid(tklabel(tt,text="     "),offspring_label,tklabel(tt,text="     "),offspring.button)
tkgrid(tklabel(tt,text="     "))
tkgrid.configure(offspring_label, sticky="w")
tkgrid.configure(offspring.button, sticky="w")
#Create button to run parentage script =========================================#
PressedOK <- function()
    {
reps=1                                                                          #see end of script
Adults2<- OFFSPRING                                                             #takes all 3 files, finds unique, and calculates freqs
Adults3<- DADS
Adults4<- MOMS
Adults5=rbind(Adults2,Adults3,Adults4)
Adults1=Adults5[match(unique(Adults5[,1]),Adults5[,1]),]
loci=ncol(Adults1)
Adults=Adults1[,c(2:loci)]                                                      #assumes that there is one column of id names ; could modify this as needed.
Offspring1<- OFFSPRING
Offspring=Offspring1[,c(2:loci)]
  total <- ncol(Adults)/2                                                       #For Progress bar
  if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = total, width = 300)}
afreqs=function(afreqs)                                                         #Begin Master simulation function
{
locus_name=L
wanted_reps=as.numeric(wanted_reps)
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, ceiling(locus_name/2), title=paste("Locus", ceiling(locus_name/2),"of ",total))}
vect=c(Adults[,L],Adults[,L+1])                                                 #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
alleles=data.frame(table(vect))
alleles=alleles[order(alleles[,1]),]
if (as.numeric(as.character(alleles[1,1]))<=0) {alleles=alleles[-1,]}
alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))                            #table allele frequencies
homos=(alleles2[,3])^2                                                          #create homozygote allele frequencies
homos2=cbind(as.character(alleles2[,1]),as.character(alleles2[,1]),homos)
hets=t(combn(alleles2[,3],2))                                                   #create heterozygote allele frequencies
hetfreq=2*(hets[,1]*hets[,2])
hetvals=t(combn(as.character(alleles2[,1]),2))                                  #create heterozygote allele names
hets2=cbind(hetvals,hetfreq)
gfreqs=rbind(hets2,homos2)                                                      #combine hets and homos and create genotypes
csum=cumsum(as.numeric(gfreqs[,3]))
gfreqs1=cbind(gfreqs,csum)
Nadults=length(Adults[,1])
Noffs=length(Offspring[,1])
#end locus-specific HWE genotype frequency calculations
wanted_number_of_reps=reps
alength=length(alleles2[,1])
for (Y in 1:wanted_reps) {
  positions=1:length(gfreqs[,1])
  sg3=sample(positions,Nadults,replace=TRUE,prob=gfreqs[,3])                    #sample the repeated genotype positions, by the number of adults
  sadults=gfreqs[sg3,1:2]                                                       #index gfreqs to create genotypes
  og3=sample(positions,Noffs,replace=TRUE,prob=gfreqs[,3])                      #create juvenile genotyes
  soffs=gfreqs[og3,1:2]
  soffs=cbind(as.numeric(as.character(soffs[,1])),as.numeric(as.character(soffs[,2])))
  asp=cbind(rep(locus_name,alength),as.numeric(as.character(alleles2[,1])),rep(0,alength))
  asp=rbind(cbind(locus_name,0,0),asp)
  for (i in 1:Nadults) {
    parent1=as.numeric(sadults[i,1])                                                #first allele in parent
    parent2=as.numeric(sadults[i,2])                                                #second allele in parent
    p1=soffs-parent1
    p2=soffs-parent2
    pp1=which(p1[,1]==0)
    pp2=which(p1[,2]==0)
    allele1=unique(c(pp1,pp2))
    p21=which(p2[,1]==0)
    p22=which(p2[,2]==0)
    allele2=unique(c(p21,p22))
    Out51=cbind(parent1,length(allele1))
    Out52=cbind(parent2,length(allele2))
    Out53=cbind(0,Noffs-length(unique(c(allele1,allele2))))
    Out5=rbind(Out51,Out52,Out53)
    if(parent2==parent1) {Out5=Out5[-1,]}                                           #remove 1 of alleles count if homozygous
    if(sum(Out5[,2])>Noffs) {                                                       #remove most common allele for double heterozygoutes
      diffs=sum(Out5[,2])-Noffs                                                     #remove to be more conservative!
      maxa=max(c(Out51[,2],Out52[,2]))                                              #will be removed twice if have exact same allele count!
      pos=which(Out5[,2]==maxa)
      Out5[pos,2]<-Out5[pos,2]-diffs}
      m1=match(Out5[,1],asp[,2])
      m2=asp[m1,3]+as.numeric(Out5[,2])
      asp[m1,3]<-m2
      asp<-asp
    }
  write.table(asp,file="out.sims",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
  }
}
L=ncol(Adults)
C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
if(.Platform$OS.type=="windows"){close(pb)}
#begin Sims_Read in ============================================================#
OUT<- read.table("out.sims", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
locname=unique(OUT[,1])                                                         #compile calculations for each locus
OUTALL=NULL
for (z in locname) {
  Loc1=OUT[which(OUT[,1]==z),]
  allfreqs=unique(Loc1[,2])
  OUT2=NULL
    for (x in allfreqs) {
    a1<-Loc1[which(Loc1[,2]==x),]
    a2=sum(a1[,3])
    a3=cbind(x,a2)
    OUT2 <- rbind(OUT2, a3)
  }
  OUT3=cbind(OUT2,OUT2[,2]/sum(OUT2[,2]))
  OUTALL <- rbind(OUTALL, cbind(z,OUT3))
}
NL=length(unique(OUTALL[,1]))
ngtypes=1                                                                       #20 million seems like plenty for all datasets (but may need to adjust at some point) (deprecated)
Ntotal=as.numeric(Ntotal)
inreps=50000                                                                    #was tested as 10000 for SNPS
repnumber=round(Ntotal/inreps)                                                  #this is the rep number to get 100,000 values.  can adjust accordingly
asp=cbind(0:(NL*2),rep(0,length(0:(NL*2))))
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = Ntotal , width = 300)}
for (n in 1:repnumber) {
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, n*inreps, title=paste("Genotype", n*inreps,"of ",Ntotal))}
    OUT=NULL
    for (r in (sort(rep(unique(OUTALL[,1]),2)))) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*(as.numeric(as.character(Loco[,4]))))
        findo=which(alleles3[,2]==0)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
  distm=apply(OUT, 1, function(x)sum(x == 1))
  distm2=data.frame(table(distm))
  m1=match(distm2[,1],asp[,1])
  m2=asp[m1,2]+distm2[,2]
  asp[m1,2]<-m2
  asp<-asp
 }
if(.Platform$OS.type=="windows"){close(pb)}
d2=asp
d3=cbind(d2,d2[,2]/sum(d2[,2]))
#Create plot of exclusionary power.  Also, some necessary data formatting ======#
Adults<- DADS
Offs<- OFFSPRING
Nads=length(Adults[,1])
Noffs=length(Offs[,1])
NL2=NL*2
asp=d3
asp=cbind(asp,NL2-as.numeric(as.character(asp[,1])))                            #first column represents the number of loci that mismatch, thus reverse sorting equals number of loci that match
asp=cbind(asp,cumsum(asp[,2]))
asp=cbind(asp,asp[,5]/max(asp[,5]))
distm=cbind(asp,asp[,6]*Nads*Noffs)                                             #calc Nloci to let mismatch
#Begin exlusion to parent-pairs ================================================#
Adults1<- MOMS
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Numberloci=ncol(Adults)
Offspring1<- DADS
Offspring=Offspring1[,c(2:a)]
Juvs<- OFFSPRING
Juv=Juvs[,c(2:a)]
#Begin combining of alleles between broodstock pairs (grandparents)
Anames=Adults1[,1]
Onames=Offspring1[,1]
names2=cbind(as.character(Anames),as.character(Onames))
categories=ncol(Adults)
Aindivids=length(Adults[,1])
Oindivids=length(Offspring[,1])
matches=function(matches)
{
A=Adults[,z]
B=Adults[,(z+1)]
C=Offspring[,(z+1)]
D=Offspring[,z]
bale=cbind(z,A,B,C,D)
write.table(bale,file="pooled.txt",row.names=FALSE,col.names=F,sep="\t",append=T)
}
z=ncol(Adults)
C1= for(z in (2*(unique(round((1:(z-2))/2)))+1)) lapply(z,matches)
Observed<- read.table("pooled.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
a=unique(Observed[,1])
U=data.frame(rep(1,length(Adults1[,1])))
for (i in a) {
u=Observed[Observed[,1]==i,-1]
U=cbind(U,u)
}
gall=cbind(names2,U[,-1])
erase=""
write.table(erase,file="pooled.txt",row.names=FALSE,col.names=F,sep="\t",append=F)
Jnames=Juvs[,1]
Juvs2=Juvs[,-1]
gnames=gall[,1:2]
gall2=gall[,-c(1,2)]
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "progress bar", min = 0, max = length(Juvs[,1]) , width = 300)}
OUT3=NULL
for  (i in 1:length(Juvs[,1])){
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, i, title=paste("Exclusion comparison", i,"of ",length(Juvs[,1])))}
j1=Juvs2[i,]
  OUT=data.frame(rep(Jnames[i],length(gall2[,1])))
  for (l in seq(from=1,to=ncol(j1),by=2)) {
    j2=j1[,c(l,l+1)]                                 #isolate gtypes by locus
    #gindex=seq(from=1,to=ncol(gall2),by=4)
    start1=(l*2)-1
    stop1=start1+3
    gall3=gall2[,start1:stop1]

    j3=as.numeric(j2[1])
    j4=as.numeric(j2[2])
    #first allele
    #match to mom?
    m11=abs(gall3[,1]-j3)
    m12=abs(gall3[,2]-j3)
    #match to dad?
    d11=abs(gall3[,3]-j3)
    d12=abs(gall3[,4]-j3)
    #second allele
    #match to mom?
    m21=abs(gall3[,1]-j4)

    m22=abs(gall3[,2]-j4)
    #match to dad?
    d21=abs(gall3[,3]-j4)
    d22=abs(gall3[,4]-j4)
    #tests
    t1=which((m11+d21)==0)
    t2=which((m11+d22)==0)
    t3=which((m12+d21)==0)
    t4=which((m12+d22)==0)
    t5=which((m21+d11)==0)
    t6=which((m22+d11)==0)
    t7=which((m21+d12)==0)
    t8=which((m22+d12)==0)
    results=sort(unique(c(t1,t2,t3,t4,t5,t6,t7,t8)))
    #data formatting
    dt1=1:length(gall3[,1])
    m1=match(dt1,results)
    r1=which(is.na(m1))
    m1[r1]<-0
    m1[-r1]<-1
    if((length(r1)>0)==FALSE) {m1=rep(1,length(m1))}   #deals with all match
      asp=cbind(m11,m12,d11,d12,m21,m22,d21,d22)        #adds 0.5s
      asp2=apply(asp,1,prod)
      asp3=which(asp2==0)
      asp4=which(m1==1)
      a1=match(asp4,asp3)
      asp3=asp3[-a1]
    m1[asp3]<-0.5
    out=data.frame(m1)
    OUT=cbind(OUT,out)
    }
OUT2=cbind(gnames,OUT)
colnames(OUT2)<-c("Mothers","Fathers","Offspring")
#OUT3=rbind(OUT3,OUT2)
write.table(OUT2,file="AllComparisons_Output.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=T)
}
if(.Platform$OS.type=="windows"){close(pb)}
unlink("pooled.txt")
Alldata<- read.table("AllComparisons_Output.txt", header=F, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) #make sure has header
Sorts=cbind(Alldata,rowSums(Alldata[,4:ncol(Alldata)]))
#Calculate phi for each number mismatching loci ================================#
if(.Platform$OS.type=="windows"){pb <- winProgressBar(title = "Posterior Processing", min = 0, max = 6 , width = 300)}
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 1, title=paste("Posterior Processing."))}
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 2, title=paste("Calculating prior"))}
Putative<- Sorts
obss1=(Putative[,ncol(Putative)])*2
observed=data.frame(table(obss1))
zerom=seq(0,NL2,1)                                                              #this chunk adds 0s for mismatches where there were no observed putative pairs
zerom2=which(is.na(match(zerom,observed[,1])))
if (length(zerom2>0)) {observed=observed[,2]
for(i in 1:length(zerom2))  {observed <- append(observed, 0.000000001, after=(zerom2[i]-1))}    #not really 0, to prevent divide by 0
}   else {observed=observed[,2]}
observed=rev(observed)
expected=distm[,7]                                                              #using cumulatinve sum   (more conservative)
#expected=distm[,3]*Nads*Noffs                                                  #not using cumulative sum
phi=expected/observed
phi=replace(phi,which(phi>=1),1)
Offs<- OFFSPRING
actualTrue=length(grep("Off",Offs[,1]))
info=cbind(actualTrue,expected,observed,phi)
#calculate phi and index========================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 3, title=paste("Creating distriubtion of falsely-shared alleles"))}
#phibase=phi[min(which(phi[]==1))-1]                                             #remove all phis after first 1 is observed (conservative)
#observed=observed[min(which(phi[]==1))-1]                                       #do the same for observed
testob=which(observed==0.000000001)
phi2=cbind(1:length(phi),phi)
if (length(testob>0)) {phi4=phi2[-testob,]} else {phi4=phi2}
nmismatch=min(which(phi4[,2]==1))-1                                             #takes loci before the first 1
index=phi4[1:nmismatch,1]                                                       #only perform analyses where phi<1
index=index[which(index>-1)]
if (length(index)>1) {
  if((index[length(index)]-index[length(index)-1])>5) {index=index[-(length(index))]}}    #removes last index if it is more than 5 mismatched loci away from next to last locus
phi=phi[index]
index=index-1
#Create plot ===================================================================#
pdf(file="Output_Dataset_Power.pdf")
x=0:(length(info[,1])-1)
y1=info[,3]
y2=info[,2]
p1=which(y2==0)
y2=replace(y2,p1,.000000001)
p1=which(y1<y2)
y3=replace(y1,p1,y2[p1])
y2=y2+1
y3=y3+1
par(mar=c(2,4,1,4)+.1,mfrow=c(2,1),mai=c(0.4,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,log10(y3),xlab="",ylab="Number of Pairs",cex=0.000000000000000000000000001,yaxt="n",ylim=c(min(c(log10(y3),log10(y2))),max(c(log10(y3),log10(y2)))))
ats=c(0,1,2,3,4,5,6,7,8,9)
labs=c(1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000)
axis(side=2,ats,labs)
lines(x,log10(y3),lwd=2)
lines(x,log10(y2),lwd=2)
points(x,log10(y3),pch=21,bg="green",cex=2)
points(x,log10(y2),pch=21,bg="blue",cex=2)
legend("bottomright",c("Observed pairs", "Expected false pairs"), pch = c(21,21),pt.bg=c("green","blue"))
yphi=y2/y3
par(mar=c(2,4,1,4)+.1,new=FALSE,mai=c(.8,1,0.2,1),cex.lab=.99,cex=1.05,lwd=2)
plot(x,yphi,xlab="",ylab=expression(Pr(phi)),cex=2,ylim=c(0,1),pch=21,bg="gray")
lines(x,yphi,pch=21,bg="blue",lty=2,lwd=2,,col="darkgray")
points(x,yphi,cex=2,pch=21,bg="gray")
mtext("Number of Mismatching Loci",side=1,line=1.94)
dev.off()
info2=cbind(x,info[,4])
colnames(info2)<-c("Number of Mismatching Loci", "Pr(Phi)")
write.table(info2, file="Output_Pr(Phi)_Bayesian Prior.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
#begin caluclation of distribution for lamda|phi ===============================#
ngtypes=20000000                                                                #20 million seems like plenty for all datasets (but may need to adjust at some point)
inreps=10000
repnumber=round(100000/(inreps))                                                #this is the rep number to get 100,000 values.  can adjust accordingly
for (n in 1:repnumber) {
    OUT=NULL
    for (r in (sort(rep(unique(OUTALL[,1]),2)))) {
        Loco=OUTALL[which(OUTALL[,1]==r),]
        alleles3=cbind(Loco,ngtypes*as.numeric(as.character(Loco[,4])))
        findo=which(alleles3[,2]==0)                                            #replace 0 with 1 (obsolete if removing 0 works, 2 lines down)
        findo2=replace(alleles3[,4],findo,1)
        alleles3=cbind(alleles3,findo2)
        alleles3=alleles3[-which(alleles3[,2]==0),]
        gtrue=sample(alleles3[,6],inreps,prob=alleles3[,4],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        }
    for (i in index) {                                                              #loop over numbers of mismatched loci
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                    #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distp=cbind(i,DIST)
    write.table(distp,file="False_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
    }
}
#Begin calculation of observed shared freqs  (empirical obs used in lamda|phi) and actual alleles==#
Adults1<- MOMS
a=ncol(Adults1)
Adults=Adults1[,c(2:a)]
Numberloci=ncol(Adults)
Offspring1<- DADS
Offspring=Offspring1[,c(2:a)]
Juvs<- OFFSPRING
Juv=Juvs[,c(2:a)]
ncolP=ncol(Putative)
Putative=cbind(Putative,trunc(NL2-(Putative[,ncolP]*2)))
OUT9=NULL
for (n in index){
  Putative2=Putative[which(Putative[,ncolP+1]==n),]
  m1=match(Putative2[,3],Juvs[,1])
  Putative3=Juvs[m1,]
  write.table(cbind(Putative2[,1:3],Putative3), file="Putative.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
      #if (length(distp[,1])>1000) {                                             #need at least 1000 values , else assigned 0 (may be redundant now)
      OBS=NULL
      afreqs=function(afreqs) {
        PutL=cbind(Putative3[,L+1],Putative3[,L+2])
        alleles2=OUTALL[which(OUTALL[,1]==L),]
        alfreq1=alleles2[match(PutL[,1],alleles2[,2]),4]
        alfreq2=alleles2[match(PutL[,2],alleles2[,2]),4]
        Puts5=cbind(alfreq1,alfreq2)                                            #compare head(cbind(Puts3,Puts4,Puts5)) to alleles 2 as a check on the above
        R1=replace(Puts5[,1],which(is.na(Puts5[,1])==TRUE),1)                   #if a mismatch, every column should be a "1"  (thus probability unaffected)
        R2=replace(Puts5[,2],which(is.na(Puts5[,2])==TRUE),1)
        Puts6=cbind(R1,R2)
        write.table(cbind(L,PutL),file="Putpair.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
        Put_share3=Puts6
        OBS <<- cbind(OBS,Put_share3)
    }
    L=NL2
    C1=for(L in (2*(unique(round((1:(L-2))/2)))+1)) lapply(L,afreqs)
    obsp=apply(OBS, 1, prod)
  #}  else obsp=rep(0,(length(Putative2[,1]))/2)
OUT9 <- rbind(OUT9,cbind(n,obsp))                                               #shared alleles (by freq of chance of sharing an allele).  empirical obs used in lamda|phi
}
#calculate actual shared alleles (empirical) and straight-up allele freqs (used in lamda|phic) ==#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 4, title=paste("Calculating Posterior Component 1"))}
for (n in index) {
Putative2=Putative[which(Putative[,ncol(Putative)]==n),]
m1=match(Putative2[,3],Juvs[,1])
OBA3=Juvs[m1,-1]
for (yyy in 1:10) {                                                               #now set at 100,000 [same as for false pairs)
      OUT=NULL
      for (r in (sort(rep(unique(OUTALL[,1]),2)))) {                            #This first section calculates products of parental alleles (True distribution)
        vect=c(Adults[,r],Adults[,r+1])                                         #currently is only calculating allele frequencies from the adults (could be good if unequal reproductive success)
        alleles=data.frame(table(vect))
        alleles=alleles[order(alleles[,1]),]
        if (as.numeric(as.character(alleles[1,1]))<=0) {alleles=alleles[-1,]}
        alleles2=cbind(alleles,alleles[,2]/sum(alleles[,2]))
        gtrue=sample(alleles2[,3],10000,prob=alleles2[,3],replace=TRUE)
        OUT <- cbind(OUT,gtrue)
        if(yyy==1) {    for (i in 1:length(OBA3[,1])) {                           #this inset finds the frequency of the shared allele    (actual freq as aoppossed to shared freq)
                            mm=alleles2[match(OBA3[i,r],alleles2[,1]),3]
                            mn=alleles2[match(OBA3[i,r+1],alleles2[,1]),3]
                            mm=cbind(r,mm)
                            mn=cbind(r+1,mn)
                            outs=rbind(mm,mn)
                            write.table(outs,file="Shared_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)   #NA situation will make this slightly more conservative
                            }
                }
        }
    for (i in index) {
    if (i==0) {DIST=as.data.frame(apply(OUT, 1, prod))} else {
    DIST=NULL                                                                   #sample distribution by appropriate number of loci
    distp2=as.matrix(OUT)
    a1=NULL
    a2=NULL
    for (z in 1:length(distp2[,1])) {a1=rbind(a1,sample(1:NL,i,replace=F))}     #prevents same locus being sampled twice   (ramdom sampling assumes equal prob of errors)
    for (p in 1:i) {a2=rbind(a2,cbind(1:length(distp2[,1]),a1[,p]))}            #deals with formatting
    distp2[a2]<-1
    a3=apply(distp2, 1, prod)
    DIST<-as.data.frame(a3)
    }
    distt<-cbind(i,DIST)
    write.table(distt,file="True_allele_freqs.txt",row.names=FALSE,col.names=FALSE,sep="\t",append=TRUE)
  }
}
}
#Calculate lamdaphi ============================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 5, title=paste("Calculating Posterior Component 2"))}
Putative3<- read.table("Putative.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Putadults=Putative3[,-(1:4)]                                                    #Combine putative pairs alleles into a single row
Names=Putative3[,1:3]
empirical=cbind(Names,OUT9)                                                     #where OUT9 equals observed freqs (really shared freqs)
distp <- read.table("False_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
P2=NULL
for (i in index) {                                                              #loop over numbers of mismatched loci
  empirical2=empirical[which(as.numeric(as.character(empirical[,4]))==i),]
  if(length(empirical2)==4) {empirical2=t(empirical2)}                          #deals with one indiviudal formatting
  if (empirical2[1,4]==0) {P=empirical2}   else{                                #deals with not enough reps
     a3=distp[which(distp[,1]==i),2]
     DIST<-as.data.frame(a3)
    P=NULL
    for (b in 1:length(empirical2[,1])) {
      p1=length(which(DIST[,1] <= as.numeric(empirical2[b,5]) ))
      if (p1==0) {p1=0.00001}
      p2=cbind(as.character(empirical2[b,1]),as.character(empirical2[b,2]),as.character(empirical2[b,3]),p1)
      p3=cbind(p2,as.numeric(p2[,4])/length(DIST[,1]))
      P <- rbind(P,p3)
      }
    }
  P=cbind(i,P)
  if ((length(P2))>0) {colnames(P2) <- c("1","2","3","4","5","6")}
  colnames(P) <- c("1","2","3","4","5","6")
  P2<-rbind(P2,P)
}
lamdaphi=as.numeric(as.character(P2[,6]))
#Calculate lamda|phic ==========================================================#
if(.Platform$OS.type=="windows"){setWinProgressBar(pb, 6, title=paste("Calculating Posterior Component 3 and Finishing Calculations"))}
lamdaphic_dist<- read.table("True_allele_freqs.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
Observed<- read.table("Shared_allele_freqs.txt", header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
mm1=which(is.na(Observed[,2])==TRUE)                                            #replace NAs with 1
Observed[mm1,2] <- 1                                                            #replace NAs with 1
a=unique(Observed[,1])
U=NULL
for (i in a) {
u=Observed[Observed[,1]==i,2]
U=cbind(U,u)
}
U=U[1:(length(U[,1])/2),]
lamdaphic=apply(U, 1, prod)
l1=length(which(OUT9[,2]==0))
if (l1>0) lamdaphic=c(rep(0,l1),lamdaphic)                                      #match up p-values (not the best way, could get messy with 0'ss)
P3=cbind(P2,lamdaphic)
for (i in index) {                                                              #loop over numbers of mismatched loci
   e2=P3[which(as.numeric(as.character(P3[,1]))==i),]
   if(length(e2)==6) {e2=t(e2)}                                                 #deals with one indiviudal formatting
   if (e2[1,6]==0) {write.table(e2[,6], file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE) }   else{   #deals with not enough reps
    a3=lamdaphic_dist[which(lamdaphic_dist[,1]==i),2]
    DIST<-as.data.frame(a3)
    }
  for (b in 1:length(e2[,1])) {                                                 #calculate p values
    p1=length(which(DIST[,1] <= e2[b,7]))
    p2=p1/length(DIST[,1])
    write.table(p2, file="lamdaphic.txt",sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    }
  }
#}
lamdaphic<- read.table("lamdaphic.txt", header=FALSE, sep="\t", na.strings="390.5", dec=".", strip.white=TRUE)
#Put it all together with Bayes theorem! =======================================#
vals=cbind(P2,lamdaphic[,1])
philength=cbind(0:(length(phi)-1),phi,table(vals[,1]))                          #add phi values to vals
phis=rep(philength[,2],philength[,3])
vals=cbind(vals,phis)
colnames(vals)<-c("Nmismatch","Parent1","Parent2","Offspring","ignore","lamdaphi","lamdaphic","phi")
phi=as.numeric(as.character(vals[,8]))
lamdaphi=as.numeric(as.character(vals[,6]))
lamdaphic=as.numeric(as.character(vals[,7]))
lamdaphi=replace(lamdaphi,which(lamdaphi==0),1)
lamdaphic=replace(lamdaphic,which(lamdaphic==0),1)
pval=(lamdaphi*phi)/((lamdaphi*phi)+(lamdaphic*(1-phi)))                        #pval=replace(pval,which(pval=="NaN"),"< 0.001")
pval=cbind(as.character(vals[,2]),as.character(vals[,3]),as.character(vals[,4]),vals[,1],pval)
pval=pval[order(as.numeric(pval[,5])),]
colnames(pval) <- c("Adult1","Adult2","Offpsring","NL_mismatch","Pvalue")
write.table(pval, file="Output_SOLOMON_Posterior_Probabilities.txt",sep="\t",append=TRUE,col.names=TRUE,row.names=FALSE)
if(.Platform$OS.type=="windows"){close(pb)}
unlink("False_allele_freqs.txt")                                                #clear all sims files
unlink("True_allele_freqs.txt")
unlink("Shared_allele_freqs.txt")
unlink("lamdaphic.txt")
unlink("Putative.txt")
unlink("True_shared_freqs.txt")
unlink("Output_genotypes.txt")
unlink("*.sims")
unlink("pooled.txt")
unlink("Alldata.txt")
unlink("Putpair.txt")
rm(list=ls())
gc()
}
label.bayes <- tklabel(tt, text="Press 'Run' to perform Bayesian parentage analysis:",font=fontTextLabel)
run.button <- tkbutton(tt, text = "Run", command = PressedOK)
tkgrid(tklabel(tt,text="     "),label.bayes,tklabel(tt,text="     "),run.button,tklabel(tt,text="     "))		# Place the button on the window
tkfocus(tt)
tkgrid.configure(label.bayes, sticky="w")
tkgrid.configure(run.button, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
}
label.BAYES2 <- tklabel(tt, text="Perform Bayesian parentage analysis with known parent-pairs:",font=fontTextLabel)
OK.but.BAYES2 <- tkbutton(tt, text = "BAYES", command = PressedBAYES2)
tkgrid(tklabel(tt,text="     "),label.BAYES2,tklabel(tt,text="     "),OK.but.BAYES2,tklabel(tt,text="     "))		# Place the button on the window
tkgrid.configure(label.BAYES2, sticky="w")
tkgrid.configure(OK.but.BAYES2, sticky="w")
tkgrid(tklabel(tt,text="     "))
tkgrid(tklabel(tt,text="     "))
tkfocus(tt)
}


environment(solomon) <- solomon.env


#end#################################################################################################################################################################################################################


