#library(gplots)
#Dir=getwd()
#controlfile="Control_Probe_Profile_12_Samples_Repeat_011111.txt"
#sampfile="Samples_Table_12_Samples_Repeat_011111.txt"

###############################################
###############   Import data   ###############
cat("Please define the following three string variables:\n1. directory path which include the control profile file and sample table file;\n2. controlprofile file name;\n3. sampletable file name:\n")
 controlName=c("Bisulfite_Control","Extension_Control","Hybridization_Control","Target_Removel_Background","Negative_Control","Non-Polymorphic_Control","Stain_Control","Specificity_Control")
 colorList=list()
 colorList[[controlName[1]]]= colors()[c(139,28,34,84)]
 colorList[[controlName[2]]]= colors()[c(34,84,28,139)]
 colorList[[controlName[3]]]= colors()[c(153,28,139)]
 colorList[[controlName[4]]]= colors()[139]
 colorList[[controlName[5]]]= colors()[c(34,76,146,142,85,47,43,29,30,62,99,452,392,411,61,163)]
 colorList[[controlName[6]]]= colors()[c(34,84,28,139)]
 colorList[[controlName[7]]]= colors()[c(34,84,139,28)]
 colorList[[controlName[8]]]= colors()[c(34,84,139,28)]

#Function ImportData(Dir,controlfile,sampfile)
ImportData <- function(Dir,controlfile,sampfile) {
rHomedir =R.home(component="home");
BeadName =paste(rHomedir,"\\library\\Meth27QC\\datafile\\BeadtypeIDs.txt",sep=""); BeadName

setwd(Dir)
dir.create("QCresults")
WorkFiles <- list.files(path = ".",pattern = ".*.txt")


#Ctrl <- WorkFiles[grep("Control", WorkFiles)]
#SamplesName <- WorkFiles[grep("Sample", WorkFiles)]
Ctrl <- paste(Dir,controlfile,sep="/")
SamplesName <- paste(Dir,sampfile,sep="/")

try(Discarder <- WorkFiles[grep("Discard", WorkFiles)], silent=TRUE)

control <- read.table(Ctrl, header = TRUE, sep = "\t",check.names=FALSE)
#AvBeta <- read.table(AverageBeta, header = TRUE, sep = "\t", flush=T)
samps <- read.table(SamplesName, header = TRUE, sep = "\t",check.names=FALSE)
BeadIds <- read.table(BeadName, header = TRUE, sep = " ",check.names=FALSE)
#DiscarderII <- ""
#try(DiscarderII <- read.table(Discarder, header = FALSE, sep = "\t"), silent=TRUE)
#
newdir=paste(Dir,"/QCresults",sep="")
setwd(newdir)
samps$SampleLabel <-  samps[["Sample ID"]]


colnames(samps)[1] <- "Index"
colnames(samps)[2] <- "SampleID"
sampfn= "sample.txt"
write.table(samps,sampfn, sep="\t",row.names=FALSE)
samps_mod <- read.table(sampfn, header = TRUE, sep = "\t")

samps2 <- samps[with(samps, order(SampleID)), ]
nsample <- length(samps$Index)
 fname= "Sample.pdf"
pdf(file=fname,paper="a4r", fonts="Times")
SamplePDF <- data.frame(samps_mod$Index[1:(length(samps_mod$Index)/2)], samps_mod$SampleID[1:(length(samps_mod$SampleID)/2)],rep(" ",(length(samps_mod$SampleID)/2)),
samps_mod$Index[((length(samps_mod$Index)/2)+1):length(samps_mod$Index)], samps_mod$SampleID[((length(samps_mod$SampleID)/2)+1):length(samps_mod$SampleID)])


colnames(SamplePDF)<-c("Index","SampleID","","Index","SampleID")
textplot(SamplePDF, halign="center", valign="center", show.rownames = FALSE)
title("Sample List")
dev.off()

return(list(ctrl=control,beadIDs =BeadIds,samples=samps,Ctrl=Ctrl, samps_mod=samps_mod, nsample=nsample, samps2=samps2))
}


####################################################
###############   Internal Control   ###############
QCRep <- function(Dir,controlfile,sampfile) {
	
dataFiles <- ImportData(Dir,controlfile,sampfile)
control <- dataFiles$ctrl
samps <- dataFiles$samples
BeadIds <- dataFiles$beadIDs
BeadIds =BeadIds[,-grep("Evaluate",names(BeadIds))]
#####1.SampelTable analysis
#---------------------------

sampNames = names(samps);sampNames
namIds=grep("(detect.*gene)|(sig.*Average)",sampNames,ignore.case=TRUE,);namIds
name2Plot=sampNames[namIds];name2Plot
sampIds= samps$SampleID
fName ="Sample Intensity Report.pdf"
pdf(file=fName,paper="a4r",width =15,height=7, fonts="Times")
for(n in 1:length(namIds))  {

plotData=samps[[name2Plot[n]]];
par(xpd=TRUE,mar=c(5,4,4,8))                      #
plot(plotData,ann=FALSE,tcl=-0.3,cex.axis=0.7)
axislen=length(plotData)+10
axis(1,at=0:axislen,lab=0:axislen,cex.axis =0.7,tcl=-0.2)
titlename1=paste("Intensity of ",name2Plot[n],sep="")         #
 # Create a title with a red, bold/italic font
title(main=titlename1, col.main="red", font.main=4,)
title(xlab= name2Plot[n], col.lab=rgb(0,0.5,0),line=3)    #
title(ylab= "Intensity", col.lab=rgb(0,0.5,0),line=3)
#Flag the outliers, which is lower than mean -sd
summa =summary(plotData)
IQr=summa[5]-summa[2]
down =summa[2] -1.5*IQr
up = summa[5] +1.5*IQr
outlierids=which((plotData>up)|(plotData<down))
if (length(outlierids) >0) {
  points(outlierids,plotData[outlierids], pch=20,col="red")
  text(outlierids,plotData[outlierids],sampIds[outlierids], cex=0.5, pos=4, col="blue")
  legend("bottomleft","outliers",fill="red",cex=0.7,bty="n")
 }
boxplot(plotData,add =TRUE,at=length(plotData)+7,axes =FALSE,col ="yellow",varwidth =TRUE, boxwex =3)
}

dev.off()
##2. Compare Detection.Pval and expected value,Bead_Type_ID
##-------------------------------------------------------------
Bead_Type_ID <- control$ProbeID
conPvalue <- control[,grep("Pval",names(control))]
tPvalue = cbind(Bead_Type_ID,conPvalue)
#mergTbl=merge(tPvalue,BeadIds,by="Bead_Type_ID")
mergTbl=merge(tPvalue,BeadIds,by="Bead_Type_ID")
pvalCols=grep("Pval",names(mergTbl))
len =length(pvalCols)
expectedId =grep("Expected",names(mergTbl))

#Define the subset of the not match table(NoMatchtbl)
indexofID =grep("ID",names(mergTbl))
indexofDesp=grep("Description",names(mergTbl))
 listData =list()
 mergTbl[[names(mergTbl)[expectedId]]]=as.character(mergTbl[[names(mergTbl)[expectedId]]])
 mergTbl[[names(mergTbl)[indexofDesp]]]=as.character(mergTbl[[names(mergTbl)[indexofDesp]]])
 mergTbl[[names(mergTbl)[indexofID]]]=as.character(mergTbl[[names(mergTbl)[indexofID]]])

#parse every P_value in the merged table, change the value to High(<0.05), Mediumn, Low, and Backgraound according the P_value
for (n in 1:len)
   {   i=pvalCols[n]
        #print(i)
        noMism=TRUE
        proID =NULL;
        exactpValue =NULL
        descp=NULL
        expValue =NULL
        sampleName=""
      for (j in 1:nrow(mergTbl))
         {
            #print(j)
             cellv=mergTbl[j,i]    #Pvalue of ith Sample, j ProbeID,parse down through the same sample first
             expv=mergTbl[j,expectedId]
                if (cellv <= 0.05)                 #Pvalue correspond to High,Medium and Low
             {  #compare the expected pvalue of the jth row,if is Background, add to mismatch table
                  if (expv =="Background")
                  {
                           noMism=FALSE
                          namsplit= unlist(strsplit(names(mergTbl)[i],".",fixed=TRUE))
                          sampleName=namsplit[1]
                              pid= mergTbl[j,indexofID];#print(pid)
                              
                         desid=mergTbl[j,indexofDesp];#print(desid)
                         proID =c(proID,pid)
                         descp=c(descp,desid)
                         exactpValue=c(exactpValue,cellv)
                         expValue =c(expValue,expv)
                                       }
             }
            else                             #Pvalue correspond to Background
            {
                if  (expv !="Background")    #if not background, add to mismatch table
                {
                          noMism=FALSE
                          namsplit= unlist(strsplit(names(mergTbl)[i],".",fixed=TRUE))
                          sampleName=namsplit[1]
                        # print("Sample ")
                          #print(sampleName)

                      #print("Not match at")
                         pid= mergTbl[j,indexofID];#print(pid);
                         desid=mergTbl[j,indexofDesp];#print(desid)
                         proID =c(proID,pid)
                         descp=c(descp,desid)
                         exactpValue=c(exactpValue,cellv)
                          expValue =c(expValue,expv)

                }
            }
          #print("This cell is finished")
         }
   if (noMism==FALSE){listData[[sampleName]]=data.frame(sampleID=sampleName,probeID =proID,Description =descp, Detected_Pvalue =exactpValue,Expected_Intensity=expValue)  }

   }

#Build the mismatch table(NoMatchtbl)
outTbl=data.frame(sampleID="",probeID ="",Description ="", Detected_Pvalue =NA,Expected_Intensity="")
tblName="pvalue_misMatchtable.txt"
for (k in 1:length(listData)) {
        dataf=as.data.frame(listData[[k]])
        outTbl=rbind(outTbl,dataf)
          }
# with(outTbl,order(colnames=()))
 write.table(outTbl,file=tblName,append=FALSE,quote=FALSE,sep="\t\t");
 totalIDlist=samps$SampleID
 misMatchids=match(names(listData),totalIDlist)
 matchSample=as.character(totalIDlist[-misMatchids])
#print(c("P_value mismatch table is saved to file \'pvalue_misMatchtable.txt\' in the subdirectory \'QCresults\'."))
#print("The match sample ID is:" )
#print(matchSample)
#####2.SampelTable analysis
#---------------------------


### 3.QUALITY CHECK
#Quality Check- Bisulfite Conversion
#--------------------------------------------
##Bisulfite_grn
targetNames =c("BC conversion_C1","BC conversion_C2","BC conversion_U1","BC conversion_U2")
Bisulfite_Sig1 <- control[control$ProbeID == "5270706",grep("_Grn",names(control))]
##NonPolymorphic.Ag <- control[control$ProbeID == "1740025",seq(4,length(control),3)]
Bisulfite_Sig1 <- as.data.frame(t(Bisulfite_Sig1))
Bisulfite_Sig2 <- control[control$ProbeID == "4670278",grep("_Grn",names(control))]
Bisulfite_Sig2 <- as.data.frame(t(Bisulfite_Sig2))
Bisulfite_Bckg1 <- control[control$ProbeID == "5290048",grep("_Grn",names(control))]
Bisulfite_Bckg1 <- as.data.frame(t(Bisulfite_Bckg1))
Bisulfite_Bckg2 <- control[control$ProbeID == "4670484",grep("_Grn",names(control))]
Bisulfite_Bckg2 <- as.data.frame(t(Bisulfite_Bckg2))

#Combine columns into a dataframe for green control
cg =cbind(Bisulfite_Sig1,Bisulfite_Sig2,Bisulfite_Bckg1,Bisulfite_Bckg2)

##Bis_Red
Bis_red_Sig1<- control[control$ProbeID == "5270706",grep("_Red",names(control))]
Bis_red_Sig1 <- as.data.frame(t(Bis_red_Sig1))
Bis_red_Sig2 <- control[control$ProbeID == "4670278",grep("_Red",names(control))]
Bis_red_Sig2 <- as.data.frame(t(Bis_red_Sig2))
Bis_red_Bckg1 <- control[control$ProbeID == "5290048",grep("_Red",names(control))]
Bis_red_Bckg1 <- as.data.frame(t(Bis_red_Bckg1))
Bis_red_Bckg2 <- control[control$ProbeID == "4670484",grep("_Red",names(control))]
Bis_red_Bckg2 <- as.data.frame(t(Bis_red_Bckg2))
#Combine columns into a dataframe for red control
cr =cbind(Bis_red_Sig1,Bis_red_Sig2,Bis_red_Bckg1,Bis_red_Bckg2)
names(cg)=targetNames
names(cr)=targetNames
bt=controlName[1]   #"BisulfiteControl Conversion"
#Output plot to pdf file
write2Pdf(cg,cr,bt)



##### QUALITY CHECK - Extension
ExNames =c("Extension(A)","Extension(T)","Extension(G)","Extension(C)")
Extension.Ag <- control[control$ProbeID == "360446",grep("_Grn",names(control))]
Extension.Ar <- control[control$ProbeID == "360446",grep("_Red",names(control))]
Extension.Tg <- control[control$ProbeID == "520537",grep("_Grn",names(control))]
Extension.Tr <- control[control$ProbeID == "520537",grep("_Red",names(control))]
Extension.Gg <- control[control$ProbeID == "1190050",grep("_Grn",names(control))]
Extension.Gr <- control[control$ProbeID == "1190050",grep("_Red",names(control))]
Extension.Cg <- control[control$ProbeID == "2630184",grep("_Grn",names(control))]
Extension.Cr <- control[control$ProbeID == "2630184",grep("_Red",names(control))]
#####Combine the A,T,C,G columns to a dataframe for green and Red respectively
Eg=as.data.frame(cbind(t(Extension.Ag),t(Extension.Tg),t(Extension.Gg),t(Extension.Cg)))
Er=as.data.frame(cbind(t(Extension.Ar),t(Extension.Tr),t(Extension.Gr),t(Extension.Cr)))
names(Eg)=ExNames
names(Er)=ExNames
et=controlName[2]         # "Extension Control"
write2Pdf(Eg,Er,et)

##### QUALITY CHECK - Hybridization
HybName =c("Hyb_(Low)","Hyb_(Medium)","Hyb_(High)")
HybL.g <- t(control[control$ProbeID == "2450040",grep("_Grn",names(control))])
HybL.r <- t(control[control$ProbeID == "2450040",grep("_Red",names(control))])
HybM.g <- t(control[control$ProbeID == "5690110",grep("_Grn",names(control))])
HybM.r <- t(control[control$ProbeID == "5690110",grep("_Red",names(control))])
HybH.g <- t(control[control$ProbeID == "5690072",grep("_Grn",names(control))])
HybH.r <- t(control[control$ProbeID == "5690072",grep("_Red",names(control))])

#plot control intensity to pdf file
Hybg=as.data.frame(cbind(HybL.g,HybM.g,HybH.g))
Hybr=as.data.frame(cbind(HybL.r,HybM.r,HybH.r))
names(Hybg) = HybName
names(Hybr) = HybName
ht=controlName[3]              #"Hybridization Control"
#write to pdf file
write2Pdf(Hybg,Hybr,ht)
#####
#### QUALITY CHECK - Target Removal
TMname="Target Removel"
TM <- as.data.frame(t(control[control$ProbeID == "580035",grep("_Grn",names(control))]))
TMg =TM
TMr =as.data.frame(t(control[control$ProbeID == "580035",grep("_Red",names(control))]))

#plot control intensity to pdf file
names(TMg) = TMname
names(TMr) = TMname
tmt=controlName[4]        #"Target Removel Background"
#write to pdf file
write2Pdf(TMg,TMr,tmt)
#######################
#### QUALITY CHECK - Negative
num =c(1:16)
NegName=paste("Negative",as.character(num),sep="")
probeIDlist=c("50110","360079","430114","460494","540577","610692","610706","670750","1190458","1500059", "1500167","1500398","1660097","1770019", "1940364","1990692")
##Import green control
Ng = as.data.frame(t(control[match(probeIDlist,control$ProbeID),grep("_Grn",names(control))]))

####Import RED control
Nr = as.data.frame(t(control[match(probeIDlist,control$ProbeID),grep("_Red",names(control))]))
#plot control intensity to pdf file
names(Ng) = NegName
names(Nr) = NegName
nt=controlName[5]          #"Negative Control"
#write to pdf file
write2Pdf(Ng,Nr,nt)
#end ##################


##### QUALITY CHECK - Non-Polymorphic (A+T / G+C)
NPname =c("NP(A)","NP(T)","NP(G)","NP(C)")
NpIDlist=c("1740025","2480348","110184","2810035")
NPg = as.data.frame(t(control[match(NpIDlist,control$ProbeID),grep("_Grn",names(control))]))
NPr = as.data.frame(t(control[match(NpIDlist,control$ProbeID),grep("_Red",names(control))]))
names(NPg) = NPname
names(NPr) = NPname
npt=controlName[6]           #"Non-Polymorphic Control"
#write to pdf file
write2Pdf(NPg,NPr,npt)
#end#


##### QUALITY CHECK - Staining DNP
#Plot staining control intensity######
Stainname =c("DNP(med)","DNP(bgnd)","Biotin(med)","Biotin(bgnd)")
StainIDlist=c("4200736","5340168","4570020","5050601")
Staing = as.data.frame(t(control[match(StainIDlist,control$ProbeID),grep("_Grn",names(control))]))
Stainr = as.data.frame(t(control[match(StainIDlist,control$ProbeID),grep("_Red",names(control))]))
#
names(Staing) = Stainname
names(Stainr) = Stainname
st=controlName[7]         #"Stain Control"
#write to pdf file
write2Pdf(Staing,Stainr,st)
#end#

##### QUALITY CHECK - Specificity (mismatch 1)
##Plot staining control intensity######
Smisname =c("GT mismatch1(PM)","GT mismatch1(MM)","GT mismatch2(PM)","GT mismatch2(MM)")
SmisIDlist=c("4610725","4610400","3800154","3800086")
Smisg = as.data.frame(t(control[match(SmisIDlist,control$ProbeID),grep("_Grn",names(control))]))
Smisr = as.data.frame(t(control[match(SmisIDlist,control$ProbeID),grep("_Red",names(control))]))
###
names(Smisg) = Smisname
names(Smisr) =Smisname
smt=controlName[8]      #"Specificity Control"
#write to pdf file
write2Pdf(Smisg,Smisr,smt)
#end#

}

######Main Function###
Meth27QC <- function(Dir,controlfile,sampfile) {
QCRep(Dir,controlfile,sampfile)
print("Analyzing is finished! All results were saved in subdirectory named QCresults")
}



####Custom Function 1 - write2Pdf #########
write2Pdf <- function(dg,dr,tt){

t1=paste(tt,"_Green",sep="")
t2=paste(tt,"_Red",sep="")
#Output plot to pdf file
fileName=paste(tt,".pdf",sep="")
pdf(file=fileName,paper="a4r",width =15,height=7, fonts="Times")

if (length(rownames(dg)) <=12) {par(mfrow=c(1,2),pty="m")}
par(xpd=TRUE,mar=c(5,4,4,8))
max_y=max(max(dg),max(dr))
QCplot(dg,t1,max_y,tt)
QCplot(dr,t2,max_y,tt)

###write plotted data to pdf file
splitandplot(dg,t1)
splitandplot(dr,t2)
dev.off()
}
#######Custom Function  2 - QCplot() ################
QCplot <- function(data4Plot,titlename,max_y,tt){
lenCol=length(names(data4Plot));lenRow=length(rownames(data4Plot))

#par(xpd=T,mar=c(5,4,4,11))

plot_colors =colorList[[tt]]
plot(data4Plot[,1], pch=20,col=plot_colors[1], ylim=c(0,max_y), axes=FALSE, ann=FALSE)
axis(1, at=1:lenRow, lab=1:lenRow,cex.axis =0.7)
# Make y axis with horizontal labels that display ticks at 
# every 4 marks. 4*0:max_y is equivalent to c(0,4,8,12).
axis(2, las=1,cex.axis =0.8)
# Create box around plot
box()
if (lenCol>1){
          for(i in 2:lenCol)
          {
           points(data4Plot[,i], pch=20,col=plot_colors[i])
          }
        }
if (lenRow>=40) legnum=3 else legnum=1;
legend(lenRow+legnum,max_y*2/3.5,names(data4Plot),fill=plot_colors,cex=0.7,bty="n")

# Create a title with a red, bold/italic font
title(main=titlename, col.main="red", font.main=4)

# Label the x and y axes with dark green text
title(xlab= "Samples", col.lab=rgb(0,0.5,0),line=2)
title(ylab= "Intensity", col.lab=rgb(0,0.5,0),line=3)
}
###End of QCplot()############
###Function splitdataf(x)
 splitdataf <- function(x) {

 len=nrow(x);len
 framelist=list()
 nm=30
 subnum=floor(len/nm);
 
if(subnum >=1) {
   if (floor(len/2) ==ceiling(len/2)) subs =subnum else subs=subnum+1;

     for (k in 1:subs){
        from = (k-1)*nm +1
        if (k > subnum) to =len else to= k*nm
        framelist[[k]]= x[(from:to),]
     }

   } else framelist[[1]]=x

 return (framelist)
 }

###Function. splitandplot(x,t)
 splitandplot <- function(x,t1) {

#par(cex=2)
 sampID =rownames(x)
 x=cbind(sampID,x)
 listd =splitdataf(x);listd
 lenlist =length(listd)

###splilt columns if column length >=9

            if (ncol(x) >= 9 ) {
                 sapply(listd,function(x){
                                 unitlist= x
                                 outlist1 = unitlist[1:9]
                                 outlist2 = unitlist[c(1,10:length(unitlist))]
                                 textplot(outlist1, halign="center", valign="center",cex=0.8, show.rownames = FALSE)
                                 title(t1)
                                 textplot(outlist2, halign="center", valign="center",cex=0.8, show.rownames = FALSE)
                                 }
                        )
             } else{
                 for (i in 1:lenlist)
                 {
                  textplot(listd[[i]], halign="center", valign="center",cex=1, show.rownames = FALSE)
                  title(t1)
                 }
              }
 }
 ####End of function splitandplot()