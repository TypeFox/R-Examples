# tables.R 
library(SEERaBomb)
library(reshape2)
library(XLConnect) 
# rnD=function(d,n=1) { 
#   d[sapply(d,is.numeric)]=round(d[sapply(d,is.numeric)],n)
#   d
# }

unlink(f<-"~/Results/amlMDS/Tables.xlsx") 
wb <- loadWorkbook(f,create=T) 
tab1="table1_validation"
createSheet(wb, tab1)
# setSheetColor(wb, tab1, XLC$COLOR.GREEN)
# As code validation, we will reproduce results in Lancet Oncol 2011; 12: 353-60
# As a representative example, for Table 1, we will use prostate cancer treatment induced
# bladder cancer risk and breast cancer treatment induced lung cancer risk
# take canc and popsa from here to avoid redefining trt and AML (to include APL)
system.time(load("~/Results/amlMDS/pm20.RData")) # 2 secs to load. 
system.time(load("~/Results/amlMDS/pf20.RData")) # 3 secs to load. 
mkExcel(pm20,"b5",outDir="~/Results/amlMDS") # makes male excel file for all 1sts and 2nds.
mkExcel(pf20,"b5",outDir="~/Results/amlMDS") # makes female file (sheets are 1sts, rows 2nds). 
# read in the two excel files just made
wbM <- loadWorkbook("~/Results/amlMDS/pMs20e85b5.xlsx") # aka validationMales.xlsx
wbF <- loadWorkbook("~/Results/amlMDS/pFs20e85b5.xlsx") # aka validationFemales.xlsx
(dm=readWorksheet(wbM,"prostate")) # and pull off sheets of interest
(df=readWorksheet(wbF,"breast"))

mkrow=function (df) {cbind(
  df%>%filter(X2nd.cancer=="bladder")%>%select(noRad2bladder=X.5.100..after.noRad),
  df%>%filter(X2nd.cancer=="lung")%>%select(noRad2lung=X.5.100..after.noRad),
  df%>%filter(X2nd.cancer=="bladder")%>%select(rad2bladder=X.5.100..after.rad),
  df%>%filter(X2nd.cancer=="lung")%>%select(rad2lung=X.5.100..after.rad))}
D=rbind(mkrow(df),mkrow(dm))
writeWorksheet(wb,"Table 1. Software validation based on RR* of second bladder- and lung cancers in SEER-9",
               sheet = tab1,header=FALSE)
# writeWorksheet(wb,cbind(first=c("breast","prostate"),D), sheet = tab1,rownames=1, startRow=2)
writeWorksheet(wb,cbind(First.Cancer=c("breast","prostate"),D), sheet = tab1,startRow=2)
for (j in 1:(dim(D)[2]+1)) setColumnWidth(wb,sheet = tab1, column = j, width = 5300)
# save workbook of Tables, including Table 1, only after making all the sheets. 

##### setup for other Tables
system.time(load("~/Results/amlMDS/pm.RData")) # 4 secs to load. 
system.time(load("~/Results/amlMDS/pf.RData")) # 4 secs to load. 
names(pm)
(pm$secondS)
names(pm$L)
#the following heme malignancies (HM) first cancer are excluded from Tables 2 and 3. 
HM=c("AML","MDS","CMML","CML","MPN","ALL","CLL","SLL","HCL","OL","NHL","MM","hodgkin") 

# the following myeloid first cancers are excluded from Tables 4-8
ML=c("AML","MDS","CMML","CML","MPN") # remove MyeLoid firsts. Interest in age correlations with de novo AML

dm=mkDF(pm,"b0_0.25_1_12")
df=mkDF(pf,"b0_0.25_1_12")
d=rbind(cbind(df,Sex="Female"),cbind(dm,Sex="Male"))
d=d%>%filter(!cancer1%in%HM)%>%group_by(trt,cancer2,Sex,int)%>%summarize(O=sum(O),E=sum(E),t=weighted.mean(t,py,na.rm=T))
d=d%>%mutate(RR=O/E, L=qchisq(.025,2*O)/(2*E),U=qchisq(.975,2*O+2)/(2*E))
d=d%>%mutate(RR=paste0(sprintf("%.2f",RR)," (",sprintf("%.2f",L),", ",sprintf("%.2f",U),")"))%>%select(-L,-U,-t)
D1=d%>%filter(int=="(1,12]")%>%select(-int)
D2=d%>%filter(int=="(12,100]")%>%select(-int)
D1$interval="1-12"
D2$interval=">12"
D=rbind(D1,D2)
D=as.data.frame(D)

tab="table2"
createSheet(wb, name = tab)
writeWorksheet(wb, "Table 2. AML/MDS RR after any non-hematological 1st cancer", sheet = tab,header=FALSE)
writeWorksheet(wb, D, sheet = tab,startRow=2)
for (i in 1:5) setColumnWidth(wb,sheet = tab, column = i, width = 2500)
setColumnWidth(wb,sheet = tab, column = 6, width = 4000)
setColumnWidth(wb,sheet = tab, column = 7, width = 3000)

####################### TABLE 3 ##############
mkDCSA=function(df,dm,N=30,intv,trtv) { # cancer specific with ages  # intv= "(1,12]";N=30;trtv="rad
  d=rbind(cbind(df,Sex="Female"),cbind(dm,Sex="Male"))
  d=d%>%filter(!cancer1%in%ML,trt==trtv,O>=N,int==intv,(RR>1&rrL>1)|(RR<1&rrU<1))%>%arrange(cancer2,desc(RR))%>%
    select(cancer1,cancer2,Sex,Age=ageO,O,E,RR)
  d=d%>%mutate(RR=O/E, L=qchisq(.025,2*O)/(2*E),U=qchisq(.975,2*O+2)/(2*E))
  d%>%mutate(RR=paste0(sprintf("%.2f",RR)," (",sprintf("%.2f",L),", ",sprintf("%.2f",U),")"))%>%select(-L,-U)
}

tab="table3"
createSheet(wb, name = tab) 
(D1=mkDCSA(df,dm,30,"(1,12]","rad"))
(D2=mkDCSA(df,dm,50,"(1,12]","noRad"))
D1$radiation="Yes"
D2$radiation="No"
D=rbind(D1,D2)
ct=cor.test(D$Age,D$O/D$E,method="spearman")
(rho=paste0("rho=",round(ct$estimate,3),"; P = ",format.pval(ct$p.value,digits=2)))
writeWorksheet(wb,"Table 3. AML/MDS RR 1-12 years after non-myeloid first cancers",
               sheet = tab,header=FALSE)
writeWorksheet(wb, D, startRow=2,sheet = tab)
writeWorksheet(wb, rho, startCol=9,startRow=3,sheet = tab, header=FALSE)
setcols=function(wb,tab) {
  for (j in 1:6) setColumnWidth(wb,sheet = tab, column = j, width = 2500)
  setColumnWidth(wb,sheet = tab, column = 7, width = 4500)
  setColumnWidth(wb,sheet = tab, column = 8, width = 2000)
#   setColumnWidth(wb,sheet = tab, column = 8, width = 4500) # for rhos
}
setcols(wb,tab)

tab="table4"  # first and last time points
createSheet(wb, name = tab) 
(D2r=mkDCSA(df,dm,5,"(12,100]","rad"))
(D2nr=mkDCSA(df,dm,10,"(12,100]","noRad"))
(D1r=mkDCSA(df,dm,20*1.9/5,"(0,0.25]","rad"))
(D1nr=mkDCSA(df,dm,N=20,"(0,0.25]","noRad"))
D1r$radiation="Yes"
D1nr$radiation="No"
D2r$radiation="Yes"
D2nr$radiation="No"
D1=rbind(D1nr,D1r)
D2=rbind(D2nr,D2r)
D1$interval="<0.25"
D2$interval=">12"
D=rbind(D1,D2)%>%select(-Age)
writeWorksheet(wb,"Table 4. AML/MDS initial and Steady State RR after non-myeloid first cancers",
               sheet = tab,header=FALSE)
writeWorksheet(wb, D, sheet = tab,startRow=2)
setcols(wb,tab)
for (j in 1:5) setColumnWidth(wb,sheet = tab, column = j, width = 2500)
setColumnWidth(wb,sheet = tab, column = 6, width = 4500)
setColumnWidth(wb,sheet = tab, column = 7, width = 3000)
setColumnWidth(wb,sheet = tab, column = 8, width = 3000)
saveWorkbook(wb)


####### this makes an excel file with summary sheets/Tables of SEER cases by year and ICD-O3
if (1) {
  library(XLConnect) 
  library(dplyr) 
  library(scales)    
  cf=function (x) comma_format(width=17,justify="right")(x)
  
  load("~/data/SEER/mrgd/cancDef.RData")
  load("~/data/SEER/mrgd/popsae.RData")
  head(popsae)
  PY=(popsae%>%group_by(year)%>%summarize(PY=sum(py)))[[2]]
  unlink(f<-"~/Results/amlMDS/SEERdataOverView.xlsx") 
  wb <- loadWorkbook(f,create=T) 
  
  createSheet(wb, name = "yearXcancer") 
  yrca=table(canc$yrdx,canc$cancer)
  yrTot=apply(yrca,1,sum)
  Total=apply(yrca,2,sum)
  gtot=sum(Total)
  D=data.frame(year=c(1973:2013,NA),as.data.frame(rbind(yrca,Total)),total=c(yrTot,gtot),PY=cf(c(PY,sum(PY))))
  writeWorksheet(wb,D, sheet = "yearXcancer")
  setColumnWidth(wb,sheet = "yearXcancer", column = dim(D)[2], width = 3700)
  createFreezePane(wb,"yearXcancer",2,2) #sheet colSplit rowSplit
  
  createSheet(wb, name = "histo3Xcancer") 
  total=table(canc$cancer)
  h3ca=table(canc$histo3,canc$cancer)
  D=as.data.frame(cbind(hist03=c(as.numeric(row.names(h3ca)),NA),rbind(h3ca,total)))
  writeWorksheet(wb,D, sheet = "histo3Xcancer")
  createFreezePane(wb, "histo3Xcancer",2,2) 
  
  createSheet(wb, name = "yearXhisto3") 
  total=table(canc$histo3)
  yrh3=table(canc$yrdx,canc$histo3)
  D=as.data.frame(cbind(year=c(1973:2013,NA),rbind(yrh3,total)))
  writeWorksheet(wb,D, sheet = "yearXhisto3")
  createFreezePane(wb, "yearXhisto3",2,2) 
  
  saveWorkbook(wb)
}

