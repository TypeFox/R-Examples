# mkExcelDemo.R
rm(list=ls()) 
library(SEERaBomb) 
if (1) {
  load("~/data/SEER/mrgd/cancDef.RData") 
  load("~/data/SEER/mrgd/popsae.RData") 
  canc=canc%>%select(-reg,-COD,-radiatn,-histo3,-ICD9)
  canc=canc%>%filter(cancer!="benign")
  popsa=popsae%>%group_by(db,race,sex,age,year)%>%summarize(py=sum(py)) # sum on regs
  m=seerSet(canc,popsa,Sex="male",ageStart=0,ageEnd=100) 
  f=seerSet(canc,popsa,Sex="female",ageStart=0,ageEnd=100) 
  m=mk2D(m) 
  f=mk2D(f) 
  brks=c(0,0.5,1,2,3,10)
  m=tsd(m,brks=brks,trts=c("rad","noRad")) 
  f=tsd(f,brks=brks,trts=c("rad","noRad"))
  system.time(save(m,f,file="~/Results/amlMDS/mfExcel.RData")) #~10 seconds 
} else {
  load("~/Results/amlMDS/mfExcel.RData") 
}

mkExcel(m,"b0_0.5_1_2_3_10",outDir="~/Results/amlMDS",outName="males") #out filenames are otherwise coded. 
mkExcel(f,"b0_0.5_1_2_3_10",outDir="~/Results/amlMDS",outName="females")
mkExcel(m,"b0_0.5_1_2_3_10",outDir="~/Results/amlMDS",outName="males",flip=T)
mkExcel(f,"b0_0.5_1_2_3_10",outDir="~/Results/amlMDS",outName="females",flip=T)
