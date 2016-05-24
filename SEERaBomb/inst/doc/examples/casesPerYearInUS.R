# the point here is to map incidences into estimates of how many americans get the cancer each year.
.seerHome="~/data/SEER" 
rm(list=ls()) # note: dot variables defined above persist through such cleanings
getUS=function(d,pop) {
  (pyM=pop$py[1])
  (pyF=pop$py[2])
  d$agerec=as.numeric(cut(d$agedx,breaks=c(0,1,seq(5,85,5),130)))
  (d<-ddply(d, .(sex), summarise,cases=length(agerec))) 
  (M=d$cases[1])
  (F=d$cases[2])
  (crudeIm=M/pyM)
  (crudeIf=F/pyF)
  (malesPerYear=151*crudeIm)
  (femalesPerYear=178*crudeIm)
  (malesPerYear+femalesPerYear)*1e6
} 

# NOTICE: you must
if (0) {  # switch this to 1 (i.e. run this chunk) if you never ran it before
  library(SEERaBomb)
  (df=getFields())
  (df=pickFields(df))
  mkSEERold(df,dataset="00")  
}

library(RSQLite)
m=dbDriver("SQLite")
con=dbConnect(m,dbname=file.path(.seerHome,"00/all.db"))
dbListTables(con)
pops=dbGetQuery(con, "SELECT * from pops")
library(plyr)
(pop<-ddply(pops, .(popsex), summarise,py=sum(population)))

dbListFields(con,"lymyleuk")

d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='CML'") 
getUS(d,pop) # 4600 CML/year
d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='CML' and COD=35022") 
getUS(d,pop) # 600/year=CML deaths

d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='CMML'") 
getUS(d,pop) # 1400 CMML/year
d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='CMML' and COD=35022") 
getUS(d,pop) # 116/year=CMML deaths.
d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='CMML' and COD=35021") 
getUS(d,pop) # 146/year go down as AML
d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='CMML' and COD=35023") 
getUS(d,pop) # 253/year go down as other myeloid cancers

d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='AML'") 
getUS(d,pop) # 12.5k AML/year
d=dbGetQuery(con, "SELECT * from lymyleuk where cancer='AML' and COD=35021") 
getUS(d,pop) # 6k/year=AML deaths

d=dbGetQuery(con, "SELECT * from other where cancer='MDS'") 
getUS(d,pop) # 15k MDS/year


# d=dbGetQuery(con, "SELECT * from lymyleuk where histo3=9863") # 3634/year CML not coded as 9875 or 9876  
# d=dbGetQuery(con, "SELECT * from lymyleuk where histo3=9875") # 834/year BCR-ABL+ CML 
# d=dbGetQuery(con, "SELECT * from lymyleuk where histo3=9876") # 49/year BCR-ABL neg CML 
# d=dbGetQuery(con, "SELECT * from lymyleuk where ICD9=2051") # 5924/year =CML+CMML
# d=dbGetQuery(con, "SELECT * from lymyleuk where histo3==9945") # 1371/year =CMML
d=dbGetQuery(con, "SELECT * from lymyleuk where ICD9=2050") # 12437/year = AML
# d=dbGetQuery(con, "SELECT * from lymyleuk where ICD9=2070") # 260/year = erythroleukemia
# d=dbGetQuery(con, "SELECT * from lymyleuk where ICD9=2072") # 131/year = Megakaryocytic leukemia
d=dbGetQuery(con, "SELECT * from lymyleuk where ICD9=2050 and COD>=20010 and COD<=37000")#8510/year 
d=dbGetQuery(con, "SELECT * from lymyleuk where ICD9=2050 and COD=35021") # 5929/year = AML deaths
d=dbGetQuery(con, "SELECT * from lymyleuk where COD=35021") # 7465/year = AML deaths 
length(unique(d$casenum)) # 9388/9541 = 98.3% unique
d=dbGetQuery(con, "SELECT * from lymyleuk where COD=35023") # 607/year = other myeloid deaths 
d=dbGetQuery(con, "SELECT * from other where COD=35021") # 1379/year AML death  other e.g. MDS incidence 
d=dbGetQuery(con, "SELECT * from other where COD=35022") # 109/year CML death other incidence 
# I'm not seeing a COD code for MDS in http://seer.cancer.gov/codrecode/1969+_d09172004/index.html
getUS(d,pop)



d=dbGetQuery(con, "SELECT * from other where cancer='MDS' and COD=35022") 
getUS(d,pop) # 44/year deaths

