###############################
### Préparation des données ###
###############################

### Ce placer dans le repertoire testsDev
if(getwd()=="C:/Documents and Settings/Christophe/Mes documents"){
    setwd("C:/Documents and Settings/Christophe/Mes documents/Recherche/R et SAS/r2lBiv/trunk/Xolh/testsDev")
}else{}
if(getwd()=="C:/Documents and Settings/Administrator/My Documents"){
    setwd("C:/Documents and Settings/Administrator/My Documents/Recherche/R et SAS/yeap/trunk/Xolh/testsDev")
}else{}

if(!file.exists("tmp")){dir.create("tmp")}
setwd("tmp")

load("../../data/AIDSknowledge.rda")
data <- AIDSknowledge

data$studyLevel <- ordered(data$studyLevel,levels=c("L1","L2","L3","M1","M2","Thesis"))
data$personnalSituation <- ordered(data$personnalSituation,levels=c("Single","Single, occasional partner","Couple, not live together","Couple, live together"))
data$howLongTogether <- ordered(data$howLongTogether,levels=c("Less than a week","Less than a month","Less than a year","More thant a year"))
data$religiousPractice <- ordered(data$religiousPractice,levels=c("None","Major ceremonies","Once or twice a month","Once a week","More than once a week"))
data$nbPartnersWho <- ordered(data$nbPartnersWho,levels=c("Long relationships","Several times","All"))
data$useCondom <- ordered(data$useCondom,levels=c("Never","Very rarely","Rarely","Occasionally","Sometimes","Often","Frequently","Very frequently","Always"))
data$surveyIntimacy <- ordered(data$surveyIntimacy,levels=c("Banal","Slithly intimate","Moderately intimate","Very intimate","Too intimate","Shocking"))

for(i in c(17,18,21,36,38:40)+3){data[,i] <- ordered(data[,i],levels=c("No","Do not know","Yes"))}
for(i in c(19,20,22,23,24,25,26,35,37)+3){data[,i] <- ordered(data[,i],levels=c("Yes","Do not know","No"))}
for(i in c(42:45)+3){data[,i] <- ordered(data[,i],levels=c("Never","Rarely","Sometimes","Frequently","Always"))}



#################################
### Univarié, single variable ###
#################################

if(!file.exists("r2lhUnivSingle")){dir.create("r2lhUnivSingle")}
setwd("r2lhUnivSingle")

### logical
r2lu(data$gender,"univLogical.tex",graphDir="graphUniv1")

### Factor
# Small
r2lu(data$sexWithRiskWhy,"univFactorSmall.tex",graphDir="graphUniv2")
# Big
r2lu(data$nativeTown,"univFactorBig.tex",graphDir="graphUniv3")

### Ordered
# Small
r2lu(data$studyLevel,"univOrderedSmall.tex",graphDir="graphUniv4")
# Big
r2lu(data$useCondom,"univOrderedBig.tex",graphDir="graphUniv5")

### Discrete
# Small
r2lu(data$brothersSisters,"univDiscreteSmall.tex",graphDir="graphUniv6")
# big
r2lu(data$firstRelation,"univDiscreteBig.tex",limDiscrete=50,graphDir="graphUniv7")

### Continuous
# Small
r2lu(data$brothersSisters,"univContinuousSmall.tex",limDiscrete=3,graphDir="graphUniv8")
# big
r2lu(data$age,"univContinuousBig.tex",graphDir="graphUniv9")

### Main
r2lMainFile(text="
\\newpage
\\section{Logical}
\\input{univLogical}

\\newpage
\\section{Factor}
\\subsection{Small}
\\input{univFactorSmall}
\\subsection{Big}
\\input{univFactorBig}

\\newpage
\\section{Ordered}
\\subsection{Small}
\\input{univOrderedSmall}
\\subsection{Big}
\\input{univOrderedBig}

\\newpage
\\section{Discrete}
\\subsection{Small}
\\input{univDiscreteSmall}
\\subsection{Big}
\\input{univDiscreteBig}

\\newpage
\\section{Continuous}
\\subsection{Small}
\\input{univContinuousSmall}
\\subsection{Big}
\\input{univContinuousBig}
")

setwd("..")



############################
### Univarié, data.frame ###
############################

if(!file.exists("r2lhUnivDataFrame")){dir.create("r2lhUnivDataFrame")}
setwd("r2lhUnivDataFrame")

### data.frame
text <- paste("\\newpage\n\\section{",names(data)[-1],"}\n",sep="")
r2lu(data[,-1],"univAnalysis.tex",textBefore=text)

### Main
r2lMainFile()

setwd("..")



################################
### Bivarié, single variable ###
################################

if(!file.exists("r2lhBivSingle")){dir.create("r2lhBivSingle")}
setwd("r2lhBivSingle")



######################
### LOGICAL

### Logical~logical
r2lb(data$firstTimeCondom~data$gender,fileOut="ll.tex",graphDir="ll")

### Logical~Factor
# Small-wide
#r2lb(data$firstTimeCondom~data$area,fileOut="lf1.tex",graphDir="lf1")
# Small-long
r2lb(data$firstTimeCondom~data$area,fileOut="lf2.tex",graphDir="lf2",displayStyle="long")

# Big-wide
#r2lb(data$firstTimeCondom~data$nativeTown,fileOut="lf3.tex",graphDir="lf3",displayStyle="wide")
# Big-long
r2lb(data$firstTimeCondom~data$nativeTown,fileOut="lf4.tex",graphDir="lf4",displayStyle="long")

### Logical~Ordered
# Small-long
r2lb(data$firstTimeCondom~data$transWithCondom,fileOut="lo2.tex",graphDir="lo2",displayStyle="long")
# Big-long
r2lb(data$firstTimeCondom~data$useCondom,fileOut="lo4.tex",graphDir="lo4",displayStyle="long")

### Logical~Discrete
# Small
r2lb(data$firstTimeCondom~data$brothersSisters,fileOut="ld2.tex",graphDir="ld2",displayStyle="long")
# Big
r2lb(data$firstTimeCondom~data$firstRelation,fileOut="ld4.tex",graphDir="ld4",displayStyle="long",limDiscreteX=20)

### Logical~Continuous
# Small
r2lb(data$firstTimeCondom~data$brothersSisters,fileOut="lc2.tex",graphDir="lc2",displayStyle="long",limDiscreteX=3)
# Big
r2lb(data$firstTimeCondom~data$firstRelation,fileOut="lc4.tex",graphDir="lc4",displayStyle="long",limDiscreteX=3)



######################
### FACTOR

### Factor~logical
r2lb(data$area~data$gender,fileOut="fl1.tex",graphDir="fl1",displayStyle="long")
r2lb(data$nativeTown~data$gender,fileOut="fl2.tex",graphDir="fl2",displayStyle="long")

### Factor~Factor
# Big-long
r2lb(data$area~data$nativeTown,fileOut="ff4.tex",graphDir="ff4",displayStyle="long")
# Big-big
r2lb(data$sexWithRiskWhy~data$nativeTown,fileOut="ff5.tex",graphDir="ff5",displayStyle="long")

### Factor~Ordered
# Small-long
r2lb(data$area~data$transWithCondom,fileOut="fo2.tex",graphDir="fo2",displayStyle="long")
# Big-long
r2lb(data$area~data$useCondom,fileOut="fo4.tex",graphDir="fo4",displayStyle="long")

### Factor~Discrete
# Small
r2lb(data$area~data$brothersSisters,fileOut="fd2.tex",graphDir="fd2",displayStyle="long")
# Big
r2lb(data$area~data$firstRelation,fileOut="fd4.tex",graphDir="fd4",displayStyle="long",limDiscreteX=20)

### Factor~Continuous
# Small
r2lb(data$area~data$brothersSisters,fileOut="fc2.tex",graphDir="fc2",displayStyle="long",limDiscreteX=3)
# Big
r2lb(data$area~data$firstRelation,fileOut="fc4.tex",graphDir="fc4",displayStyle="long",limDiscreteX=3)



######################
### ORDERED

### Ordered~logical
r2lb(data$howLongTogether~data$gender,fileOut="ol.tex",graphDir="ol")

### Ordered~Factor
# Small-long
r2lb(data$howLongTogether~data$area,fileOut="of2.tex",graphDir="of2",displayStyle="long")
# Big-long
r2lb(data$howLongTogether~data$nativeTown,fileOut="of4.tex",graphDir="of4",displayStyle="long")

### Ordered~Ordered
# Small-long
r2lb(data$howLongTogether~data$transWithCondom,fileOut="oo2.tex",graphDir="oo2",displayStyle="long")
# Big-long
r2lb(data$howLongTogether~data$useCondom,fileOut="oo4.tex",graphDir="oo4",displayStyle="long")

### Ordered~Discrete
# Small
r2lb(data$howLongTogether~data$brothersSisters,fileOut="od2.tex",graphDir="od2",displayStyle="long")
# Big
r2lb(data$howLongTogether~data$firstRelation,fileOut="od4.tex",graphDir="od4",displayStyle="long",limDiscreteX=20)

### Ordered~Continuous
# Small
r2lb(data$howLongTogether~data$brothersSisters,fileOut="oc2.tex",graphDir="oc2",displayStyle="long",limDiscreteX=3)
# Big
r2lb(data$howLongTogether~data$firstRelation,fileOut="oc4.tex",graphDir="oc4",displayStyle="long",limDiscreteX=3)



######################
### DISCRETE

### Discrete~logical
r2lb(data$brothersSisters~data$gender,fileOut="dl.tex",graphDir="dl")

### Discrete~Factor
# Small-long
r2lb(data$brothersSisters~data$area,fileOut="df2.tex",graphDir="df2",displayStyle="long")
# Big-long
r2lb(data$brothersSisters~data$nativeTown,fileOut="df4.tex",graphDir="df4",displayStyle="long")

### Discrete~Ordered
# Small-long
r2lb(data$brothersSisters~data$transWithCondom,fileOut="do2.tex",graphDir="do2",displayStyle="long")
# Big-long
r2lb(data$brothersSisters~data$useCondom,fileOut="do4.tex",graphDir="do4",displayStyle="long")

### Discrete~Discrete
# Small
r2lb(data$brothersSisters~data$firstRelation,fileOut="dd2.tex",graphDir="dd2",displayStyle="long",limDiscreteX=20)
# Big
r2lb(data$firstRelation~data$brothersSisters,fileOut="dd4.tex",graphDir="dd4",displayStyle="long",limDiscreteY=20)

### Discrete~Continuous
# Small
r2lb(data$brothersSisters~data$age,fileOut="dc2.tex",graphDir="dc2",displayStyle="long",limDiscreteX=3)
# Big
r2lb(data$brothersSisters~data$firstRelation,fileOut="dc4.tex",graphDir="dc4",displayStyle="long",limDiscreteX=3)



######################
### CONTINUOUS

### Continuous~logical
r2lb(data$age~data$gender,fileOut="cl.tex",graphDir="dl")

### Continuous~Factor
# Small-long
r2lb(data$age~data$area,fileOut="cf2.tex",graphDir="cf2",displayStyle="long")
# Big-long
r2lb(data$age~data$nativeTown,fileOut="cf4.tex",graphDir="cf4",displayStyle="long")

### Continuous~Ordered
# Small-long
r2lb(data$age~data$transWithCondom,fileOut="co2.tex",graphDir="co2",displayStyle="long")
# Big-long
r2lb(data$age~data$useCondom,fileOut="co4.tex",graphDir="co4",displayStyle="long")

### Continuous~Discrete
# Small
r2lb(data$age~data$firstRelation,fileOut="cd2.tex",graphDir="cd2",displayStyle="long",limDiscreteX=20)
# Big
r2lb(data$age~data$brothersSisters,fileOut="cd4.tex",graphDir="cd4",displayStyle="long",limDiscreteX=20)

### Continuous~Continuous
# Small
r2lb(data$age~data$brothersSisters,fileOut="cc2.tex",graphDir="cc2",displayStyle="long",limDiscreteX=3)
# Big
r2lb(data$age~data$firstRelation,fileOut="cc4.tex",graphDir="cc4",displayStyle="long",limDiscreteX=3)





### Main
r2lMainFile(text="
\\tableofcontents

\\newpage
\\section{Logical}

\\subsection{Logical$\\sim$Logical}
\\input{ll}

\\newpage
\\subsection{Logical$\\sim$Factor}
\\input{lf4}

\\newpage
\\subsection{Logical$\\sim$Ordered}
\\input{lo4}

\\newpage
\\subsection{Logical$\\sim$Discrete}
\\input{ld4}

\\newpage
\\subsection{Logical$\\sim$Continuous}
\\input{lc4}



\newpage
\\section{Factor (3 et +)}
\\subsection{Factor$\\sim$Logical}
\\input{fl1}

\\newpage
\\subsection{Factor$\\sim$Factor}
\\input{ff4}

\\newpage
\\subsection{Factor$\\sim$Ordered}
\\input{fo4}

\\newpage
\\subsection{Factor$\\sim$Discrete}
\\input{fd4}

\\newpage
\\subsection{Factor$\\sim$Continuous}
\\input{fc4}



\\newpage
\\section{Ordered}

\\subsection{Ordered$\\sim$Logical}
\\input{ol}

\\newpage
\\subsection{Ordered$\\sim$Factor}
\\input{of2}

\\newpage
\\subsection{Ordered$\\sim$Ordered}
\\input{oo2}

\\newpage
\\subsection{Ordered$\\sim$Discrete}
\\input{od2}

\\newpage
\\subsection{Ordered$\\sim$Continuous}
\\input{oc4}



\\newpage
\\section{Discrete}

\\subsection{Discretel$\\sim$Logical}
\\input{dl}

\\newpage
\\subsection{Discrete$\\sim$Factor}
\\input{df2}

\\newpage
\\subsection{Discrete$\\sim$Ordered}
\\input{do2}

\\newpage
\\subsection{Discrete$\\sim$Discrete}
\\input{dd2}

\\newpage
\\subsection{Discrete$\\sim$Continuous}
\\input{dc2}



\\newpage
\\section{Continuous}

\\subsection{Continuousl$\\sim$Logical}
\\input{cl}

\\newpage
\\subsection{Continuous$\\sim$Factor}
\\input{cf2}

\\newpage
\\subsection{Continuous$\\sim$Ordered}
\\input{co2}

\\newpage
\\subsection{Continuous$\\sim$Discrete}
\\input{cd2}

\\newpage
\\subsection{Continuous$\\sim$Continuous}
\\input{cc4}




\\newpage
\\section{Too much...}
\\subsection{Factor $\\sim$ Logical}
\\input{fl2}

\\newpage
\\subsection{Factor $\\sim$ Factor}
\\input{ff5}

\\newpage
\\subsection{Ordered $\\sim$ Factor}
\\input{of4}

")

setwd("..")



