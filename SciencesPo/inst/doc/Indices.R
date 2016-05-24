## ----highestAverages1, message=FALSE, comment=NA-------------------------
library("SciencesPo")

# The d'Hondt will give the same results as Jefferson's method
highestAverages(parties=names(lijphart),
                votes=lijphart,
                seats = 6, 
                method = "dh") 

## ----highestAverages2, echo=TRUE, message=FALSE, comment=NA--------------
# The Sainte-Laguë will give the same results as the Webster's method (wb)
highestAverages(parties=names(lijphart),
                votes=lijphart,
                seats = 6,
                method = "sl") 

## ----highestAverages3, echo=TRUE, message=FALSE, comment=NA--------------
highestAverages(parties=names(lijphart),
                votes=lijphart, 
                seats = 6,
                method = "msl") 

## ----highestAverages4, echo=TRUE, message=FALSE, comment=NA--------------
highestAverages(parties=names(lijphart), 
                votes=lijphart, 
                seats = 6, 
                method = "danish") 

## ----highestAverages5, echo=TRUE, message=FALSE, comment=NA--------------
highestAverages(parties=names(lijphart),
                votes=lijphart, 
                seats = 6, 
                method = "hsl") 

## ----highestAverages6, echo=TRUE, message=FALSE, comment=NA--------------
highestAverages(parties=names(lijphart),
                votes=lijphart, 
                seats = 6, 
                method = "wb") 

## ----highestAverages7, echo=TRUE, message=FALSE, comment=NA--------------
highestAverages(parties=names(lijphart),
                votes=lijphart, 
                seats = 6, 
                method = "imperiali") 

## ----highestAverages8, echo=TRUE, message=FALSE, comment=NA--------------
Bruges=c("CD&V/N-VA"=32092, "SP.A/Spirit"=20028, 
         "Flemish Interest"=13408, "Open VLD/Vivant"=9520,
         "Green!"=5328, "Other"=2207)

highestAverages(parties=names(Bruges),
                votes=Bruges, 
                seats = 47,
                method = "imperiali") 

## ----highestAverages9, echo=TRUE, message=FALSE, comment=NA--------------
highestAverages(parties=names(lijphart),
                votes=lijphart, 
                seats = 6, method = "hh") 

## ----highestAverages10, echo=TRUE, message=FALSE, comment=NA-------------

const <- c("A"=100, "B"=150,"C"=300, "D"=400, "E"=50)

highestAverages(parties=names(const),
                votes=const,
                seats = 3, method = "dh",
                threshold = 7/100) 

## ----Valencia, echo=TRUE, message=FALSE, comment=NA----------------------
# Valencia returned 15 members
highestAverages(parties=names(Valencia),
                votes=Valencia, 
                seats=15, method = "dh",
                threshold = 3/100)

## ----largestRemainders1, eval=FALSE, echo=TRUE, message=FALSE, comment=NA----
#  largestRemainders(parties=names(lijphart),
#                    votes=lijphart,
#                    seats = 8, method = "hare")

## ----largestRemainders2, eval=FALSE, echo=TRUE, message=FALSE, comment=NA----
#  largestRemainders(parties=names(lijphart),
#                    votes=lijphart,
#                    seats = 8, method = "droop")

## ----largestRemainders3, eval=FALSE, echo=TRUE, message=FALSE, comment=NA----
#  largestRemainders(parties=names(lijphart),
#                    votes=lijphart,
#                    seats = 8, method = "hagb")

## ----data-Italy, eval=FALSE, echo=TRUE, message=FALSE--------------------
#  # The 1946 Italian Constituent Assembly election results: parties and unspoilt votes
#  
#  Italy = data.frame(party=c("DC", "PSIUP", "PCI", "UDN", "UQ", "PRI",
#                              "BNL", "PdA", "MIS", "PCd'I", "CDR",
#                             "PSd'Az", "MUI", "PCS", "PDL", "FDPR"),
#                     votes=c(8101004, 4758129, 4356686, 1560638,	1211956,
#                             1003007, 637328, 334748, 171201, 102393,
#                             97690, 78554, 71021, 51088, 40633, 21853))

## ----largestRemainders4, eval=FALSE, echo=TRUE, message=FALSE, comment=NA----
#  with(Italy, largestRemainders(parties=party,
#                                votes=votes, seats = 556,
#                                method = "imperiali.q") )

## ----highestAverages11, echo=TRUE, message=FALSE, comment=NA-------------
mytable = highestAverages(parties=names(Ceara), 
                          votes=Ceara,
                          seats = 42, method = "dh") 

library(knitr)

kable(mytable, align=c("l","c","c"), caption="Outcome under d'Hondt")

## ----largestRemainders5, eval=TRUE, echo=TRUE, message=FALSE, fig.width=7, fig.height=4.5, fig.align="center", fig.cap= "2014 Legislative Election in Ceará (M=42)"----

out1 = highestAverages(parties=names(Ceara), votes=Ceara, 
                seats = 42, method = "dh")
out2 = highestAverages(parties=names(Ceara), votes=Ceara, 
                seats = 42, method = "imperiali") 
out3 = highestAverages(parties=names(Ceara), votes=Ceara, 
                seats = 42, method = "sl")

# add the method:
out1$Method = "d'Hondt"
out2$Method = "imperiali"
out3$Method = "Saint-Laguë"


data <- rbind(out1, out2, out3)

p = ggplot(data=data, aes(x=reorder(Party, -Seats), y=Seats, fill=Method)) +
    geom_bar(stat="identity",position=position_dodge()) +
   labs(x="", y="Seats")
p + scale_fill_fte() + 
  theme_fte(legend = "top") 

## ----largestRemainders6, eval=FALSE, echo=TRUE, message=FALSE, fig.width=7, fig.height=4.5, fig.align="center", fig.cap= "2014 Legislative Election in Ceara (M=42)"----
#  
#  #2014 Federal elections, 30 seats to be returned in the state of Parana, Brazil.
#  
#  PR=c("PSDB/DEM/PR/PSC/PTdoB/PP/SD/PSD/PPS"=2601709,
#       "PT/PDT/PRB/PTN/PCdoB"=1109905,
#       "PSDC/PEN/PTB/PHS/PMN/PROS"=501148,
#       "PV/PPL"=280767)
#  
#  2014 Federal elections, 70 seats to be returned in the state of Sao Paulo, Brazil.
#  SP=c("PSDB/DEM/PPS"=5537630, "PT/PCdoB"=3170003,
#       "PMDB/PROS/PP/PSD"=2384740, "PSOL/PSTU"=462992,
#       "PSL/PTN/PMN/PTC/PTdoB"=350186, "PHS/PRP"=252205)
#  
#  

## ----politicalDiversity1, echo=TRUE, message=FALSE, comment=NA-----------
# The 2004 presidential election in the US (vote share):

US2004 <- c("Democratic"=0.481, "Republican"=0.509, 
            "Independent"=0.0038, "Libertarian"=0.0032, 
            "Constitution"=0.0012, "Green"=0.00096,
            "Others"=0.00084)

print(US2004)

## ----politicalDiversity2, echo=TRUE, message=FALSE, comment=NA-----------
politicalDiversity(US2004); # ENEP (laakso/taagepera) method 

## ----politicalDiversity3, echo=TRUE, message=FALSE, comment=NA-----------
politicalDiversity(US2004, index= "golosov");

## ----politicalDiversity4, echo=TRUE, message=FALSE, comment=NA-----------
politicalDiversity(US2004, index= "herfindahl");

## ----Helsinki-election, echo=TRUE, message=FALSE-------------------------
# Helsinki's 1999

Helsinki <- data.frame(
  votes = c(68885, 18343, 86448, 21982, 51587,
            27227, 8482, 7250, 365, 2734, 1925,
            475, 1693, 693, 308, 980, 560, 590, 185),
  seats.SL=c(5, 1, 6, 1, 4, 2, 1, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0),
  seats.dH=c(5, 1, 7, 1, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0))

## ----politicalDiversity5, echo=TRUE, message=FALSE, comment=NA-----------
# politicalDiversity(Helsinki$votes); #ENEP Votes

politicalDiversity(Helsinki$seats.SL); #ENP for Saint-Lague

politicalDiversity(Helsinki$seats.dH); #ENP for D'Hondt

## ----proportionality1, echo=TRUE, message=FALSE, comment=NA--------------
with(Queensland, gallagher(pvotes, pseats))

with(Quebec, gallagher(pvotes, pseats))

## ----proportionality2, echo=TRUE, message=FALSE, comment=NA--------------
with(Queensland, lijphart(pvotes, pseats))

with(Quebec, lijphart(pvotes, pseats))

## ----proportionality3, echo=TRUE, message=FALSE, comment=NA--------------
with(Queensland, grofman(pvotes, pseats))

with(Quebec, grofman(pvotes, pseats))

## ----proportionality4, echo=TRUE, message=FALSE, comment=NA--------------
with(Queensland, farina(pvotes, pseats))

with(Quebec, farina(pvotes, pseats))

## ----proportionality5, echo=TRUE, message=FALSE, comment=NA--------------
with(Queensland, cox.shugart(pvotes, pseats))

with(Quebec, cox.shugart(pvotes, pseats))

## ----proportionality6, echo=TRUE, message=FALSE, comment=NA--------------
with(Queensland, inv.cox.shugart(pvotes, pseats))

with(Quebec, inv.cox.shugart(pvotes, pseats))

## ----eval=FALSE, echo=FALSE, message=FALSE, comment=NA-------------------
#  sessionInfo()

