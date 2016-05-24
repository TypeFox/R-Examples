#####################################################################
#                           CIPW norm                               #
#                       (Verma, 2002, 2003)                         #
#####################################################################

CIPW<-
function (filename, Volcanic = TRUE, AdjTAS = TRUE, Cancrinite = FALSE, Calcite = FALSE, digits = 3){
#####################################################################
#                                                                   #
#       Calculate major elements data on an anhydrous basis         #
#####################################################################
x <- AdjRock (filename, Volcanic = Volcanic, AdjTAS= AdjTAS, Cancrinite = Cancrinite, Calcite = Calcite, digits=digits)
x[is.na(x)] <- 0 ##Change values NA->0

######################################################################
#                ##Environment for dataframes##                      #
######################################################################

Weight <- new.env()
data(AtomWeight, envir = Weight)
data(MinWeight, envir = Weight)
data(OxiWeight, envir = Weight)

######################################################################
#                     ##Mole computation##                           #
######################################################################
######################################################################
NORRRM.y<-c("sio2","tio2","al2o3","fe2o3","feo","mno","mgo","cao","na2o","k2o","p2o5","xmno","xfeo")
y<-rep(0,length(NORRRM.y))
names(y)<-NORRRM.y
y<-as.list(y)

y$sio2 <- (x$SiO2.adj/Weight$OxiWeight["SiO2","Rounded"])
y$tio2 <- (x$TiO2.adj/Weight$OxiWeight["TiO2","Rounded"])
y$al2o3 <- (x$Al2O3.adj/Weight$OxiWeight["Al2O3","Rounded"])
y$mno <- (x$MnO.adj/Weight$OxiWeight["MnO","Rounded"]) 
y$feo <- ((x$FeO.adj/Weight$OxiWeight["FeO","Rounded"]) + (x$MnO.adj/Weight$OxiWeight["MnO","Rounded"]))
y$fe2o3 <- (x$Fe2O3.adj/Weight$OxiWeight["Fe2O3","Rounded"])
y$mgo <- (x$MgO.adj/Weight$OxiWeight["MgO","Rounded"])
y$cao <- (x$CaO.adj/Weight$OxiWeight["CaO","Rounded"])
y$na2o <- (x$Na2O.adj/Weight$OxiWeight["Na2O","Rounded"]) 
y$k2o <- (x$K2O.adj/Weight$OxiWeight["K2O","Rounded"]) 
y$p2o5 <- (x$P2O5.adj/Weight$OxiWeight["P2O5","Rounded"]) 
y$xmno <- (y$mno/y$feo) 
y$xfeo <- ((x$FeO.adj/Weight$OxiWeight["FeO","Rounded"])/y$feo) 
#test out <- round (y$xmno + y$xfeo, digits = 3)== 1
y$co2 <- (x$CO2.adj/Weight$OxiWeight["CO2","Rounded"]) 

######################################################################
#                         Main functions                             #
######################################################################
NORRRM.minerals.mole<-c("Q","C","Or","Ab","An","Ne","Lc","Kp","Ac","Ns","Ks","Di","Wo","Hy","Ol","Cs","Csp","Mt","Il","Hm","Tn","Pf","Ru","Y","Ap","Nc","Cc","DEFSIO2","FREE_P2O5","FREE_CO2","Orp","Abp","Lcp","Dip","Olp","Olpp","Wop","Hyp","Hypp","D","D1","D2","D3","D4","D5","D6","D7","Tnp")
m<-rep(0,length(NORRRM.minerals.mole))
names(m)<-NORRRM.minerals.mole
m<-as.list(m)

###############################################################
#                    Apatite/FREE_P2O5                        #
###############################################################
m$Ap <- ifelse (y$cao >= (10 / 3)*y$p2o5, y$p2o5, 3/10 * y$cao)
y$cao <- ifelse (y$cao >= (10 / 3)*y$p2o5, y$cao-(m$Ap * 10 / 3), 0)
y$p2o5 <- ifelse (y$cao == 0, y$p2o5-m$Ap, 0)
m$FREE_P2O5 <- ifelse (y$p2o5 > 0, y$p2o5, 0)

###############################################################
#               Sodium carbonate/calcite                      #
###############################################################
Cancrinite [is.na(Cancrinite)] <- FALSE ##Change values NA->FALSE
m$Nc <- ifelse (Cancrinite == TRUE & y$na2o >= y$co2, y$co2, 0)
m$Nc <- ifelse (Cancrinite == TRUE & y$na2o < y$co2, y$na2o, m$Nc)
y$na2o <- ifelse (Cancrinite == TRUE & y$na2o >= y$co2, y$na2o-y$co2, y$na2o)
y$na2o <- ifelse (Cancrinite == TRUE & y$na2o < y$co2, 0, y$na2o)
m$FREE_CO2 <- ifelse (Cancrinite == TRUE & y$na2o < y$co2, y$co2-m$Nc, 0)

###############################################################
Calcite [is.na(Calcite)] <- FALSE ##Change values NA->FALSE
m$Cc <- ifelse (Calcite== TRUE & y$cao >= y$co2, y$co2, 0)
m$Cc <- ifelse (Calcite == TRUE & y$cao < y$co2, y$cao, m$Cc)
y$cao <- ifelse (Calcite == TRUE & y$cao >= y$co2, y$cao-m$Cc, y$cao)
y$cao <- ifelse (Calcite == TRUE & y$cao < y$co2, 0, y$cao)
y$co2 <- ifelse (Calcite == TRUE & y$cao < y$co2, y$co2-m$Cc, 0)
m$FREE_CO2 <- ifelse (Calcite == TRUE & y$cao < y$co2, y$co2-m$Cc, m$FREE_CO2)

###############################################################
#                        Ilmenite                             #
###############################################################
m$Il <- ifelse (y$feo >= y$tio2, y$tio2, y$feo)
y$feo <- ifelse (y$feo >= y$tio2, y$feo-m$Il, 0)
y$tio2 <- ifelse (y$feo ==0, y$tio2-m$Il, 0)

###############################################################
#             Orthoclase/Potassium metasilicate               #
###############################################################       
m$Orp <- ifelse (y$al2o3 >= y$k2o, y$k2o, y$al2o3)
y$al2o3 <- ifelse (y$al2o3 >= y$k2o, y$al2o3-m$Orp, 0)
y$k2o <- ifelse (y$al2o3 >= y$k2o, 0, y$k2o-m$Orp)
m$Ks <- ifelse (y$k2o > 0, y$k2o, 0)
m$Y <- ifelse (m$Ks > 0, (6*m$Orp)+m$Ks, 6*m$Orp)

###############################################################
#                       Albite                                #
############################################################### 
m$Abp <- ifelse (y$al2o3 >= y$na2o, y$na2o, y$al2o3)
y$al2o3 <- ifelse (y$al2o3 >= y$na2o, y$al2o3-m$Abp, 0)
y$na2o <- ifelse (y$al2o3 == 0, y$na2o-m$Abp, 0)
m$Y<- m$Y+(6*m$Abp)

###############################################################
#                Acmite/Sodium metasilicate                   #
############################################################### 
m$Ac <- ifelse (y$na2o >= y$fe2o3, y$fe2o3, y$na2o)
y$na2o <- ifelse (y$na2o >= y$fe2o3, y$na2o-m$Ac, 0)
y$fe2o3 <- ifelse (y$na2o == 0, y$fe2o3-m$Ac, 0)
m$Ns <- ifelse (y$na2o > 0, y$na2o, 0)
m$Y <- ifelse (m$Ns > 0, (m$Ac*4)+m$Ns+m$Y, (m$Ac*4)+m$Y)

###############################################################
#                    Anorthite/Corundum                       #
############################################################### 
m$An <- ifelse (y$al2o3 >= y$cao, y$cao, y$al2o3)
y$al2o3 <- ifelse (y$al2o3 >= y$cao, y$al2o3-m$An, 0)
y$cao <- ifelse (y$al2o3 == 0, y$cao-m$An, 0)
m$C <- ifelse (y$al2o3 > 0 , y$al2o3, 0)
m$Y<- (m$An*2)+m$Y

###############################################################
#                      Sphene/Rutile                          #
############################################################### 
m$Tnp <- ifelse(y$cao >= y$tio2, y$tio2, y$cao)
y$cao <- ifelse(y$cao >= y$tio2, y$cao-m$Tnp, 0)
y$tio2 <- ifelse(y$cao == 0, y$tio2-m$Tnp, 0)
m$Ru <- ifelse(y$tio2 > 0 , y$tio2, 0)
m$Y<- (m$Tnp)+m$Y

###############################################################
#                  Magnetite/Hematite                         #
############################################################### 
m$Mt <- ifelse(y$fe2o3 >= y$feo, y$feo, y$fe2o3)
y$feo <- ifelse(y$fe2o3 >= y$feo, 0, y$feo-m$Mt)
y$fe2o3 <- ifelse(y$feo == 0, y$fe2o3-m$Mt, 0)
m$Hm <- ifelse(y$fe2o3 > 0, y$fe2o3, 0)


######################################################################
#       Subdivision of some normative minerals (end-members)         #
######################################################################
y$femg <- y$feo + y$mgo

y$xfer <- y$feo / (y$feo + y$mgo)
y$xmgr <- y$mgo / (y$mgo + y$feo)

###############################################################
#             Diopside/Wollastonite/Hypersthene               #
############################################################### 
m$Dip <- ifelse(y$cao >= y$femg, y$femg, y$cao)
y$cao <- ifelse(y$cao >= y$femg, y$cao-m$Dip, 0)
y$femg <- ifelse(y$cao  == 0, y$femg-m$Dip, 0)
m$Wop <- ifelse(y$cao > 0, y$cao, 0)
m$Hyp <- ifelse(y$femg > 0, y$femg, 0)
m$Y <- ifelse(m$Wop > 0 ,(m$Dip*2)+m$Wop+m$Y, (m$Dip*2)+m$Hyp+m$Y)

###############################################################
#                Quartz/Undersatured Quartz                   #
############################################################### 
m$Q <- ifelse(y$sio2 >= m$Y, y$sio2-m$Y, 0)
m$D <- ifelse(y$sio2 < m$Y, m$Y-y$sio2, 0)
saturated <- ifelse(y$sio2 >= m$Y,TRUE,FALSE)#
unsaturated <- ifelse(y$sio2 < m$Y,TRUE,FALSE)

###############################################################
#                   Olivine/Hypersthene                       #
############################################################### 
m$Olp <- ifelse(m$D < (m$Hyp/2), m$D, m$Hyp/2)
m$Olp <- ifelse(unsaturated, m$Olp, 0)
m$Hy <- ifelse(m$D < (m$Hyp/2), m$Hyp-(2*m$D), m$Hy<-0)
m$Hy <- ifelse(unsaturated, m$Hy, m$Hyp)
m$D1 <- ifelse(m$D >= (m$Hyp/2), m$D-(m$Hyp/2), 0)
m$D1 <- ifelse(unsaturated, m$D1, 0)
unsaturated <- ifelse(m$D1 > 0,TRUE,FALSE)#

###############################################################
#                    Sphene/Perovskite                        #
############################################################### 
m$Pf <- ifelse(m$D1 < m$Tnp, m$D1, m$Tnp)
m$Pf <- ifelse(unsaturated, m$Pf, 0)
m$Tn <- ifelse(m$D1 < m$Tnp, m$Tnp-m$D1, 0)
m$D2 <- ifelse(m$D1 >= m$Tnp, m$D1-m$Tnp, 0)
unsaturated <- ifelse(m$D2 > 0,TRUE,FALSE)#

###############################################################
#                    Nepheline/Albite                         #
############################################################### 
m$Ne <- ifelse(m$D2 < m$Abp*4, m$D2/4, m$Abp)
m$Ne <- ifelse(unsaturated, m$Ne, 0)
m$Abpp <- ifelse(m$D2 < m$Abp*4, m$Abp-(m$D2/4), 0)
m$Ab <- ifelse(unsaturated, m$Abpp, m$Abp)
m$D3 <- ifelse(m$D2 >= m$Abp*4, m$D2-(m$Abp*4), 0)
unsaturated <- ifelse(m$D3 > 0,TRUE,FALSE)#

###############################################################
#                    Leucite/Orthoclase                       #
############################################################### 
m$Lcp <- ifelse(m$D3 < m$Orp*2, m$D3/2, m$Orp)
m$Lcp <- ifelse(unsaturated, m$Lcp, 0)
m$Orpp <- ifelse(m$D3 < m$Orp*2, m$Orp-(m$D3/2), 0)
m$Or <- ifelse(unsaturated, m$Orpp, m$Orp)
m$D4 <- ifelse(m$D3 >= m$Orp*2, m$D3-(m$Orp*2), 0)
unsaturated <- ifelse(m$D4 > 0,TRUE,FALSE)#
                    
###############################################################
#            Dicalcium silicate/Wollastonite                  #
###############################################################     
m$Cs <- ifelse(m$D4 < m$Wop/2, m$D4, m$Wop/2)
m$Cs <- ifelse(unsaturated, m$Cs, 0)
m$Wopp <- ifelse(m$D4 < m$Wop/2, m$Wop-(m$D4*2), 0)
m$Wo <- ifelse(unsaturated, m$Wopp, m$Wop)
m$D5 <- ifelse(m$D4 >= m$Wop/2, m$D4-(m$Wop/2), 0)
unsaturated <- ifelse(m$D5 > 0,TRUE,FALSE)#

###############################################################
#         Adjust Diopside/Dicalcium silicate/Olivine          #
############################################################### 
m$Cs <- ifelse(m$D5 < m$Dip, m$Cs+(m$D5/2), m$Cs+m$Dip/2)
m$Olpp <- ifelse(m$D5 < m$Dip, m$Olp+(m$D5/2), m$Olp+m$Dip/2)
m$Ol <- ifelse(unsaturated, m$Olpp, m$Olp)
m$Dipp <- ifelse(m$D5 < m$Dip, m$Dip-(m$D5), 0)
m$Di <- ifelse(unsaturated, m$Dipp, m$Dip)
m$D6 <- ifelse(m$D5 >= m$Dip, m$D5-(m$Dip), 0)
unsaturated <- ifelse(m$D6 > 0,TRUE,FALSE)#

###############################################################
#                Adjust Kaliophilite/Leucite                  #
############################################################### 
m$Kp <- ifelse(m$Lcp >= m$D6/2, m$D6/2, m$Lcp)
m$Lcpp <- ifelse(m$Lcp >= m$D6/2, m$Lcp-(m$D6/2), 0)
m$Lc <- ifelse(unsaturated, m$Lcpp, m$Lcp)
m$DEFSIO2  <- m$D6-(2*m$Kp)

######################################################################
#                         Print Minerals                             #
######################################################################
NORRRM.minerals <- c("Q","C","Or","Ab","An","Ne","Lc","Kp","Ac","Ns","Ks","Di.Mg","Di.Fe","Wo","En","Fs","Fo","Fa","Cs","Mt","Il","Hm","Tn","Pf","Ru","Ap","Nc","Cc","FREE_CO2","DEFSIO2","FREE_P2O5","Sums.min")
Normative.minerals <- matrix(data = NA, ncol = length(NORRRM.minerals), nrow = nrow(x), byrow = FALSE, dimnames = NULL)
colnames(Normative.minerals) <- NORRRM.minerals

Normative.minerals[,"Q"] <- (m$Q*Weight$MinWeight["Quartz",3])
Normative.minerals[,"C"] <- (m$C*Weight$MinWeight["Corundum",3])
Normative.minerals[,"Or"] <- (m$Or*Weight$MinWeight["Orthoclase",3])
Normative.minerals[,"Ab"] <- (m$Ab*Weight$MinWeight["Albite",3])
Normative.minerals[,"An"] <- (m$An*Weight$MinWeight["Anorthite",3])
Normative.minerals[,"Ne"] <- (m$Ne*Weight$MinWeight["Nepheline",3])
Normative.minerals[,"Lc"] <- (m$Lc*Weight$MinWeight["Leucite",3])
Normative.minerals[,"Kp"] <- (m$Kp*Weight$MinWeight["Kaliophilite",3])
Normative.minerals[,"Ac"] <- (m$Ac*Weight$MinWeight["Acmite",3])
Normative.minerals[,"Ns"] <- (m$Ns*Weight$MinWeight["Sodium metasilicate",3])
Normative.minerals[,"Ks"] <- (m$Ks*Weight$MinWeight["Potassium metasilicate",3])
Normative.minerals[,"Di.Mg"] <- ((m$Di*y$xmgr) * Weight$MinWeight["Diopside-Mg",3])
Normative.minerals[,"Di.Fe"] <- ((m$Di*y$xfer) * ((Weight$OxiWeight["FeO","OWeight"]*y$xfeo) + (Weight$OxiWeight["MnO","OWeight"]*y$xmno) + Weight$MinWeight["Diopside-Fe",4]))
Normative.minerals[,"Wo"] <- (m$Wo*Weight$MinWeight["Wollastonite",3])
Normative.minerals[,"En"] <- ((m$Hy*y$xmgr) * Weight$MinWeight["Hypersthene-Mg",3])
Normative.minerals[,"Fs"] <- ((m$Hy*y$xfer) * ((Weight$OxiWeight["FeO","OWeight"]*y$xfeo) + (Weight$OxiWeight["MnO","OWeight"]*y$xmno) + Weight$MinWeight["Hypersthene-Fe",4]))
Normative.minerals[,"Fo"] <- ((m$Ol*y$xmgr) * Weight$MinWeight["Olivine-Mg",3])
Normative.minerals[,"Fa"] <- ((m$Ol*y$xfer) * ((((Weight$OxiWeight["FeO","OWeight"])*y$xfeo*2) + ((Weight$OxiWeight["MnO","OWeight"]*y$xmno*2))) + Weight$MinWeight["Olivine-Fe",4]))
Normative.minerals[,"Cs"] <- (m$Cs*Weight$MinWeight["Dicalcium silicate",3])
Normative.minerals[,"Mt"] <- (m$Mt*((Weight$OxiWeight["FeO","OWeight"]*y$xfeo) + (Weight$OxiWeight["MnO","OWeight"]*y$xmno) + Weight$MinWeight["Magnetite",4]))
Normative.minerals[,"Il"] <- (m$Il*((Weight$OxiWeight["FeO","OWeight"]*y$xfeo) + (Weight$OxiWeight["MnO","OWeight"]*y$xmno) + Weight$MinWeight["Ilmenite",4]))
Normative.minerals[,"Hm"] <- (m$Hm*Weight$MinWeight["Hematite",3])
Normative.minerals[,"Tn"] <- (m$Tn*Weight$MinWeight["Sphene",3])
Normative.minerals[,"Pf"] <- (m$Pf*Weight$MinWeight["Perovskite",3])
Normative.minerals[,"Ru"] <- (m$Ru*Weight$MinWeight["Rutile",3])
Normative.minerals[,"Ap"] <- (m$Ap*Weight$MinWeight["Apatite-Ca",3])
Normative.minerals[,"Nc"] <- (m$Nc*Weight$MinWeight["Sodium Carbonate",3])
Normative.minerals[,"Cc"] <- (m$Cc*Weight$MinWeight["Calcite",3])
Normative.minerals[,"DEFSIO2"] <- m$DEFSIO2*Weight$OxiWeight["SiO2",2] * (-1)
Normative.minerals[,"FREE_P2O5"] <- (m$FREE_P2O5*Weight$OxiWeight["P2O5",2])
Normative.minerals[,"FREE_CO2"] <- (m$FREE_CO2*Weight$OxiWeight["CO2",2])
output.minerals <- round (data.frame (Normative.minerals,  stringsAsFactors = FALSE), digits = digits)
output.minerals[,"Sums.min"] <- round (rowSums(output.minerals, na.rm = any(!is.na(output.minerals))), digits = digits)


######################################################################
#                         Print Parameters                           #
######################################################################
NORRRM.parameters<-c("Salic","Femic","C.I.","D.I.","S.I.","A.R.","Mg.number","Fe.index","MALI","ACNK","ANK","AI","FSSI","Density.rock")
P<-rep(0,length(NORRRM.parameters))
names(P)<-NORRRM.parameters
P<-as.list(P)

P$Salic <- round (Normative.minerals[,"Q"] + Normative.minerals[,"Or"] + Normative.minerals[,"Ab"] + Normative.minerals[,"An"], digits = digits)
P$Femic <- Normative.minerals[,"Di.Mg"] + Normative.minerals[,"Di.Fe"] + Normative.minerals[,"En"] + Normative.minerals[,"Fs"] + Normative.minerals[,"Fo"] + Normative.minerals[,"Fa"] + Normative.minerals[,"Mt"] + Normative.minerals[,"Il"] + Normative.minerals[,"Hm"]
P$C.I. <- Normative.minerals[,"An"] + (2.1570577*Normative.minerals[,"Di.Mg"]) + Normative.minerals[,"Fo"] + (0.7007616*Normative.minerals[,"Fs"])
P$D.I. <-  Normative.minerals[,"Q"] + Normative.minerals[,"Or"] + Normative.minerals[,"Ab"] + Normative.minerals[,"Ne"] + Normative.minerals[,"Lc"] 
P$S.I. <- (100*x$MgO.adj)/(x$MgO.adj + x$Fe2O3.adj + x$FeO.adj + x$Na2O.adj + x$K2O.adj)
P$A.R. <- ifelse (x$SiO2.adj > 0 & x$K2O.adj/x$Na2O.adj > 1.0 & x$K2O.adj/x$Na2O.adj < 2.5, (x$Al2O3.adj + x$CaO.adj + (2*x$Na2O.adj))/(x$Al2O3.adj + x$CaO.adj - (2*x$Na2O.adj)), (x$Al2O3.adj + x$CaO.adj + x$Na2O.adj + x$K2O.adj)/(x$Al2O3.adj + x$CaO.adj - x$Na2O.adj - x$K2O.adj) )
P$Mg.number <- 100 * ((Weight$AtomWeight["Mg",2]/Weight$OxiWeight["MgO",2]) * x$MgO.adj)/(((Weight$AtomWeight["Mg",2]/Weight$OxiWeight["MgO",2]) * x$MgO.adj) + ((Weight$AtomWeight["Fe",2]/Weight$OxiWeight["FeO",2]) * x$FeO.adj))
P$Fe.index <- x$FeO.adj / (x$FeO.adj + x$MgO.adj)
P$MALI <-  round(x$Na2O.adj + x$K2O.adj - x$CaO.adj,digits=digits)
P$ACNK <- ((Weight$AtomWeight["Al",2]/Weight$OxiWeight["Al2O3",2]) * x$Al2O3.adj) / ((((Weight$AtomWeight["Ca",2]/Weight$OxiWeight["CaO",2]) * x$CaO.adj) - (1.67 * ((Weight$AtomWeight["P",2]/Weight$OxiWeight["P2O5",2]) * x$P2O5.adj))) * ((Weight$AtomWeight["Na",2]/Weight$OxiWeight["Na2O",2]) * x$Na2O.adj) * ((Weight$AtomWeight["K",2]/Weight$OxiWeight["K2O",2]) * x$K2O.adj))
P$ANK <- ((Weight$AtomWeight["Al",2]/Weight$OxiWeight["Al2O3",2]) * x$Al2O3.adj) / (((Weight$AtomWeight["Na",2]/Weight$OxiWeight["Na2O",2]) * x$Na2O.adj) * ((Weight$AtomWeight["K",2]/Weight$OxiWeight["K2O",2]) * x$K2O.adj))
P$AI <- ((Weight$AtomWeight["Al",2]/Weight$OxiWeight["Al2O3",2]) * x$Al2O3.adj) - (((Weight$AtomWeight["K",2]/Weight$OxiWeight["K2O",2]) * x$K2O.adj) - ((Weight$AtomWeight["Na",2]/Weight$OxiWeight["Na2O",2]) * x$Na2O.adj))
P$FSSI <- Normative.minerals[,"Q"] - (Normative.minerals[,"Lc"]  + (2*(Normative.minerals[,"Ne"]  + Normative.minerals[,"Kp"] ))) /100


#####################################################################
#                                                                   #
#            Calculate the theoretical density of rocks             #
#####################################################################

DensityRock<-
function (Normative.minerals, na.rm = TRUE,normsum=FALSE){

######################################################################
#                     ##Matrix computation##                         #
######################################################################

######################################################################
#                      Vol. prop. of mineral                         #
######################################################################
NORRRM.Vol<- c("Q","C","Or","Ab","An","Ne","Lc","Kp","Ac","Ns","Ks","Di.Mg","Di.Fe","Wo","En","Fs","Fo","Fa","Cs","Mt","Il","Hm","Tn","Pf","Ru","Ap","Nc","Cc")
Volume <- matrix(data = NA, ncol = length(NORRRM.Vol), nrow = nrow(Normative.minerals), byrow = FALSE, dimnames = NULL)
colnames(Volume) <- NORRRM.Vol

Volume[,"Q"]  <- round ( Normative.minerals[,"Q"] / Weight$MinWeight ["Quartz",5], digits = digits)         
Volume[,"C"] <- Normative.minerals[,"C"] / Weight$MinWeight["Corundum",5]
Volume[,"Or"] <- Normative.minerals[,"Or"] / Weight$MinWeight["Orthoclase",5]
Volume[,"Ab"] <- Normative.minerals[,"Ab"] / Weight$MinWeight["Albite",5]
Volume[,"An"] <- Normative.minerals[,"An"] / Weight$MinWeight["Anorthite",5]
Volume[,"Ne"] <- Normative.minerals[,"Ne"] / Weight$MinWeight["Nepheline",5]
Volume[,"Lc"] <- Normative.minerals[,"Lc"] / Weight$MinWeight["Leucite",5]
Volume[,"Kp"] <- Normative.minerals[,"Kp"] / Weight$MinWeight["Kaliophilite",5]
Volume[,"Ac"] <- Normative.minerals[,"Ac"] / Weight$MinWeight["Acmite",5]
Volume[,"Ns"] <- Normative.minerals[,"Ns"] / Weight$MinWeight["Sodium metasilicate",5]
Volume[,"Ks"] <- Normative.minerals[,"Ks"] / Weight$MinWeight["Potassium metasilicate",5]
Volume[,"Ap"] <- Normative.minerals[,"Ap"] / Weight$MinWeight["Apatite-Ca",5]
Volume[,"Di.Mg"] <- Normative.minerals[,"Di.Mg"] / Weight$MinWeight["Diopside-Mg",5]
Volume[,"Di.Fe"] <- Normative.minerals[,"Di.Fe"] / Weight$MinWeight["Diopside-Fe",5]
Volume[,"Wo"] <- Normative.minerals[,"Wo"] / Weight$MinWeight["Wollastonite",5]
Volume[,"En"] <- Normative.minerals[,"En"] / Weight$MinWeight["Hypersthene-Mg",5]
Volume[,"Fs"]  <- Normative.minerals[,"Fs"] /  Weight$MinWeight["Hypersthene-Fe",5]
Volume[,"Fo"]  <- Normative.minerals[,"Fo"] /  Weight$MinWeight["Olivine-Mg",5]
Volume[,"Fa"]  <- Normative.minerals[,"Fa"]  / Weight$MinWeight["Olivine-Fe",5]#
Volume[,"Cs"]  <- Normative.minerals[,"Cs"] / Weight$MinWeight["Dicalcium silicate",5]
Volume[,"Mt"] <- Normative.minerals[,"Mt"] / Weight$MinWeight["Magnetite",5]
Volume[,"Il"] <- Normative.minerals[,"Il"] / Weight$MinWeight["Ilmenite",5]
Volume[,"Hm"] <- Normative.minerals[,"Hm"] / Weight$MinWeight["Hematite",5]
Volume[,"Tn"] <- Normative.minerals[,"Tn"] / Weight$MinWeight["Sphene",5]
Volume[,"Pf"] <- Normative.minerals[,"Pf"] / Weight$MinWeight["Perovskite",5]
Volume[,"Ru"] <- Normative.minerals[,"Ru"] / Weight$MinWeight["Rutile",5]
Volume[,"Nc"] <- Normative.minerals[,"Nc"] / Weight$MinWeight["Sodium Carbonate",5]
Volume[,"Cc"] <- Normative.minerals[,"Cc"] / Weight$MinWeight["Calcite",5]

Sums.Density<- rowSums(Volume, na.rm = any(!is.na(Volume)))


######################################################################
#                      Density of minerals                           #
######################################################################
NORRRM.Density<- c("Q","C","Or","Ab","An","Ne","Lc","Kp","Ac","Ns","Ks","Di.Mg","Di.Fe","Wo","En","Fs","Fo","Fa","Cs","Mt","Il","Hm","Tn","Pf","Ru","Ap","Nc","Cc")
Dens <- matrix(data = NA, ncol = length(NORRRM.Density), nrow = nrow(Normative.minerals), byrow = FALSE, dimnames = NULL)
colnames(Dens) <- NORRRM.Density

Dens[,"Q"] <- ((Volume[,"Q"] * (100/Sums.Density)) * Weight$MinWeight ["Quartz",5]) / 100  
Dens[,"C"] <- ((Volume[,"C"] * (100/Sums.Density)) * Weight$MinWeight ["Corundum",5]) / 100  
Dens[,"Or"] <- ((Volume[,"Or"] * (100/Sums.Density)) * Weight$MinWeight["Orthoclase",5]) / 100
Dens[,"Ab"] <- ((Volume[,"Ab"] * (100/Sums.Density)) * Weight$MinWeight["Albite",5]) / 100
Dens[,"An"] <- ((Volume[,"An"] * (100/Sums.Density)) * Weight$MinWeight["Anorthite",5]) / 100
Dens[,"Ne"] <- ((Volume[,"Ne"] * (100/Sums.Density)) * Weight$MinWeight["Nepheline",5]) / 100
Dens[,"Lc"] <- ((Volume[,"Lc"] * (100/Sums.Density)) * Weight$MinWeight["Leucite",5]) / 100
Dens[,"Kp"] <- ((Volume[,"Kp"] * (100/Sums.Density)) * Weight$MinWeight["Kaliophilite",5]) / 100
Dens[,"Ac"] <- ((Volume[,"Ac"] * (100/Sums.Density)) * Weight$MinWeight["Acmite",5]) / 100
Dens[,"Ns"] <- ((Volume[,"Ns"] * (100/Sums.Density)) * Weight$MinWeight["Sodium metasilicate",5]) / 100
Dens[,"Ks"] <- ((Volume[,"Ks"] * (100/Sums.Density)) * Weight$MinWeight["Potassium metasilicate",5]) / 100
Dens[,"Ap"] <- ((Volume[,"Ap"] * (100/Sums.Density)) * Weight$MinWeight["Apatite-Ca",5]) / 100
Dens[,"Di.Mg"] <- ((Volume[,"Di.Mg"] * (100/Sums.Density)) * Weight$MinWeight["Diopside-Mg",5]) / 100
Dens[,"Di.Fe"] <- ((Volume[,"Di.Fe"] * (100/Sums.Density)) * Weight$MinWeight["Diopside-Fe",5]) / 100
Dens[,"Wo"] <- ((Volume[,"Wo"] * (100/Sums.Density)) * Weight$MinWeight["Wollastonite",5]) / 100
Dens[,"En"] <- ((Volume[,"En"] * (100/Sums.Density)) * Weight$MinWeight["Hypersthene-Mg",5]) / 100
Dens[,"Fs"] <- ((Volume[,"Fs"] * (100/Sums.Density)) * Weight$MinWeight["Hypersthene-Fe",5]) / 100
Dens[,"Fo"] <- ((Volume[,"Fo"] * (100/Sums.Density)) * Weight$MinWeight["Olivine-Mg",5]) / 100
Dens[,"Fa"] <- ((Volume[,"Fa"] * (100/Sums.Density)) * Weight$MinWeight["Olivine-Fe",5]) / 100
Dens[,"Cs"] <- ((Volume[,"Cs"] * (100/Sums.Density)) * Weight$MinWeight["Dicalcium silicate",5]) / 100
Dens[,"Mt"] <- ((Volume[,"Mt"] * (100/Sums.Density)) * Weight$MinWeight["Magnetite",5]) / 100
Dens[,"Il"] <- ((Volume[,"Il"] * (100/Sums.Density)) * Weight$MinWeight["Ilmenite",5]) / 100
Dens[,"Hm"] <- ((Volume[,"Hm"] * (100/Sums.Density)) * Weight$MinWeight["Hematite",5]) / 100
Dens[,"Tn"]<- ((Volume[,"Tn"] * (100/Sums.Density)) * Weight$MinWeight["Sphene",5]) / 100
Dens[,"Pf"] <- ((Volume[,"Pf"] * (100/Sums.Density)) * Weight$MinWeight["Perovskite",5]) / 100
Dens[,"Ru"] <- ((Volume[,"Ru"] * (100/Sums.Density)) * Weight$MinWeight["Rutile",5]) / 100
Dens[,"Nc"] <- ((Volume[,"Nc"] * (100/Sums.Density)) * Weight$MinWeight["Sodium Carbonate",5]) / 100
Dens[,"Cc"] <- ((Volume[,"Cc"] * (100/Sums.Density)) * Weight$MinWeight["Calcite",5]) / 100

Density.Rock  <- round ( rowSums(Dens, na.rm = any(!is.na(Dens))), digits = digits)###Total Density

}

## End(Not run)

P$Density.rock <- DensityRock (Normative.minerals)

######################################################################
#                         Clean Environment                          #
######################################################################
rm(Weight)

######################################################################
#                         Define Output                              #
######################################################################

##Define output##
output.major <- x
output.parameters <- round (data.frame (P,  stringsAsFactors = FALSE), digits = digits)
output <- cbind(output.major,"",output.minerals,"",output.parameters)#
row.names(output) <- row.names(x)
output [output == 0] <- NA  ##Change values 0->NA  
return(output)

}
## End(Not run)
