#####################################################################
#                         CIPW trace norm                           #
#                       (Verma, 2002, 2003)                         #
#####################################################################

CIPW.trace<-
function (filename, Volcanic = TRUE, AdjTAS = TRUE, Cancrinite = FALSE, Calcite = FALSE, digits = 3){
######################################################################
#                ##Environment for dataframes##                      #
######################################################################
Weight <- new.env()
data(AtomWeight, envir = Weight)
data(MinWeight, envir = Weight)
data(OxiWeight, envir = Weight)


#####################################################################
#                        Create dataframe                           #
#                                                                   #
#####################################################################
DATA<-data.frame(filename)

x.maj <- data.frame(DATA[,'SiO2'],DATA[,'TiO2'],DATA[,'Al2O3'],DATA[,'Fe2O3'],DATA[,'FeO'],DATA[,'MnO'],DATA[,'MgO'],DATA[,'CaO'],DATA[,'Na2O'],DATA[,'K2O'],DATA[,'P2O5'],DATA[,'CO2.'],DATA[,'LOI']) 
x.maj [is.na(x.maj)] <- 0 ##Change values NA->0
colx.maj<-c("sio2","tio2","al2o3","fe2o3","feo","mno","mgo","cao","na2o","k2o","p2o5","co2","loi")
colnames(x.maj)<-colx.maj

x.min <- data.frame(DATA[,'Cr2O3'],DATA[,'NiO'],DATA[,'SO3'],DATA[,'Ba'],DATA[,'Cl'],DATA[,'Co'],DATA[,'Cr'],DATA[,'Cs'],DATA[,'F.'],DATA[,'Li'],DATA[,'Ni'],DATA[,'Rb'],DATA[,'S'],DATA[,'Sr'],DATA[,'V'],DATA[,'Zr'])
x.min [is.na(x.min)] <- 0 ##Change values NA->0
colx.min<-c("cr2o3","nio","so3","ba","cl","co","cr","cs","f","li","ni","rb","s","sr","v","zr")
colnames(x.min)<-colx.min

x.min$cr2o3 <- ifelse(x.min$cr2o3 > 0, x.min$cr2o3, (x.min$cr*(Weight$OxiWeight["Cr2O3","OWeight"] / Weight$AtomWeight["Cr","AWeight"]))/10000)   
x.min$nio <- ifelse(x.min$nio > 0, x.min$nio, (x.min$nio*(Weight$OxiWeight["NiO","OWeight"] / Weight$AtomWeight["Ni","AWeight"]))/10000) 
x.min$bao <- (x.min$ba*(Weight$OxiWeight["BaO","OWeight"] / Weight$AtomWeight["Ba","AWeight"]))/10000
x.min$cl <- x.min$cl/10000
x.min$coo <- (x.min$co*(Weight$OxiWeight["CoO","OWeight"] / Weight$AtomWeight["Co","AWeight"]))/10000
x.min$cs2o <- (x.min$cs*(Weight$OxiWeight["Cs2O","OWeight"] / Weight$AtomWeight["Cs","AWeight"]))/10000
x.min$f <- x.min$f/10000
x.min$li2o <- (x.min$li*(Weight$OxiWeight["Li2O","OWeight"] / Weight$AtomWeight["Li","AWeight"]))/10000
x.min$rb2o <- (x.min$rb*(Weight$OxiWeight["Rb2O","OWeight"] / Weight$AtomWeight["Rb","AWeight"]))/10000
x.min$s <- x.min$s/10000
x.min$sro <- (x.min$sr*(Weight$OxiWeight["SrO","OWeight"] / Weight$AtomWeight["Sr","AWeight"]))/10000
x.min$v2o3 <- (x.min$v*(Weight$OxiWeight["V2O3","OWeight"] / Weight$AtomWeight["V","AWeight"]))/10000
x.min$zr2o3 <- (x.min$v*(Weight$OxiWeight["ZrO2","OWeight"] / Weight$AtomWeight["Zr","AWeight"]))/10000

######################################################################
#             Create dataframe and define loi and co2                #
######################################################################
x.tot <- data.frame(x.maj,x.min)
x <- x.tot[,-c(17,19:21,23:25,27:29)] #remove ppm
x$co2 <- ifelse (Cancrinite == TRUE | Calcite == TRUE, x$co2, 0)
x$loi <- x$loi-x$co2

#####################################################################
#                                                                   #
#                  Calculate Fe-oxidation ratio                     #
#####################################################################
fe <- AdjRock (filename, Volcanic = Volcanic, AdjTAS= AdjTAS, Cancrinite = Cancrinite, Calcite = Calcite, digits=digits)

######################################################################
#                Define Fe.output (Pruseth, 2009)                    #
######################################################################
T.Fe <- x$fe2o3 + (x$feo*1.111348)
x$feo <- round((T.Fe/(1+(fe$Fe.ratio *0.899809)))*0.899809,  digits)#
x$fe2o3 <- round(T.Fe-(T.Fe/(1+(fe$Fe.ratio *0.899809))), digits)#

######################################################################
#     ##Recalculate major elements data on an anhydrous basis##      #
######################################################################
######################################################################
x <- x[,-13] #remove loi
Sums.m <- rowSums(x)
adj<-function (x) round ((100*x)/(Sums.m), digits)##function of anhydrous adjust
x <- adj (x)#

######################################################################
#                          ##Changenames##                           #
######################################################################
col.x<-c("SiO2.adj","TiO2.adj","Al2O3.adj","Fe2O3.adj","FeO.adj","MnO.adj","MgO.adj","CaO.adj","Na2O.adj","K2O.adj","P2O5.adj","CO2.adj","Cr2O3.adj","NiO.adj","SO3.adj","Cl.adj","F.adj","S.adj","BaO.adj","CoO.adj","Cs2O.adj","Li2O.adj","Rb2O.adj","SrO.adj","V2O3.adj","Zr2O.adj")
colnames(x)<-col.x

######################################################################
#                     ##Mole computation##                           #
######################################################################
######################################################################
NORRRM.y<-c("sio2","tio2","al2o3","fe2o3","feo","mno","mgo","cao","na2o","k2o","p2o5","co2","cr2o3","so3","nio","coo","sro","bao","rb2o","cs2o","li2o","zr2o","v2o3","f","cl","s","cao.corr","ca.corr","na2o.corr","na.corr","k2o.corr","cr2o3.corr","xmno","xfeo","feo.corr","xnio","xcoo","xcao","xsro","xbao","xk2o","xrb2o","xcs2o","xna2o","xli2o","xcr2o3","xv2o3")
y<-rep(0,length(NORRRM.y))
names(y)<-NORRRM.y
y<-as.list(y)

y$sio2 <- (x$SiO2.adj/Weight$OxiWeight["SiO2","OWeight"])
y$tio2 <- (x$TiO2.adj/Weight$OxiWeight["TiO2","OWeight"])
y$al2o3 <- (x$Al2O3.adj/Weight$OxiWeight["Al2O3","OWeight"])
y$fe2o3 <- (x$Fe2O3.adj/Weight$OxiWeight["Fe2O3","OWeight"])
y$mgo <- (x$MgO.adj/Weight$OxiWeight["MgO","OWeight"])
y$p2o5 <- (x$P2O5.adj/Weight$OxiWeight["P2O5","OWeight"]) 
y$co2 <- (x$CO2.adj/Weight$OxiWeight["CO2","OWeight"])
y$zr2o <- (x$Zr2O.adj/Weight$OxiWeight["ZrO2","OWeight"])
y$so3 <- (x$SO3.adj/Weight$OxiWeight["SO3","OWeight"])
y$f <- (x$F.adj/Weight$AtomWeight["F","AWeight"])
y$cl <- (x$Cl.adj/Weight$AtomWeight["Cl","AWeight"])
y$s <- (x$S.adj/Weight$AtomWeight["S","AWeight"])

######################################################################
#              ##Corrected molecular weight FeO##                    #
######################################################################
y$mno <- (x$MnO.adj/Weight$OxiWeight["MnO","OWeight"]) 
y$nio <- (x$NiO.adj/Weight$OxiWeight["NiO","OWeight"])
y$coo <- (x$CoO.adj/Weight$OxiWeight["CoO","OWeight"])
y$feo <- ((x$FeO.adj/Weight$OxiWeight["FeO","OWeight"]) + y$mno + y$nio + y$coo)#test out (y$xmno + y$xfeo + y$xnio + y$xmno)== 1

y$xfeo <- ((x$FeO.adj/Weight$OxiWeight["FeO","OWeight"])/y$feo) 
y$xfeo[is.na(y$xfeo)] <- 0 ##Change values NA->0
y$xmno <- (y$mno/y$feo) 
y$xmno[is.na(y$xmno)] <- 0 ##Change values NA->0
y$xnio <- (y$nio/y$feo)
y$xnio[is.na(y$xnio)] <- 0 ##Change values NA->0 
y$xcoo <- (y$coo/y$feo) 
y$xcoo[is.na(y$xcoo)] <- 0 ##Change values NA->0
y$feo.corr <- (Weight$OxiWeight["FeO","OWeight"]*y$xfeo) + (Weight$OxiWeight["MnO","OWeight"]*y$xmno)+ (Weight$OxiWeight["NiO","OWeight"]*y$xnio)+ (Weight$OxiWeight["CoO","OWeight"]*y$xcoo)
y$fe.corr <- y$feo.corr - Weight$AtomWeight["O","AWeight"]

######################################################################
#              ##Corrected molecular weight CaO##                    #
######################################################################
y$sro <- (x$SrO.adj/Weight$OxiWeight["SrO","OWeight"])
y$bao <- (x$BaO.adj/Weight$OxiWeight["BaO","OWeight"])
y$cao <- (x$CaO.adj/Weight$OxiWeight["CaO","OWeight"] + y$sro + y$bao)

y$xcao <- ((x$CaO.adj/Weight$OxiWeight["CaO","OWeight"])/y$cao) 
y$xcao[is.na(y$xcao)] <- 0 ##Change values NA->0
y$xsro <- (y$sro/y$cao) 
y$xsro[is.na(y$xsro)] <- 0 ##Change values NA->0
y$xbao <- (y$bao/y$cao) 
y$xbao[is.na(y$xbao)] <- 0 ##Change values NA->0
y$cao.corr <- (Weight$OxiWeight["CaO","OWeight"]*y$xcao) + (Weight$OxiWeight["SrO","OWeight"]*y$xsro) + (Weight$OxiWeight["BaO","OWeight"]*y$xbao)
y$ca.corr <- y$cao.corr - Weight$AtomWeight["O","AWeight"]

######################################################################
#              ##Corrected molecular weight K2O##                    #
######################################################################
y$rb2o <- (x$Rb2O.adj/Weight$OxiWeight["Rb2O","OWeight"])
y$cs2o <- (x$Cs2O.adj/Weight$OxiWeight["Cs2O","OWeight"])
y$k2o <- (x$K2O.adj/Weight$OxiWeight["K2O","OWeight"] +y$rb2o +y$cs2o) 

y$xk2o <- ((x$K2O.adj/Weight$OxiWeight["K2O","OWeight"])/y$k2o) 
y$xk2o[is.na(y$xk2o)] <- 0 ##Change values NA->0
y$xrb2o <- (y$rb2o/y$k2o) 
y$xrb2o[is.na(y$xrb2o)] <- 0 ##Change values NA->0
y$xcs2o <- (y$cs2o/y$k2o) 
y$xcs2o[is.na(y$xcs2o)] <- 0 ##Change values NA->0
y$k2o.corr <- (Weight$OxiWeight["K2O","OWeight"]*y$xk2o) + (Weight$OxiWeight["Rb2O","OWeight"]*y$xrb2o) + (Weight$OxiWeight["Cs2O","OWeight"]*y$xcs2o)

######################################################################
#              ##Corrected molecular weight Na2O##                   #
######################################################################
y$li2o <- (x$Li2O.adj/Weight$OxiWeight["Li2O","OWeight"])
y$na2o <- (x$Na2O.adj/Weight$OxiWeight["Na2O","OWeight"] + y$li2o) 

y$xna2o <- ((x$Na2O.adj/Weight$OxiWeight["Na2O","OWeight"])/y$na2o) 
y$xna2o[is.na(y$xna2o)] <- 0 ##Change values NA->0
y$xli2o <- (y$li2o/y$na2o) 
y$xli2o[is.na(y$xli2o)] <- 0 ##Change values NA->0
y$na2o.corr <- (Weight$OxiWeight["Na2O","OWeight"]*y$xna2o) + (Weight$OxiWeight["Li2O","OWeight"]*y$xli2o) 
y$na.corr <- (y$na2o.corr - Weight$AtomWeight["O","AWeight"])/2

######################################################################
#             ##Corrected molecular weight Cr2O3##                   #
######################################################################
y$v2o3 <- (x$V2O3.adj/Weight$OxiWeight["V2O3","OWeight"])
y$cr2o3 <- (x$Cr2O3.adj/Weight$OxiWeight["Cr2O3","OWeight"] + y$v2o3)
y$xcr2o3 <- ((x$Cr2O3/Weight$OxiWeight["Cr2O3","OWeight"])/y$cr2o3) 
y$xcr2o3[is.na(y$xcr2o3)] <- 0 ##Change values NA->0
y$xv2o3 <- (y$v2o3/y$cr2o3) 
y$xv2o3[is.na(y$xv2o3)] <- 0 ##Change values NA->0
y$cr2o3.corr <- (Weight$OxiWeight["Cr2O3","OWeight"]*y$xcr2o3) + (Weight$OxiWeight["V2O3","OWeight"]*y$xv2o3) 

######################################################################
#                         Main functions                             #
######################################################################
NORRRM.minerals.mole<-c("Q","C","Or","Ab","An","Ne","Lc","Kp","Nc","FREE_Cl","FREE_SO3","FREE_Cr2O3","f","FREE_CO2","Ac","Z","Ns","Ks","Di","Wo","Hy","Ol","Ds","Th","Cm","Cc","Mt","Il","Hm","Tn","Pf","Ru","Y","Ap","DEFSIO2","FREE_P2O5","Fr","Pr","FREEO_12b","FREEO_12c","FREEO_13","FREEO_14","FREEO_15","Orp","Abp","Lcp","Dip","Olp","Wop","D","D1","D2","D3","D4","D5","D6","Tnp","Hl","ApCaO","FREE_F","ApCaF")
m<-rep(0,length(NORRRM.minerals.mole))
names(m)<-NORRRM.minerals.mole
m<-as.list(m)

###############################################################
#                          Zircon                             #
###############################################################
m$Z <- y$zr2o
m$Y<- m$Z

###############################################################
#                    Apatite/FREE_P2O5                        #
###############################################################
m$Ap <- ifelse (y$cao >= (10 / 3)*y$p2o5, y$p2o5, 3/10 * y$cao)
y$cao <- ifelse (y$cao >= ((10 / 3)*y$p2o5), y$cao-(m$Ap * 10 / 3), 0)
y$cao <- ifelse (y$p2o5 > 0, y$cao, (x$CaO.adj/Weight$OxiWeight["CaO","OWeight"] + y$sro + y$bao))
m$FREE_P2O5 <- ifelse (y$cao == 0, y$p2o5-m$Ap, 0)
m$FREE_P2O5 <- ifelse (y$p2o5 > 0, m$FREE_P2O5, 0)
m$f <- ifelse (y$f >= ((2 / 3)*m$Ap), ((2 / 3)*m$Ap), y$f)
y$f <- ifelse (y$f >= ((2 / 3)*m$Ap), y$f - ((2 / 3)*m$Ap), 0)
m$ApCaF <- ifelse (m$f > 0, 1.5*m$f, 0)
m$ApCaO <- m$Ap - m$ApCaF #test out (m$ApCaF + m$ApCaO)== m$Ap

##############################################################
m$FREEO_12b <- ifelse ((x$F.adj/Weight$AtomWeight["F","AWeight"]) >= (2 / 3)*m$Ap, (m$Ap * 1/3),0)
m$FREEO_12c <- ifelse ((x$F.adj/Weight$AtomWeight["F","AWeight"]) >= (2 / 3)*m$Ap, 0,(x$F.adj/Weight$AtomWeight["F","AWeight"])*0.5)
###############################################################
#                        Fluorite                             #
###############################################################
m$Fr <- ifelse (y$cao >= 0.5*y$f, 0.5*y$f, y$cao) 
y$cao <- ifelse (y$cao >= 0.5*y$f, y$cao - (0.5*y$f), 0)
m$FREE_F <- ifelse (y$cao >= 0.5*y$f, 0, y$f - (m$Fr*2))
m$FREEO_13 <- m$Fr
      
###############################################################
#                         Halite                              #
###############################################################
m$Hl <- ifelse (y$na2o >= 2*y$cl, y$cl, y$na2o)
m$FREEO_14 <- m$Hl*0.5
m$FREE_Cl <- ifelse (y$na2o >= 2*y$cl, 0,y$cl-(m$FREEO_14*0.25))
m$FREE_Cl <- ifelse (y$na2o > 0, m$FREE_Cl,y$cl)
y$na2o <- ifelse (y$na2o >= 2*y$cl, y$na2o-(m$Hl*0.5), 0)

###############################################################
#                       Thenardite                            #
###############################################################
m$Th <- ifelse (y$na2o >= y$so3, y$so3, y$na2o)
y$na2o <- ifelse (y$na2o >= y$so3, y$na2o - y$so3, y$na2o-m$Th)
m$FREE_SO3 <- ifelse (y$na2o >= y$so3, 0, y$so3-m$Th)

###############################################################
#                         Pyrite                              #
###############################################################
m$Pr <- ifelse (y$feo >= (0.5*y$s), y$s*0.5, y$feo)
m$Pr <- ifelse (y$s>0, m$Pr, 0)
y$feo <- ifelse (y$feo >= (0.5*y$s), y$feo - (y$s*0.5), 0)
m$FREE_S <- ifelse (((x$FeO.adj/Weight$OxiWeight["FeO","OWeight"]) + y$mno + y$nio + y$coo) >= (0.5*y$s), 0, y$s-(m$Pr*2))
m$FREE_S <- ifelse (y$s>0, m$FREE_S, 0)
y$feo <- ifelse (y$s>0, y$feo, ((x$FeO.adj/Weight$OxiWeight["FeO","OWeight"]) + y$mno + y$nio + y$coo))
m$FREEO_15 <- m$Pr

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
#                        Chromite                             #
###############################################################
m$Cm <- ifelse (y$feo >= y$cr2o3, y$cr2o3, y$feo)
y$feo <- ifelse (y$feo >= y$cr2o3, y$feo-y$cr2o3, 0)
m$FREE_Cr2O3 <- ifelse (y$feo >= y$cr2o3, 0, y$cr2o3 - m$Cm)
                                                                                                                       
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
m$Y <- ifelse (m$Ks > 0, m$Y+(6*m$Orp)+m$Ks, m$Y+(6*m$Orp))

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
m$Ds <- ifelse(m$D4 < m$Wop/2, m$D4, m$Wop/2)
m$Ds <- ifelse(unsaturated, m$Ds, 0)
m$Wopp <- ifelse(m$D4 < m$Wop/2, m$Wop-(m$D4*2), 0)
m$Wo <- ifelse(unsaturated, m$Wopp, m$Wop)
m$D5 <- ifelse(m$D4 >= m$Wop/2, m$D4-(m$Wop/2), 0)
unsaturated <- ifelse(m$D5 > 0,TRUE,FALSE)#

###############################################################
#         Adjust Diopside/Dicalcium silicate/Olivine          #
############################################################### 
m$Ds <- ifelse(m$D5 < m$Dip, m$Ds+(m$D5/2), m$Ds+m$Dip/2)
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
NORRRM.minerals.t <- c("Q","C","Or","Ab","An","Ne","Lc","Kp","Z","Ac","Ns","Ks","Di.Mg","Di.Fe","Wo","En","Fs","Fo","Fa","Dcs","Mt","Il","Hm","Tn","Pf","Ru","Ap-F","Ap-Ca","Fr","NaCl","Th","Pr","Nc","Cc","Cm","DEFSIO2","FREE_P2O5","FREE_F","FREE_Cl","FREE_SO3","FREE_S","FREE_CO2","FREE_Cr2O3","FREEO","Sums.min")
Normative.minerals.t <- matrix(data = NA, ncol = length(NORRRM.minerals.t), nrow = nrow(x), byrow = FALSE, dimnames = NULL)
colnames(Normative.minerals.t) <- NORRRM.minerals.t

Normative.minerals.t[,"Z"] <- m$Z*Weight$MinWeight["Zircon",3]
Normative.minerals.t[,"C"] <- m$C*Weight$MinWeight["Corundum",3]
Normative.minerals.t[,"Q"] <- m$Q*Weight$MinWeight["Quartz",3]
Normative.minerals.t[,"En"] <- (m$Hy*y$xmgr) * Weight$MinWeight["Hypersthene-Mg",3]
Normative.minerals.t[,"Hm"] <- m$Hm*Weight$MinWeight["Hematite",3]
Normative.minerals.t[,"Ru"] <- m$Ru*Weight$MinWeight["Rutile",3]
Normative.minerals.t[,"An"] <- m$An*(y$cao.corr + Weight$MinWeight["Anorthite",4])
Normative.minerals.t[,"Di.Mg"] <- (m$Di*y$xmgr) * (y$cao.corr + Weight$MinWeight["Diopside-Mg",4])
Normative.minerals.t[,"Wo"] <- m$Wo*(y$cao.corr + Weight$MinWeight["Wollastonite",4])
Normative.minerals.t[,"Dcs"] <- m$Ds*((y$cao.corr*2) + Weight$MinWeight["Dicalcium silicate",4])
Normative.minerals.t[,"Tn"] <- m$Tn*(y$cao.corr + Weight$MinWeight["Sphene",4])
Normative.minerals.t[,"Pf"] <- m$Pf*(y$cao.corr + Weight$MinWeight["Perovskite",4])
Normative.minerals.t[,"Ap-F"] <- m$ApCaF*((3*y$cao.corr) + ((1/3)*(y$ca.corr)) + Weight$MinWeight["Apatite-F",4])
Normative.minerals.t[,"Ap-Ca"] <- m$ApCaO*((10/3)*(y$cao.corr) + (Weight$MinWeight["Apatite-Ca",4]))
Normative.minerals.t[,"Cc"] <- m$Cc*(y$cao.corr + Weight$MinWeight["Calcite",4])
Normative.minerals.t[,"Fs"] <- (m$Hy*y$xfer) * (y$feo.corr + Weight$MinWeight["Hypersthene-Fe",4]) 
Normative.minerals.t[,"Fa"] <- ((m$Ol*y$xfer) * ((y$feo.corr*2) + Weight$MinWeight["Olivine-Fe",4]))#
Normative.minerals.t[,"Fo"] <- ((m$Ol*y$xmgr) * Weight$MinWeight["Olivine-Mg",3])
Normative.minerals.t[,"Mt"] <- m$Mt*(y$feo.corr + Weight$MinWeight["Magnetite",4])
Normative.minerals.t[,"Il"] <- m$Il*(y$feo.corr + Weight$MinWeight["Ilmenite",4])
Normative.minerals.t[,"Ab"] <- m$Ab*(y$na2o.corr + Weight$MinWeight["Albite",4])
Normative.minerals.t[,"Ne"] <- m$Ne*(y$na2o.corr + Weight$MinWeight["Nepheline",4])
Normative.minerals.t[,"Th"] <- m$Th*(y$na2o.corr + Weight$MinWeight["Thenardite",4])
Normative.minerals.t[,"Nc"] <- m$Nc*(y$na2o.corr + Weight$MinWeight["Sodium Carbonate",4])
Normative.minerals.t[,"Ac"] <- m$Ac*(y$na2o.corr + Weight$MinWeight["Acmite",4])
Normative.minerals.t[,"Ns"] <- m$Ns*(y$na2o.corr + Weight$MinWeight["Sodium metasilicate",4])
Normative.minerals.t[,"Or"] <- m$Or*(y$k2o.corr + Weight$MinWeight["Orthoclase",4])
Normative.minerals.t[,"Lc"] <- m$Lc*(y$k2o.corr + Weight$MinWeight["Leucite",4])
Normative.minerals.t[,"Kp"] <- m$Kp*(y$k2o.corr + Weight$MinWeight["Kaliophilite",4])
Normative.minerals.t[,"Ks"] <- m$Ks*(y$k2o.corr + Weight$MinWeight["Potassium metasilicate",4])
Normative.minerals.t[,"Di.Fe"] <- (m$Di*y$xfer) * (y$cao.corr + y$feo.corr + Weight$MinWeight["Diopside-Ca",4]) 
Normative.minerals.t[,"Cm"] <- m$Cm*(y$feo.corr + y$cr2o3.corr)
Normative.minerals.t[,"NaCl"] <- m$Hl*(y$na.corr + Weight$AtomWeight["Cl","AWeight"])
Normative.minerals.t[,"Fr"] <- m$Fr*(y$ca.corr + Weight$MinWeight["Fluorite",4])
Normative.minerals.t[,"Pr"] <- m$Pr*(y$fe.corr + (Weight$AtomWeight["S","AWeight"]*2))

######################################################################
#                             Excess                                 #
######################################################################
Normative.minerals.t[,"FREE_F"] <- m$FREE_F * Weight$AtomWeight["F","AWeight"]
Normative.minerals.t[,"FREE_Cl"] <- m$FREE_Cl  * Weight$AtomWeight["Cl","AWeight"]
Normative.minerals.t[,"FREE_SO3"] <- m$FREE_SO3 * Weight$OxiWeight["SO3","OWeight"]
Normative.minerals.t[,"FREE_S"] <- m$FREE_S * Weight$AtomWeight["S","AWeight"]
Normative.minerals.t[,"FREE_CO2"] <- m$FREE_CO2 * Weight$OxiWeight["CO2","OWeight"]
Normative.minerals.t[,"FREE_Cr2O3"] <- m$FREE_Cr2O3 * y$cr2o3.corr 
Normative.minerals.t[,"FREE_P2O5"] <- m$FREE_P2O5*Weight$OxiWeight["P2O5",1]
######################################################################
#                             Free O                                 #
######################################################################
NORRRM.FreeO<-c("MWApCF2","FreeO_a","FreeO_b","FreeO_c","FreeO_d","FreeO_e")
f<-rep(0,length(NORRRM.FreeO))
names(f)<-NORRRM.FreeO
f<-as.list(f)
f$MWApCF2<-((3*y$cao.corr) + ((1/3)*(y$ca.corr)) + Weight$MinWeight["Apatite-F",4])
f$FreeO_a <- (1 + ((0.1) *((f$MWApCF2/328.8691887)-1)))* Weight$AtomWeight["O","AWeight"] * m$FREEO_12b
f$FreeO_b <- (1 + ((0.1) * (m$ApCaF/m$Ap) * ((f$MWApCF2/Weight$MinWeight["Apatite-F",3])-1))) * Weight$AtomWeight["O","AWeight"] *  m$FREEO_12c
f$FreeO_b[is.na(f$FreeO_b)] <- 0 ##Change values NA->0
f$FreeO_c <- (1 + ((y$cao.corr/56.0774)-1)) * Weight$AtomWeight["O","AWeight"] *  m$FREEO_13
f$FreeO_d <- (1 + (y$na2o.corr/Weight$OxiWeight['Na2O','OWeight'])-1) * Weight$AtomWeight["O","AWeight"] *  m$FREEO_14
f$FreeO_e <- (1 + ((y$feo.corr/Weight$OxiWeight['FeO','OWeight'])-1)) * Weight$AtomWeight["O","AWeight"] *  m$FREEO_15
Normative.minerals.t[,"FREEO"] <- f$FreeO_a + f$FreeO_b + f$FreeO_c + f$FreeO_d + f$FreeO_e
Normative.minerals.t[,"DEFSIO2"] <- m$DEFSIO2*Weight$OxiWeight["SiO2",2] * (-1)
output.minerals <- data.frame (Normative.minerals.t,  stringsAsFactors = FALSE)
output.minerals[,"Sums.min"] <- round(rowSums(output.minerals, na.rm = any(!is.na(output.minerals))),3)


######################################################################
#                         Print Parameters                           #
######################################################################
NORRRM.parameters<-c("Salic","Femic","C.I.","D.I.","S.I.","A.R.","Mg.number","Fe.index","MALI","ACNK","ANK","AI","FSSI","Density.rock")
P<-rep(0,length(NORRRM.parameters))
names(P)<-NORRRM.parameters
P<-as.list(P)

P$Salic <- Normative.minerals.t[,"Q"] + Normative.minerals.t[,"Or"] + Normative.minerals.t[,"Ab"] + Normative.minerals.t[,"An"]
P$Femic <- Normative.minerals.t[,"Di.Mg"] + Normative.minerals.t[,"Di.Fe"] + Normative.minerals.t[,"En"] + Normative.minerals.t[,"Fs"] + Normative.minerals.t[,"Fo"] + Normative.minerals.t[,"Fa"] + Normative.minerals.t[,"Mt"] + Normative.minerals.t[,"Il"] + Normative.minerals.t[,"Hm"]
P$C.I. <- Normative.minerals.t[,"An"] + (2.1570577*Normative.minerals.t[,"Di.Mg"]) + Normative.minerals.t[,"Fo"] + (0.7007616*Normative.minerals.t[,"Fs"])
P$D.I. <-  Normative.minerals.t[,"Q"] + Normative.minerals.t[,"Or"] + Normative.minerals.t[,"Ab"] + Normative.minerals.t[,"Ne"] + Normative.minerals.t[,"Lc"] 
P$S.I. <- (100*x$MgO.adj)/(x$MgO.adj + x$Fe2O3.adj + x$FeO.adj + x$Na2O.adj + x$K2O.adj)
P$A.R. <- ifelse (x$SiO2.adj > 0 & x$K2O.adj/x$Na2O.adj > 1.0 & x$K2O.adj/x$Na2O.adj < 2.5, (x$Al2O3.adj + x$CaO.adj + (2*x$Na2O.adj))/(x$Al2O3.adj + x$CaO.adj - (2*x$Na2O.adj)), (x$Al2O3.adj + x$CaO.adj + x$Na2O.adj + x$K2O.adj)/(x$Al2O3.adj + x$CaO.adj - x$Na2O.adj - x$K2O.adj) )
P$Mg.number <- 100 * ((Weight$AtomWeight["Mg",2]/Weight$OxiWeight["MgO",2]) * x$MgO.adj)/(((Weight$AtomWeight["Mg",2]/Weight$OxiWeight["MgO",2]) * x$MgO.adj) + ((Weight$AtomWeight["Fe",2]/Weight$OxiWeight["FeO",2]) * x$FeO.adj))
P$Fe.index <- x$FeO.adj / (x$FeO.adj + x$MgO.adj)
P$MALI <-  x$Na2O.adj + x$K2O.adj - x$CaO.adj
P$ACNK <- ((Weight$AtomWeight["Al",2]/Weight$OxiWeight["Al2O3",2]) * x$Al2O3.adj) / ((((Weight$AtomWeight["Ca",2]/Weight$OxiWeight["CaO",2]) * x$CaO.adj) - (1.67 * ((Weight$AtomWeight["P",2]/Weight$OxiWeight["P2O5",2]) * x$P2O5.adj))) * ((Weight$AtomWeight["Na",2]/Weight$OxiWeight["Na2O",2]) * x$Na2O.adj) * ((Weight$AtomWeight["K",2]/Weight$OxiWeight["K2O",2]) * x$K2O.adj))
P$ANK <- ((Weight$AtomWeight["Al",2]/Weight$OxiWeight["Al2O3",2]) * x$Al2O3.adj) / (((Weight$AtomWeight["Na",2]/Weight$OxiWeight["Na2O",2]) * x$Na2O.adj) * ((Weight$AtomWeight["K",2]/Weight$OxiWeight["K2O",2]) * x$K2O.adj))
P$AI <- ((Weight$AtomWeight["Al",2]/Weight$OxiWeight["Al2O3",2]) * x$Al2O3.adj) - (((Weight$AtomWeight["K",2]/Weight$OxiWeight["K2O",2]) * x$K2O.adj) - ((Weight$AtomWeight["Na",2]/Weight$OxiWeight["Na2O",2]) * x$Na2O.adj))
P$FSSI <- Normative.minerals.t[,"Q"] - (Normative.minerals.t[,"Lc"]  + (2*(Normative.minerals.t[,"Ne"]  + Normative.minerals.t[,"Kp"] ))) /100

#####################################################################
#                                                                   #
#            Calculate the theoretical density of rocks             #
#####################################################################

DensityRock<-
function (Normative.minerals.t, na.rm = TRUE,normsum=FALSE){

######################################################################
#                     ##Matrix computation##                         #
######################################################################

######################################################################
#                      Vol. prop. of mineral                         #
######################################################################
NORRRM.Vol<- c("Z","Q","C","Or","Ab","An","Ne","Lc","Kp","Ac","Ns","Ks","Di.Mg","Di.Fe","Wo","En","Fs","Fo","Fa","Cs","Mt","Il","Hm","Tn","Pf","Ru","Ap-F","Ap-Ca","Nc","Cc","Cm","Hl","Fr","Pr")
Volume <- matrix(data = NA, ncol = length(NORRRM.Vol), nrow = nrow(Normative.minerals.t), byrow = FALSE, dimnames = NULL)
colnames(Volume) <- NORRRM.Vol
Volume[,"Z"]  <- Normative.minerals.t[,"Z"] / Weight$MinWeight ["Zircon",5] 
Volume[,"Q"]  <- Normative.minerals.t[,"Q"] / Weight$MinWeight ["Quartz",5]     
Volume[,"C"] <- Normative.minerals.t[,"C"] / Weight$MinWeight["Corundum",5]
Volume[,"Or"] <- Normative.minerals.t[,"Or"] / Weight$MinWeight["Orthoclase",5]
Volume[,"Ab"] <- Normative.minerals.t[,"Ab"] / Weight$MinWeight["Albite",5]
Volume[,"An"] <- Normative.minerals.t[,"An"] / Weight$MinWeight["Anorthite",5]
Volume[,"Ne"] <- Normative.minerals.t[,"Ne"] / Weight$MinWeight["Nepheline",5]
Volume[,"Lc"] <- Normative.minerals.t[,"Lc"] / Weight$MinWeight["Leucite",5]
Volume[,"Kp"] <- Normative.minerals.t[,"Kp"] / Weight$MinWeight["Kaliophilite",5]
Volume[,"Ac"] <- Normative.minerals.t[,"Ac"] / Weight$MinWeight["Acmite",5]
Volume[,"Ns"] <- Normative.minerals.t[,"Ns"] / Weight$MinWeight["Sodium metasilicate",5]
Volume[,"Ks"] <- Normative.minerals.t[,"Ks"] / Weight$MinWeight["Potassium metasilicate",5]
Volume[,"Ap-F"] <- Normative.minerals.t[,"Ap-F"] / Weight$MinWeight["Apatite-F",5]
Volume[,"Ap-Ca"] <- Normative.minerals.t[,"Ap-Ca"] / Weight$MinWeight["Apatite-Ca",5]
Volume[,"Di.Mg"] <- Normative.minerals.t[,"Di.Mg"] / Weight$MinWeight["Diopside-Mg",5]
Volume[,"Di.Fe"] <- Normative.minerals.t[,"Di.Fe"] / Weight$MinWeight["Diopside-Fe",5]
Volume[,"Wo"] <- Normative.minerals.t[,"Wo"] / Weight$MinWeight["Wollastonite",5]
Volume[,"En"] <- Normative.minerals.t[,"En"] / Weight$MinWeight["Hypersthene-Mg",5]
Volume[,"Fs"]  <- Normative.minerals.t[,"Fs"] /  Weight$MinWeight["Hypersthene-Fe",5]
Volume[,"Fo"]  <- Normative.minerals.t[,"Fo"] /  Weight$MinWeight["Olivine-Mg",5]
Volume[,"Fa"]  <- Normative.minerals.t[,"Fa"]  / Weight$MinWeight["Olivine-Fe",5]#
Volume[,"Cs"]  <- Normative.minerals.t[,"Dcs"] / Weight$MinWeight["Dicalcium silicate",5]
Volume[,"Mt"] <- Normative.minerals.t[,"Mt"] / Weight$MinWeight["Magnetite",5]
Volume[,"Il"] <- Normative.minerals.t[,"Il"] / Weight$MinWeight["Ilmenite",5]
Volume[,"Hm"] <- Normative.minerals.t[,"Hm"] / Weight$MinWeight["Hematite",5]
Volume[,"Tn"] <- Normative.minerals.t[,"Tn"] / Weight$MinWeight["Sphene",5]
Volume[,"Pf"] <- Normative.minerals.t[,"Pf"] / Weight$MinWeight["Perovskite",5]
Volume[,"Ru"] <- Normative.minerals.t[,"Ru"] / Weight$MinWeight["Rutile",5]
Volume[,"Nc"] <- Normative.minerals.t[,"Nc"] / Weight$MinWeight["Sodium Carbonate",5]
Volume[,"Cc"] <- Normative.minerals.t[,"Cc"] / Weight$MinWeight["Calcite",5]
Volume[,"Cm"] <- Normative.minerals.t[,"Cm"] / Weight$MinWeight["Chromite",5]
Volume[,"Hl"] <- Normative.minerals.t[,"NaCl"] / Weight$MinWeight["Halite",5]
Volume[,"Fr"] <- Normative.minerals.t[,"Fr"] / Weight$MinWeight["Fluorite",5]
Volume[,"Pr"] <- Normative.minerals.t[,"Pr"] / Weight$MinWeight["Pyrite",5]


Sums.Density<- rowSums(Volume, na.rm = any(!is.na(Volume)))


######################################################################
#                      Density of minerals                           #
######################################################################
NORRRM.Density<- c("Z","Q","C","Or","Ab","An","Ne","Lc","Kp","Ac","Ns","Ks","Di.Mg","Di.Fe","Wo","En","Fs","Fo","Fa","Cs","Mt","Il","Hm","Tn","Pf","Ru","Ap-F","Ap-Ca","Nc","Cc","Cm","Hl","Fr","Pr")
Dens <- matrix(data = NA, ncol = length(NORRRM.Density), nrow = nrow(Normative.minerals.t), byrow = FALSE, dimnames = NULL)
colnames(Dens) <- NORRRM.Density
Dens[,"Z"] <- ((Volume[,"Z"] * (100/Sums.Density)) * Weight$MinWeight ["Zircon",5]) / 100  
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
Dens[,"Ap-F"] <- ((Volume[,"Ap-F"] * (100/Sums.Density)) * Weight$MinWeight["Apatite-F",5]) / 100
Dens[,"Ap-Ca"] <- ((Volume[,"Ap-Ca"] * (100/Sums.Density)) * Weight$MinWeight["Apatite-Ca",5]) / 100
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
Dens[,"Cm"] <- ((Volume[,"Cm"] * (100/Sums.Density)) * Weight$MinWeight["Chromite",5]) / 100
Dens[,"Hl"] <- ((Volume[,"Hl"] * (100/Sums.Density)) * Weight$MinWeight["Halite",5]) / 100
Dens[,"Fr"] <- ((Volume[,"Fr"] * (100/Sums.Density)) * Weight$MinWeight["Fluorite",5]) / 100
Dens[,"Pr"] <- ((Volume[,"Pr"] * (100/Sums.Density)) * Weight$MinWeight["Pyrite",5]) / 100
Density.Rock  <- round ( rowSums(Dens, na.rm = any(!is.na(Dens))), digits = digits)###Total Density
}

## End(Not run)

P$Density.rock <- DensityRock (Normative.minerals.t)
######################################################################
#                         Clean Environment                          #
######################################################################
rm(Weight)


######################################################################
#                         Define Output                              #
######################################################################

##Define output##
fe[,'Total.adj'] <- rowSums(x)
fe[,'Total.trace'] <- x$Cr2O3.adj+x$NiO.adj+x$SO3.adj+x$Cl.adj+x$F.adj+x$S.adj+x$BaO.adj+x$CoO.adj+x$Cs2O.adj+x$Li2O.adj+x$Rb2O.adj+x$SrO.adj+x$V2O3.adj+x$Zr2O.adj
output.minerals <- round (output.minerals, digits)
output.parameters <- round (data.frame (P,  stringsAsFactors = FALSE), digits )
output <- cbind(fe,"",output.minerals,"",output.parameters)#
row.names(output) <- row.names(fe)
output [output == 0] <- NA  ##Change values 0->NA  

return(output)

}
## End(Not run)
