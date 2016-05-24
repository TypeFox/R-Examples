#####################################################################
#                        Main function                              #
#       Calculate major elements data on an anhydrous basis         #
#####################################################################
# must recast to 100 %

AdjRock <- function(filename, Volcanic= TRUE, AdjTAS= TRUE, Cancrinite= FALSE, Calcite= FALSE, digits=3){

#####################################################################
#                        Create dataframe                           #
#                                                                   #
#####################################################################
DATA<-data.frame(filename)

x <- data.frame(DATA[,'SiO2'],DATA[,'TiO2'],DATA[,'Al2O3'],DATA[,'Fe2O3'],DATA[,'FeO'],DATA[,'MnO'],DATA[,'MgO'],DATA[,'CaO'],DATA[,'Na2O'],DATA[,'K2O'],DATA[,'P2O5'],DATA[,'CO2.'],DATA[,'LOI']) 

     col.x<-c("sio2","tio2","al2o3","fe2o3","feo","mno","mgo","cao","na2o","k2o","p2o5","co2","loi")
     colnames(x)<-col.x
     #row.names(x) <- Sample
     x [is.na(x)] <- 0 ##Change values NA->0
  
    x$co2 <- ifelse (Cancrinite == TRUE | Calcite == TRUE, x$co2, 0)
    x$loi <- x$loi-x$co2
    Sums.maj <- rowSums(x)
    adj<-function (m) round ((100*m)/(Sums.maj-x$loi), digits)##function of anhydrous adjust
    x[c(1:12)] <- adj (x[,1:12])#
    x$loi <- round (x$loi, digits)#

######################################################################
#          Adjustment of Fe-oxidation (Middlemost, 1989)             #
######################################################################
######################################################################

    TAS = data.frame(S = x$sio2, NK = x$k2o + x$na2o)

#require(SDMTools)#We need this package to define Fe.ratio

#Rhyolite# 
#define the polygon
Rhyolite = cbind(a = c(77,100,100,69,69),
                 b = c(0,0,25,25,8))
#create check which points fall within the polygon
in.polygon.rhy <- pnt.in.poly(TAS, Rhyolite)
rhy <- ifelse (in.polygon.rhy [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (rhy == TRUE, 0.5, 0)
x$Rock <- ifelse (rhy == TRUE, "Rhyolite", "NA")
x$Rock.Intrusive <- ifelse (rhy == TRUE, "Granite", "NA")
#polygon(Rhyolite)


#Phonolite# 
#define the polygon
Phonolite = cbind(a = c(57.6,69,30),
                  b = c(11.7,17.73,24.15))
#create check which points fall within the polygon
in.polygon.pho <- pnt.in.poly(TAS, Phonolite)
pho <- ifelse (in.polygon.pho [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (pho == TRUE, 0.5, x$R.Fe)
x$Rock <- ifelse (pho == TRUE, "Phonolite", x$Rock)
x$Rock.Intrusive <- ifelse (pho == TRUE, "Nephelinic Syenite", x$Rock.Intrusive)
#polygon(Phonolite)


#Trachyte/Trachydacite# 
#define the polygon
Trachyte.Trachydacite = cbind(a = c(63,69,69,57.6),
                              b = c(7,8,17.73,11.7))
#create check which points fall within the polygon
in.polygon.tr.tr <- pnt.in.poly(TAS, Trachyte.Trachydacite)
tr.tr <- ifelse (in.polygon.tr.tr [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (tr.tr == TRUE, 0.5, x$R.Fe)
x$Rock <- ifelse (tr.tr == TRUE, "Trachyte/Trachydacite", x$Rock)
x$Rock.Intrusive <- ifelse (tr.tr == TRUE, "Quartz Monzonite", x$Rock.Intrusive)
#polygon(Trachyte.Trachydacite)


#Dacite# 
#define the polygon
Dacite = cbind(a = c(63,77,69,63),
               b = c(0,0,8,7))
#create check which points fall within the polygon
in.polygon.dac <- pnt.in.poly(TAS, Dacite)
dac <- ifelse (in.polygon.dac [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (dac == TRUE, 0.4, x$R.Fe)
x$Rock <- ifelse (dac == TRUE, "Dacite", x$Rock)
x$Rock.Intrusive <- ifelse (dac == TRUE, "Granodiorite", x$Rock.Intrusive)
#polygon(Dacite)


#Tephriphonolite# 
#define the polygon
Tephriphonolite = cbind(a = c(53,57.6,52.5,48.4),
                        b = c(9.3,11.7,14,11.5))
#create check which points fall within the polygon
in.polygon.teph <- pnt.in.poly(TAS, Tephriphonolite)
teph <- ifelse (in.polygon.teph [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (teph == TRUE, 0.4, x$R.Fe)
x$Rock <- ifelse (teph == TRUE, "Tephriphonolite", x$Rock)
x$Rock.Intrusive <- ifelse (teph == TRUE, "Foide Syenite", x$Rock.Intrusive)
#polygon(Tephriphonolite)


#Trachyandesite# 
#define the polygon
Trachyandesite = cbind(a = c(57,63,57.6,53),
                       b = c(5.9,7,11.7,9.3))
#create check which points fall within the polygon
in.polygon.tr.an <- pnt.in.poly(TAS, Trachyandesite)
tr.an <- ifelse (in.polygon.tr.an [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (tr.an == TRUE, 0.4, x$R.Fe)
x$Rock <- ifelse (tr.an == TRUE, "Trachyandesite", x$Rock)
x$Rock.Intrusive <- ifelse (tr.an == TRUE, "Monzonite", x$Rock.Intrusive)
#polygon(Trachyandesite)


#Andesite# 
#define the polygon
Andesite = cbind(a = c(57,63,63,57),
                 b = c(0,0,7,5.9))
#create check which points fall within the polygon
in.polygon.and <- pnt.in.poly(TAS, Andesite)
and <- ifelse (in.polygon.and [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (and == TRUE, 0.35, x$R.Fe)
x$Rock <- ifelse (and == TRUE, "Andesite", x$Rock)
x$Rock.Intrusive <- ifelse (and == TRUE, "Diorite", x$Rock.Intrusive)
#polygon(Andesite)


#Phonotephrite# 
#define the polygon
Phonotephrite = cbind(a = c(49.4,53,48.4,45),
                      b = c(7.3,9.3,11.5,9.4))
#create check which points fall within the polygon
in.polygon.pho.teph <- pnt.in.poly(TAS, Phonotephrite)
pho.teph <- ifelse (in.polygon.pho.teph [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (pho.teph == TRUE, 0.35, x$R.Fe)
x$Rock <- ifelse (pho.teph == TRUE, "Phonotephrite", x$Rock)
x$Rock.Intrusive <- ifelse (pho.teph == TRUE, "Foide Syenite", x$Rock.Intrusive)
#polygon(Phonotephrite)


#Basaltictrachyandesite# 
#define the polygon
Basaltictrachyandesite = cbind(a = c(52,57,53,49.4),
                               b = c(5,5.9,9.3,7.3))
#create check which points fall within the polygon
in.polygon.ba.tr.an <- pnt.in.poly(TAS, Basaltictrachyandesite)
ba.tr.an <- ifelse (in.polygon.ba.tr.an [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (ba.tr.an == TRUE, 0.35, x$R.Fe)
x$Rock <- ifelse (ba.tr.an == TRUE, "Basaltictrachyandesite", x$Rock)
x$Rock.Intrusive <- ifelse (ba.tr.an == TRUE, "Monzodiorite", x$Rock.Intrusive)
#polygon(Basaltictrachyandesite)


#Basalticandesite# 
#define the polygon
Basalticandesite = cbind(a = c(52,57,57,52),
                         b = c(0,0,5.9,5))
#create check which points fall within the polygon
in.polygon.ba.an <- pnt.in.poly(TAS, Basalticandesite)
ba.an <- ifelse (in.polygon.ba.an [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (ba.an == TRUE, 0.3, x$R.Fe)
x$Rock <- ifelse (ba.an == TRUE, "Basalticandesite", x$Rock)
x$Rock.Intrusive <- ifelse (ba.an == TRUE, "Gabrodiorite", x$Rock.Intrusive)
#polygon(Basalticandesite)


#Tephrite/Basanite#
#define the polygon
Tephrite.Basanite = cbind(a = c(41,41,45,49.4,47),
                          b = c(6,7,9.4,7.3,6))
#create check which points fall within the polygon
in.polygon.teph.bas <- pnt.in.poly(TAS, Tephrite.Basanite)
teph.bas <- ifelse (in.polygon.teph.bas  [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (teph.bas  == TRUE, 0.3, x$R.Fe)
x$Rock <- ifelse (teph.bas  == TRUE, "Tephrite/Basanite", x$Rock)
x$Rock.Intrusive <- ifelse (teph.bas == TRUE, "Foid Monzogabbro/Diorite", x$Rock.Intrusive)
#polygon(Tephrite.Basanite)


#Trachybasalt#
#define the polygon
Trachybasalt = cbind(a = c(45,52,49.4),
                     b = c(5,5,7.3))
#create check which points fall within the polygon
in.polygon.tr.ba<- pnt.in.poly(TAS, Trachybasalt)
tr.ba <- ifelse (in.polygon.tr.ba  [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (tr.ba  == TRUE, 0.3, x$R.Fe)
x$Rock <- ifelse (tr.ba  == TRUE, "Trachybasalt", x$Rock)
x$Rock.Intrusive <- ifelse (tr.ba == TRUE, "Monzogabbro", x$Rock.Intrusive)
#polygon(Trachybasalt)


#Basalt# 
#define the polygon
Basalt = cbind(a = c(45,45,52,52),
               b = c(0,5,5,0))
#create check which points fall within the polygon
in.polygon.bas <- pnt.in.poly(TAS, Basalt)
bas<- ifelse (in.polygon.bas [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (bas == TRUE, 0.2, x$R.Fe)
x$Rock <- ifelse (bas  == TRUE, "Basalt", x$Rock)
x$Rock.Intrusive <- ifelse (bas == TRUE, "Gabbro", x$Rock.Intrusive)
#polygon(Basalt)


#Tephrite/Basanite#
#define the polygon
Basanite.Tephrite = cbind(a = c(41,41,47,45,45),
                          b = c(3,6,6,5,3))
#create check which points fall within the polygon
in.polygon.ba.teph <- pnt.in.poly(TAS, Basanite.Tephrite)
ba.teph<- ifelse (in.polygon.ba.teph [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (ba.teph == TRUE, 0.2, x$R.Fe)
x$Rock <- ifelse (ba.teph == TRUE, "Tephrite/Basanite", x$Rock)
x$Rock.Intrusive <- ifelse (ba.teph == TRUE, "Foid Monzogabbro/Diorite", x$Rock.Intrusive)
#polygon(Basanite.Tephrite)


#Picrobasalt#
#define the polygon
Picrobasalt = cbind(a = c(41,45,45,41),
                    b = c(0,0,3,3))
#create check which points fall within the polygon
in.polygon.pic <- pnt.in.poly(TAS, Picrobasalt)
pic <- ifelse (in.polygon.pic  [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (pic  == TRUE, 0.15, x$R.Fe)
x$Rock <- ifelse (pic  == TRUE, "Picrobasalt", x$Rock)
x$Rock.Intrusive <- ifelse (pic  == TRUE, "Gabbro", x$Rock.Intrusive)
#polygon(Picrobasalt)


#Foidite.1# 
#define the polygon
Foidite.1 = cbind(a = c(46,52.5,30,0,0),
                  b = c(10,14,24.15,24.15,10))
#create check which points fall within the polygon
in.polygon.foi1 <- pnt.in.poly(TAS, Foidite.1)
foi1 <- ifelse (in.polygon.foi1  [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (foi1  == TRUE, 0.4, x$R.Fe)
x$Rock <- ifelse (foi1  == TRUE, "Foidite", x$Rock)
x$Rock.Intrusive <- ifelse (foi1  == TRUE, "Foidite", x$Rock.Intrusive)
#polygon(Foidite.1)


#Foidite.2# 
#define the polygon
Foidite.2 = cbind(a = c(0,46,41,0),
                  b = c(10,10,7,7))
#create check which points fall within the polygon
in.polygon.foi2 <- pnt.in.poly(TAS, Foidite.2)
foi2 <- ifelse (in.polygon.foi2  [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (foi2  == TRUE, 0.3, x$R.Fe)
x$Rock <- ifelse (foi2  == TRUE, "Foidite", x$Rock)
x$Rock.Intrusive <- ifelse (foi2  == TRUE, "Foidite", x$Rock.Intrusive)
#polygon(Foidite.2)


#Foidite.3# 
#define the polygon
Foidite.3 = cbind(a = c(0,41,41,0),
                  b = c(7,7,3,3))
#create check which points fall within the polygon
in.polygon.foi3 <- pnt.in.poly(TAS, Foidite.3)
foi3 <- ifelse (in.polygon.foi3  [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (foi3  == TRUE, 0.2, x$R.Fe)
x$Rock <- ifelse (foi3  == TRUE, "Foidite", x$Rock)
x$Rock.Intrusive <- ifelse (foi3  == TRUE, "Foidite", x$Rock.Intrusive)
#polygon(Foidite.3)


#Foidite.4# 
#define the polygon
Foidite.4 = cbind(a = c(0,41,41,0),
                  b = c(0,0,3,3))
#create check which points fall within the polygon
in.polygon.foi4 <- pnt.in.poly(TAS, Foidite.4)
foi4 <- ifelse (in.polygon.foi4  [3] == 1, TRUE, FALSE)
x$R.Fe <- ifelse (foi4  == TRUE, 0.15, x$R.Fe)
x$Rock <- ifelse (foi4  == TRUE, "Foidite", x$Rock)
x$Rock.Intrusive <- ifelse (foi4  == TRUE, "Ijolite", x$Rock.Intrusive)
#polygon(Foidite.4)


######################################################################
#                     Define type of magma                           #
######################################################################
######################################################################
#Hipersilicic# 
#define the polygon
Hipersilicic = cbind(a = c(66,100,100,66),
                     b = c(0,0,24,24))
#create check which points fall within the polygon
in.polygon.hip <- pnt.in.poly(TAS, Hipersilicic)
hip <- ifelse (in.polygon.hip  [3] == 1, TRUE, FALSE)
x$Type <- ifelse (hip  == TRUE, "Hipersilicic", NA)


#Intermediate# 
#define the polygon
Intermediate = cbind(a = c(52,66,66,52),
                     b = c(0,0,24,24))
#create check which points fall within the polygon
in.polygon.int <- pnt.in.poly(TAS, Intermediate)
int <- ifelse (in.polygon.int  [3] == 1, TRUE, FALSE)
x$Type <- ifelse (int  == TRUE, "Intermediate", x$Type)


#Basic# 
#define the polygon
Basic = cbind(a = c(45,52,52,45),
              b = c(0,0,24,24))
#create check which points fall within the polygon
in.polygon.bas <- pnt.in.poly(TAS, Basic)
bas <- ifelse (in.polygon.bas  [3] == 1, TRUE, FALSE)
x$Type <- ifelse (bas  == TRUE, "Basic", x$Type)


#Ultrabasic# 
#define the polygon
Ultrabasic = cbind(a = c(0,45,45,0),
                   b = c(0,0,24,24))
#create check which points fall within the polygon
in.polygon.Ubas <- pnt.in.poly(TAS, Ultrabasic)
Ubas <- ifelse (in.polygon.Ubas  [3] == 1, TRUE, FALSE)
x$Type <- ifelse (Ubas  == TRUE, "Ultrabasic", x$Type)

######################################################################
#           Adjustment of Fe-oxidation (LeMaitre, 1976)              #
######################################################################
x$Volc <- ifelse (Volcanic == TRUE , TRUE, FALSE)
x$Ox.Fe <- ifelse(x$Volc == TRUE, 0.93-(0.0042 * x$sio2)-(0.022*(x$na2o + x$k2o)), 0.88-(0.0016 * x$sio2)-(0.027*(x$na2o + x$k2o)))

######################################################################
#                      Define FeOx output                            #
######################################################################
x$Magma <- x$Type[,1]
x$Rock.TAS<- ifelse (x$Volc == TRUE, x$Rock,x$Rock.Intrusive)
x$Adj <- ifelse (AdjTAS == TRUE , TRUE, FALSE)
x$Fe.ratio <- ifelse (x$Adj == TRUE , x$R.Fe, x$Ox.Fe)
x$Fe.ratio <- ifelse (x$fe2o3 == 0 | x$feo == 0, x$Fe.ratio, x$fe2o3/x$feo)
x$Fe.ratio <- round(x$Fe.ratio,digits)
######################################################################

######################################################################
#                Define Fe.output (Pruseth, 2009)                    #
######################################################################
T.Fe <- x$fe2o3 + (x$feo*1.111348)
x$feo <- round((T.Fe/(1+(x$Fe.ratio *0.899809)))*0.899809, digits = digits)#
x$fe2o3 <- round(T.Fe-(T.Fe/(1+(x$Fe.ratio *0.899809))), digits = digits)#

######################################################################
#     ##Recalculate major elements data on an anhydrous basis##      #
######################################################################
######################################################################

y <- data.frame(DATA[,'SiO2'],DATA[,'TiO2'],DATA[,'Al2O3'],x$fe2o3,x$feo,DATA[,'MnO'],DATA[,'MgO'],DATA[,'CaO'],DATA[,'Na2O'],DATA[,'K2O'],DATA[,'P2O5'],x$co2,x$loi)
col.y<-c("sio2","tio2","al2o3","fe2o3","feo","mno","mgo","cao","na2o","k2o","p2o5","co2","loi")
colnames(y)<-col.y
y [is.na(y)] <- 0 ##Change values NA->0
Sums.maj <- rowSums(y)
y[,1:12] <- adj (y[,1:12])#
y$Total.adj <- rowSums(y)-y$loi#test out x$Total == 100  

######################################################################
#                         Define Output                              #
######################################################################

output <- data.frame(x$Magma,x$Rock.TAS,x$Fe.ratio,subset(y, select = -13))

##Define output##
col.output<-c("Magma","Rock.TAS","Fe.ratio","SiO2.adj","TiO2.adj","Al2O3.adj","Fe2O3.adj","FeO.adj","MnO.adj","MgO.adj","CaO.adj","Na2O.adj","K2O.adj","P2O5.adj","CO2.adj","Total.adj")
colnames(output) <- col.output
output [output == 0] <- NA  ##Change values 0->NA  
rownames(output) <- DATA$Sample
return(output)

}
## End(Not run)
