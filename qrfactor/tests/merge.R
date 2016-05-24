#libabry(qrfactor)
#tablematch="Country"
#gismatch="COUNTRY"
#file<- system.file("external","freshwater.csv", package = "qrfactor")
#file="C:/Users/george/Documents/Rpackages/qrfactor14/inst/external/freshwater.csv"
#table=read.csv(file)
#summary(table)
#
#source<- system.file("external", package = "qrfactor")
#layer="Africawater"
#gisdata2 <- na.omit(readOGR(source, layer))
#summary(gisdata2)

#row.names(table)=table$Country
#row.names(slot(gisdata2, "data"))=gisdata2[[gismatch]]
#row.names(gisdata2)=row.names(slot(gisdata2, "data"))
#o <- match(gisdata2[[gismatch]], table[[tablematch]])
#variables1=table[o,]
#names(variables1)=paste(crossvar[2],names(variables1),sep=".")
#gisdata<- spCbind(gisdata2, variables1)
#names(gisdata)
#summary(gisdata)
#var=c("Annual.Renewable.Water.Resources..km.3.yr.","Total.Freshwater.Withdrawal",
#"Per.Capita.Withdrawal..m.3.p.yr.","Domestic.Use..m.3.p.yr." ,"Industrial.Use..m.3.p.yr.","Agricultural.Use..m.3.p.yr.",
# "Population..millions.","Per.Capita.Water.Resources..km.3.p.yr." )
#mod=qrfactor(gisdata,var=var,scale="log")
#plot(mod,cex=c("Domestic.Use..m.3.p.yr."),rowname=c("COUNTRY"),values=TRUE,pch=23)
#plot(mod,cex=c("Annual.Renewable.Water.Resources..km.3.yr."),rowname=c("COUNTRY"),values=TRUE,pch=23)

###################################################################################################
#summary(mod$gisdata)
