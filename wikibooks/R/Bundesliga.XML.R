`Bundesliga.XML` <-
function(Datei="Bundesliga.xml", Saison="all"){
	data(Bundesliga)
	cat("Daten werden geschrieben...\r")
	sink(Datei)
	if(Saison=="all"){
		saison <- as.character(c("1963/1964", "1964/1965", "1965/1966", "1966/1967", "1967/1968", "1968/1969", "1969/1970", "1970/1971", "1971/1972", "1972/1973", "1973/1974", "1974/1975", "1975/1976", "1976/1977", "1977/1978", "1978/1979", "1979/1980", "1980/1981", "1981/1982", "1982/1983", "1983/1984", "1984/1985", "1985/1986", "1986/1987", "1987/1988", "1988/1989", "1989/1990", "1990/1991", "1991/1992", "1992/1993", "1993/1994", "1994/1995", "1995/1996", "1996/1997", "1997/1998", "1998/1999", "1999/2000", "2000/2001", "2001/2002", "2002/2003", "2003/2004", "2004/2005", "2005/2006", "2006/2007"))
		}else{
		saison <- as.character(Saison)
		}
	saison.l <- length(saison)
	spielnummer<-1
	start.saison <- 1
	ende.saison <- saison.l +1
	cat("<?xml version='1.0' encoding='utf-8'?> \r")
	cat("<Bundesliga>\r")
	
	while(start.saison<ende.saison){   # gehe f<c3><bc>r jede Saison durch
		sjahr <- saison[start.saison]
		cat(sep = "","\t<saison jahr=\"",sjahr , "\">\r")
		saison.liga <- subset(Bundesliga, Bundesliga$Saison==sjahr)
		anzahl.spieltage <- length(levels(as.factor(saison.liga$Spieltag)))
		start.spieltag <- 1
		ende.spieltag <- anzahl.spieltage+1
		while(start.spieltag<ende.spieltag){ # gehe die Spieltage durch
			cat(sep = "","\t\t<spieltag nummer=\"",start.spieltag,"\">\r")
			spieltag.liga <- subset(saison.liga, saison.liga$Spieltag==start.spieltag)
			anzahl.spiele <- length(spieltag.liga$Heim)
			start.spiele <- 1
			ende.spiele <- anzahl.spiele+1
			while(start.spiele<ende.spiele){ # gehe die Begegnungen durch
				cat(sep = "","\t\t\t<spiel nummer=\"",spielnummer ,"\">\r")
				Datum   		<- as.character(spieltag.liga$Datum[start.spiele])
				Anpfiff 		<- as.character(spieltag.liga$Anpfiff[start.spiele])
				Heim			<- as.character(spieltag.liga$Heim[start.spiele])
				Gast			<- as.character(spieltag.liga$Gast[start.spiele])
				Heimtore		<- spieltag.liga$Tore.Heim[start.spiele]
				Gasttore		<- spieltag.liga$Tore.Gast[start.spiele]
				Heimhalbzeit	<- spieltag.liga$Tore.Heim.Halbzeit[start.spiele]
				Gasthalbzeit 	<- spieltag.liga$Tore.Gast.Halbzeit[start.spiele] 
				cat(sep = "","\t\t\t\t<datum>",Datum,"</datum>\r")
				cat(sep = "","\t\t\t\t<anpfiff>",Anpfiff,"</anpfiff>\r")
				cat(sep = "","\t\t\t\t<heim>",Heim,"</heim>\r")
				cat(sep = "","\t\t\t\t<gast>",Gast,"</gast>\r")
				cat(sep = "","\t\t\t\t<heimtore>",Heimtore,"</heimtore>\r")
				cat(sep = "","\t\t\t\t<gasttore>",Gasttore,"</gasttore>\r")
				cat(sep = "","\t\t\t\t<heimhalbzeit>",Heimhalbzeit,"</heimhalbzeit>\r")
				cat(sep = "","\t\t\t\t<gasthalbzeit>",Gasthalbzeit,"</gasthalbzeit>\r")
				cat(sep = "","\t\t\t</spiel>\r")
				start.spiele<- start.spiele+1
				spielnummer<-spielnummer+1
				}
			
			cat("\t\t</spieltag>\r")
			start.spieltag <- start.spieltag+1
			}
		cat("\t</saison>\r")
		start.saison <- start.saison+1	
		}
	
	cat("</Bundesliga>\r")
	
	sink()
	}

