`Bundesliga.Mannschaft` <-
function(Mannschaft, Saison="all"){
	data(Bundesliga)
	
	if(Saison=="all"){
		saison<- as.character(c("1963/1964", "1964/1965", "1965/1966", "1966/1967", "1967/1968", "1968/1969", "1969/1970", "1970/1971", "1971/1972", "1972/1973", "1973/1974", "1974/1975", "1975/1976", "1976/1977", "1977/1978", "1978/1979", "1979/1980", "1980/1981", "1981/1982", "1982/1983", "1983/1984", "1984/1985", "1985/1986", "1986/1987", "1987/1988", "1988/1989", "1989/1990", "1990/1991", "1991/1992", "1992/1993", "1993/1994", "1994/1995", "1995/1996", "1996/1997", "1997/1998", "1998/1999", "1999/2000", "2000/2001", "2001/2002", "2002/2003", "2003/2004", "2004/2005", "2005/2006", "2006/2007"))
		}else{
		saison <- as.character(Saison)
		}
	
	start <- 1
	ende  <- length(saison)+1
	team <- as.character(Mannschaft)
	Uebersicht <- data.frame("saison", "datum", "spieltag","heim", "gast", "ergebnis1", "ergebnis1", "Halbzeit","Halbzeit")
	Uebersicht <- Uebersicht[-1,]
	
	while(start<ende){
		jahr <- as.character(saison[start])
		season <- subset(Bundesliga, Bundesliga$Saison==jahr&(Bundesliga$Heim==Mannschaft|Bundesliga$Gast==Mannschaft))
		dummy <- data.frame(season$Saison, season$Datum, season$Spieltag, season$Heim, season$Gast, season$Tore.Heim, season$Tore.Gast,season$Tore.Heim.Halbzeit, season$Tore.Gast.Halbzeit)
		
		Uebersicht <- rbind(Uebersicht, dummy)
		start<-start+1
		}
	colnames(Uebersicht) <- c("Saison", "Datum", "Spieltag", "Heim", "Gast", "Tore.Heim", "Tore.Gast", "Tore.Heim.Halbzeit","Tore.Gast.Halbzeit")
	return(Uebersicht)	
	}

