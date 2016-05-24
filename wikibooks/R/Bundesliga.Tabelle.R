`Bundesliga.Tabelle` <-
function(Saison, Spieltag=1, output="Tabelle"){
	data(Bundesliga)
	saison <- as.character(Saison)

	### Anzahl der Mannschaften aenderte sich ab 1965
	if(saison=="1963/1964" | saison=="1964/1965"){
		dummy <- integer(16)
		pdummy <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
		c <- 9 # Aufschlag
		d <- 8 # Anzahl Spiele pro Spieltag
		if(Spieltag > 30) Spieltag <- 30
		}else{
		# 1991/1992 waren es 20 Mannschaften, wegen der DDR-Integration
		if(saison=="1991/1992"){
			if(Spieltag > 38) Spieltag <- 38
			dummy <- integer(20)	
			pdummy <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
			c <- 11 # Aufschlag
			d <- 10  # Anzahl Spiele pro Spieltag

			}else{
			# Ansonsten sind es immer 18 Mannschaften
			if(Spieltag > 34) Spieltag <- 34
			dummy <- integer(18)	
			pdummy <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
			c <- 10 # Aufschlag
			d <- 9  # Anzahl Spiele pro Spieltag
			}
		}
	#---------------------------------------------------
	
	### Einfuehrung der 3-Punkte-Regel ab Saison 1995/96
	punkteregel <- c("1963/1964", "1964/1965", "1965/1966", "1966/1967", "1967/1968", "1968/1969", "1969/1970", "1970/1971", "1971/1972", "1972/1973", "1973/1974", "1974/1975", "1975/1976", "1976/1977", "1977/1978", "1978/1979", "1979/1980", "1980/1981", "1981/1982", "1982/1983", "1983/1984", "1984/1985", "1985/1986", "1986/1987", "1987/1988", "1988/1989", "1989/1990", "1990/1991", "1991/1992", "1992/1993", "1993/1994", "1994/1995")
	regel.test <- is.element(saison, punkteregel)
	if(regel.test==TRUE){
		switch <- 0	  # 2-Punkte-Regel
		}else{
		switch <- 1   # 3-Punkte-Regel
		}
	#--------------------------------------------------------
	sz <- 1
	liga  <- subset(Bundesliga,Bundesliga$Saison==saison)
	liga2 <- subset(liga, liga$Spieltag==1)
	teams <- as.character(liga2$Heim)
	teams <- c(teams, as.character(liga2$Gast))
	#cat(length(teams))
	Tabelle <- data.frame(dummy,teams,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy)
	colnames(Tabelle) <- c("Platz", "Mannschaft", "Spiele", "G","U","V","t","g","Diff", "Punkte")
	if(Saison=="1993/1994"){
		Tabelle$Punkte[Tabelle$Mannschaft=="Dynamo Dresden"] <- -4 # Abzug wegen Lizenzverstoessen
		}
	if(Saison=="1999/2000"){
		Tabelle$Punkte[Tabelle$Mannschaft=="Eintracht Frankfurt"] <- -2 # Abzug wegen Lizenzverstoessen
		}
	if(Saison=="2003/2004"){
		Tabelle$Punkte[Tabelle$Mannschaft=="1. FC Kaiserslautern"] <- -3 # Abzug wegen Lizenzverstoessen
		}
	Platzierung <- data.frame(teams,pdummy)
	#cat(Platzierung)	
	while(sz < (Spieltag+1)){
		
		
		
		### Spieltag-Ergebnisse
		a <- 1
		while(a<c){
			b <- a+d
			sliga <- subset(liga, liga$Spieltag==sz)
			teams <- as.character(sliga$Heim)
			teams <- c(teams, as.character(sliga$Gast))
			heim <- teams[a]
			gast <- teams[b]
			heim.tore <- sliga$Tore.Heim[sliga$Heim==heim]
			gast.tore <- sliga$Tore.Gast[sliga$Gast==gast]
		
			Tabelle$Spiele[Tabelle$Mannschaft==heim] <- (Tabelle$Spiele[Tabelle$Mannschaft==heim] +1)
			Tabelle$Spiele[Tabelle$Mannschaft==gast] <- (Tabelle$Spiele[Tabelle$Mannschaft==gast] +1)
			Tabelle$t[Tabelle$Mannschaft==heim] <- (Tabelle$t[Tabelle$Mannschaft==heim] + heim.tore)
			Tabelle$g[Tabelle$Mannschaft==heim] <- (Tabelle$g[Tabelle$Mannschaft==heim] + gast.tore)
			Tabelle$t[Tabelle$Mannschaft==gast] <- (Tabelle$t[Tabelle$Mannschaft==gast] + gast.tore)
			Tabelle$g[Tabelle$Mannschaft==gast] <- (Tabelle$g[Tabelle$Mannschaft==gast] + heim.tore)
				
			heim.diff <- Tabelle$t[Tabelle$Mannschaft==heim] - Tabelle$g[Tabelle$Mannschaft==heim]
			gast.diff <- Tabelle$t[Tabelle$Mannschaft==gast] - Tabelle$g[Tabelle$Mannschaft==gast]
			Tabelle$Diff[Tabelle$Mannschaft==heim] <- heim.diff
			Tabelle$Diff[Tabelle$Mannschaft==gast] <- gast.diff
			Tabelle$Diff<-as.numeric(Tabelle$Diff)

			if(gast.tore==heim.tore){
			#unendschieden
				Tabelle$U[Tabelle$Mannschaft==heim] <- (Tabelle$U[Tabelle$Mannschaft==heim] +1)
				Tabelle$U[Tabelle$Mannschaft==gast] <- (Tabelle$U[Tabelle$Mannschaft==gast] +1)
				Tabelle$Punkte[Tabelle$Mannschaft==heim] <- (Tabelle$Punkte[Tabelle$Mannschaft==heim] +1)
				Tabelle$Punkte[Tabelle$Mannschaft==gast] <- (Tabelle$Punkte[Tabelle$Mannschaft==gast] +1)
				}
			if(gast.tore<heim.tore){
				#Heim gewinnt
				Tabelle$G[Tabelle$Mannschaft==heim] <- (Tabelle$G[Tabelle$Mannschaft==heim] +1)
				Tabelle$V[Tabelle$Mannschaft==gast] <- (Tabelle$V[Tabelle$Mannschaft==gast] +1)
				if(switch==0){
					Tabelle$Punkte[Tabelle$Mannschaft==heim] <- (Tabelle$Punkte[Tabelle$Mannschaft==heim] +2)
					}else{
					Tabelle$Punkte[Tabelle$Mannschaft==heim] <- (Tabelle$Punkte[Tabelle$Mannschaft==heim] +3)
					}
				}
			if(gast.tore>heim.tore){
				#Gast gewinnt
				Tabelle$V[Tabelle$Mannschaft==heim] <- (Tabelle$V[Tabelle$Mannschaft==heim] +1)
				Tabelle$G[Tabelle$Mannschaft==gast] <- (Tabelle$G[Tabelle$Mannschaft==gast] +1)
				if(switch==0){
					Tabelle$Punkte[Tabelle$Mannschaft==gast] <- (Tabelle$Punkte[Tabelle$Mannschaft==gast] +2)
					}else{
					Tabelle$Punkte[Tabelle$Mannschaft==gast] <- (Tabelle$Punkte[Tabelle$Mannschaft==gast] +3)
					}
				}
			
			a <- a+1
			}
		if(Saison=="1971/1972"){
			Tabelle$Punkte[Tabelle$Mannschaft=="Arminia Bielefeld"] <- 0 # Bundesliga-Skandal
			}
		###		die Sortierung laeuft jetzt ueber die Indexierung, warum das geht und das andere nicht 
		###		weiss ich auch nicht, aber es haengt mit den Vorzeichen zusammen. 
		###		Hatte zuerst mit *-1 rumprobiert, das war aber auch nicht gut... 
		###		check aber noch mal... 
		temp <- rev(order(Tabelle[,"Punkte"],Tabelle[,"Diff"])) 
		Tabelle <- Tabelle[temp,]
		sort <- 1
		rank <- 1
		while(sort< 2*d+1){
			sort.team <- Tabelle$Mannschaft[sort]
			Tabelle$Platz[Tabelle$Mannschaft==sort.team] <- rank
			pt <- sz+1
			Platzierung[Platzierung$teams==sort.team, pt] <- rank
			sort <- sort+1
			rank <- rank+1
			}
		sz <- sz+1
	}
	if(output=="Tabelle")     return(Tabelle)
	if(output=="Platzierung") {
		Reihen <- 1:Spieltag
		Spalten <- Platzierung[,1]
		Platzierung <- data.frame(matrix(unlist(unclass(Platzierung)),nrow=length(Platzierung),byrow=TRUE,dimnames=list(names(Platzierung),Platzierung[,0])))
		Platzierung <- Platzierung[-1,]
		colnames(Platzierung) <- Spalten
		rownames(Platzierung) <- Reihen
		return(Platzierung)
	}
}

