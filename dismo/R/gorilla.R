
.gorilla <- function(alt, P, P50, TSD) {

	groupSize <- 0
	while (TRUE) {
		groupSize <- groupSize + 1
		
		pLeaves <- min(100, 4.109 + 0.024 * alt)
		nFruits <- max(0, 179.742 - 40.094 * log(pLeaves))
		pFruits <- 100 - pLeaves
		density <- max(0, 0.023 + 0.00035 * alt - 0.21 * P50 + 0.017 * nFruits)
		if (density <= 0) { return(0) }
		
		P[P > 2500] <- 2500  # RH
		Feed = 73.266 - 0.016 * P + 6.16 * density

		Move = 10.053 + 0.245 * groupSize - 0.002 * alt
		Rest = 33.98 - 0.406 * pFruits + 30.168 * TSD
		Groom = 1.55 + 0.23 * groupSize
		
		time <- Feed + Move + Rest + Groom
		
		if (time > 100) {
			return(groupSize - 1)
		}
	}
}



#TmoSD	Pann	Pmo		%Fr		Moist	P2T		Grpsize		Partysize	Grpbm	%Fruit		%Feed		%Move
#0.776	1684.95	143.85	74.3	0.571	10.1275	25.351875	9.4565		1355.15	54.31428571	48.71947368	14.90133333
# moimomx = maximum monthly moisture index
# alt = altitude
# Tann = mean annual temperature
# TmsSD = monthly variation in temperature
# Pmo = average monthly rainfall (mm)
# Pann = mean annual rainfall (mm)
# P2T = plant productivity index (the number of months in the year in which rainfall [in mm] was more than twice the average monthly temperature (Le Houerou 1984))


.apes <- function(alt=1000, Tann=22.5, TmoSD=0.77, Pann=1680, P2T=10, forestcover=80, moimomx=0.5, species='gorilla') {
	if (species == 'chimp') {
		bm <- 40
	} else if (species == 'gorilla') {
		bm <- 120
	} else {
		stop('unkown species')
	}

	Pmo=Pann/12

	groupSize <- 0
	while (TRUE) {
		groupSize <- groupSize + 1

		if (species == 'chimp') {
			partySize <- 21.49 + 0.07 * forestcover - 0.33*Pmo + 0.0012*Pmo^2
		} else {
			partySize <- groupSize
		} 
		pFruit <- 169.43 - 50.65 * log10(bm) - 0.02 * alt - 62.02 * moimomx + 0.39 * forestcover
		pLeaves <- 100 - pFruit
		grpbm <- 4.24 * bm + 29.83 * groupSize

		feeding <- 33.09 + 0.005 * grpbm + 0.14 * bm + 0.16 * pFruit - 0.006*Pann
		moving <- 18.74 + 13.92 * TmoSD + 0.35*partySize - 4.94*(P2T) + 0.32 * P2T^2 
		resting <- -29.47+1.28 * Tann + 0.34*pLeaves + 5.95 * TmoSD 
		grooming <- 1.01+0.23 * groupSize 
	
		time <- feeding + moving + resting + grooming
		if (time > 100) {
			return(groupSize - 1)
		}	
	}
}

