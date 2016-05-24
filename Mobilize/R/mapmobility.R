mapmobility <- function(date, username,...){
	mobility <- oh.mobility.read(date=date, username=username,...);
	maps::map("county", "california");
	lines(mobility$lo, mobility$la, type="l", col="red", lwd=10);
}
