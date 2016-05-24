collate.bathy <- function(east,west){
	as.numeric(rownames(west))+360 -> rownames(west)  # new coordinate system
	rbind(east, west) -> collated                     # collate the two matrices into one
	collated[unique(rownames(collated)), ]-> collated # remove the extra antimeridian line
	class(collated)<-"bathy"                          # assign the class bathy to new matrix c
	return(collated)
}