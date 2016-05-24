#' Download an image
#' @param campaign_urn campaign id
#' @param owner owner of the image
#' @param id photo uuid
#' @param ... other parameters
#' @return path to image
#' @export
oh.image.read <- function(campaign_urn, owner, id, ...){
	tf <- oh.call("/image/read ", responseformat="file", campaign_urn=campaign_urn, owner=owner, id=id, ...);
	if(attr(tf,"Content-Type") == "image/jpeg"){
		newname <- paste("/tmp/", id,".jpg", sep="");
	} else if(attr(tf,"Content-Type") == "image/png") {
		#BUG in ohmage //newname <- paste("/tmp/", id,".png", sep="");
		newname <- paste("/tmp/", id,".jpg", sep="");
	} else {
		stop("Unknown image format: ",attr(tf,"Content-Type"),". Please update in oh.image.read.");
	}
		
	file.rename(tf, newname);
	return(newname);
}
