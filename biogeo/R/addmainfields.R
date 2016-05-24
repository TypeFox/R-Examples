addmainfields <-
function (dat,species) 
{
  if(is.na(species)){
    species <- 'Species'
    Species <- NA
  } else {
    Species <- dat[,species]
  }
  reqNames <- c("ID",species,"x","y","x_original","y_original","Correction","Modified","Exclude","Reason")
  missingNames <- reqNames[!sapply(reqNames,FUN=function(x) x%in%names(dat))]
  z <- data.frame(dat, ID=1:nrow(dat), Species, x=NA, y=NA, x_original = NA, y_original = NA, Correction = "........", 
                  Modified = "01-01-1900 12:01:01", Exclude = 0, Reason = "........", stringsAsFactors = F)
  z <- z[,c(names(dat),missingNames)]
  return(z)
}
