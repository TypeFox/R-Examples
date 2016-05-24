
clean.zipcodes <- structure(function#clean up and standardize zip codes
### Attempts to detect and clean up suspected zip codes.
### Will strip "ZIP+4" suffixes to match format of zipcode data.frame. 
### Restores leading zeros, converts invalid entries to NAs, and returns character vector.
### Note that this function does not attempt to find a matching zip code in the database,
### but rather examines formatting alone.
(zips 
### character vector of suspect entries, will be cast if non-character
){
	
	zips = as.character(zips)
	
	# first strip all whitespace:
	zips = gsub('\\s', '', zips)

	# only keep numbers before a dash, and only if followed by numbers
	zips = gsub('^([0-9]+)-[0-9]+', "\\1", zips)
	
	# delete anything with a non-digit:
	zips[grepl('[^0-9]', zips)] = NA
	
	# delete empty strings
	zips[grepl('^$', zips)] = NA
	
	# restore leading zeroes to 3, 4 digit codes:
	zips = gsub('^([0-9]{3})$', '00\\1', zips)
	zips = gsub('^([0-9]{4})$', '0\\1', zips)
	
	return(zips)
	### character vector containing cleaned zip codes with NAs for non-conforming entries
},ex=function(){
	
# given a mix of possible zip codes, including ZIP+4 and foreign postal codes,
# attempt to identify valid zip codes and return character vector:

zips = c(2061, "02142", 2043, "20210", "2061-2203", "SW1P 3JX", "210", '02199-1880')
	
clean.zipcodes(zips)
# [1] "02061" "02142" "02043" "20210" "02061" NA      "00210" "02199"	
})
