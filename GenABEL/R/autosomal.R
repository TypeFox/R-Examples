"autosomal" <-
function (data) {
	if (is(data,"gwaa.data")) 
		out <- data@gtdata@snpnames[data@gtdata@chromosome != "X" & 
						data@gtdata@chromosome != "XY" &
						data@gtdata@chromosome != "Y" &
						data@gtdata@chromosome != "mt"]
	else if (is(data,"snp.data")) 
		out <- data@snpnames[data@chromosome != "X" & 
						data@chromosome != "XY" &
						data@chromosome != "Y" &
						data@chromosome != "mt"]
	else stop("data argument should have class gwaa.data or snp.data")
	out
}
