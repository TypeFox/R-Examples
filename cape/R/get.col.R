get.col <-
function(col.name, light.dark = c("light", "dark")){
	
	possible.cols <- c("green", "purple", "orange", "blue", "brown", "yellow", "gray")		
	col.check <- match(col.name, possible.cols)

	if(is.na(col.check)){
		cat("Possible colors are:", possible.cols, sep = "\n")
		stop()
		}

	if(length(light.dark) > 1){
		light.dark <- "dark"
		}
	
	if(light.dark == "light"){
		# color order: same as in arguments, light to dark 
		all.col.ref <- matrix(
		c("#edf8e9", "#c7e9c0", "#a1d99b",
		"#f2f0f7", "#dadaeb", "#bcbddc",
		"#ffffd4", "#fee391", "#fec44f",
		"#eff3ff", "#c6dbef", "#9ecae1",
		"#f6e8c3", "#dfc27d", "#bf812d",
		"#ffffcc", "#ffeda0", "#fed976",
		"#f0f0f0", "#d9d9d9", "#bdbdbd"), nrow = 3, byrow = FALSE)
		}else{
		#color order: same as in arguments, light to dark 
		all.col.ref <- matrix(
		c("#f7fcf5", "#a1d99b", "#238b45",
		"#efedf5", "#9e9ac8", "#6a51a3",
		"#fff7bc", "#fe9929", "#cc4c02",
		"#deebf7", "#6baed6", "#08519c",
		"#f6e8c3", "#bf812d", "#8c510a",
		"#fff9cc", "#fff07f", "#ffe200",
		"#f0f0f0", "#969696", "#252525"), nrow = 3, byrow = FALSE)
		}
	
	
	colnames(all.col.ref) <- c("green", "purple", "orange", "blue", "brown", "yellow", "gray")

	col.locale <- which(colnames(all.col.ref) == col.name)
	
	return(all.col.ref[,col.locale])
	
	}
