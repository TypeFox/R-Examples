ograph <- function( # maaritellaan yleisfunktio piirtamiseen
		ovariable, 
		x, 
		y = character(), 
		type = character(), 
		other = character(),
		fill = NA, 
		...
) {
	if(class(ovariable) == "ovariable")  {
		if(nrow(ovariable@output) == 0) ovariable <- EvalOutput(ovariable)
		data <- ovariable@output
		title <- ovariable@name
		if(length(y) == 0) y <- paste(title, "Result", sep = "")
	} else {
		data <- ovariable
		title <- character()
		if(length(y) == 0) y <- "Result"
	}
	if(length(type) == 0) {
		if("Iter" %in% colnames(data)) type <- geom_boxplot() else type <- geom_bar(stat = "identity")
	}
	out <- ggplot(data, aes_string(x = x, y = y, fill = fill)) # maaritellaan kuvan sarakkeet
	out <- out + type
	out <- out + theme_grey(base_size=24) # Fontin kokoa suurennetaan
	out <- out + labs(
			title	= title,
			y = paste(unique(data[[paste(title, "Yksikk\u00f6", sep = "")]]), sep = "", collapse = ", ")
	)
	out <- out + theme(axis.text.x = element_text(angle = 90, hjust = 1)) # X-akselin tekstit kaannetaan niin etta mahtuvat
	if(length(other) != 0) out <- out + other
	return(out)
}