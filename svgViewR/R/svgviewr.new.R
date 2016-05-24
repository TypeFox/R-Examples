svgviewr.new <- function(file, window.title="SVG Viewer", animate.duration = 1, 
	animate.reverse = FALSE, animate.repeat = -1, fdir = NULL){

	# CHECK THAT FILE IS OF TYPE HTML

	if(is.null(file)) return(0)

	n <- rep(NA, 0)

	n <- c(n, "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">")
	n <- c(n, paste("<title>", window.title,"</title>", sep=""))
	n <- c(n, "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" >\n")
	
	if(is.null(fdir)) fdir <- paste0(path.package("svgViewR"), "/extdata/")

	n <- c(n, "<svg_doc xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" style=\"visibility:hidden;\" >")
	n <- c(n, "</svg_doc>\n")

	n <- c(n, "<body style=\"margin:0px;background-color:;overflow:hidden;\" >")
	n <- c(n, "\t<a id='keydown' type='checkbox' onkeydown=\"javascript:;\" ></a>")
	n <- c(n, "\t<div style=\"font-family:Bookman Old Style;font-size:12px;margin:0px;background-color:;overflow:hidden;position:fixed;width:170px;top:20px;right:20px;\" id='control_panel' >")
	n <- c(n, "\t\t<div id='control_panel_layers' style=\"width:200px;\" ></div>")
	n <- c(n, "\t</div>")
	n <- c(n, "\t<svg id=\"world\" style=\"background-color:;width:100%;height:100%;position:absolute;top:0;left:0;z-index:-1;\" onload=\"Javascript:onLoadFunctions();\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\"></svg>")
	n <- c(n, "</body>\n")
	
	n <- c(n, "<script type=\"text/javascript\" >")
	n <- c(n, "\tvar svgDocument = document.getElementById(\"world\");")
	n <- c(n, "</script>\n")

	n <- c(n, "<script type=\"text/javascript\" >\n")

	# DURATION OF ANIMATION INCLUDING REVERSE (IF SELECTED) IN SECONDS
	#	DEPENDING ON NUMBER OF SHAPES, SVG CANNOT EXCEED CERTAIN DURATION
	n <- c(n, "\t// Parameters set by user")
	n <- c(n, paste0("\tvar animation_duration = ", animate.duration, ";"))
	n <- c(n, paste0("\tvar animation_reverse = ", ifelse(animate.reverse, 1, 0), ";"))
	n <- c(n, paste0("\tvar animation_repeat = ", animate.repeat, ";"))
	n <- c(n, paste0("\tvar animation_count = ", 0, ";"))
	n <- c(n, "")

	n <- c(n, paste("\t", paste(readLines(paste0(fdir, "math.js")), collapse="\n\t"), "\n</script>", sep=""))
	n <- c(n, paste("<script type=\"text/javascript\" >\n\t", paste(readLines(paste0(fdir, "ui_functions.js")), collapse="\n\t"), "\n</script>", sep=""))
	n <- c(n, paste("<script type=\"text/javascript\" >\n\t", paste(readLines(paste0(fdir, "shape_operations.js")), collapse="\n\t"), "\n</script>", sep=""))
	n <- c(n, paste("<script type=\"text/javascript\" >\n\t", paste(readLines(paste0(fdir, "shapes.js")), collapse="\n\t"), "\n</script>", sep=""))
	n <- c(n, paste("<script type=\"text/javascript\" >\n\t", paste(readLines(paste0(fdir, "control_panel.js")), collapse="\n\t"), "\n</script>", sep=""))

	write(n, file)
	
	# SET CURRENT SVG FPATH AS GLOBAL SVG SAVE AS FPATH
	#write.SVG.current <<- file
}