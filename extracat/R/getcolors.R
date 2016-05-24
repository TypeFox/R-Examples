

getcolors <- function(N, palette, col.opt=list(), revert = FALSE){
#require(colorspace)
	
	colv <- NULL
	
	if(!("sample" %in% labels(col.opt))){
						sample <- FALSE
					}else{
						sample <- col.opt$sample
					}
	wijffelaars <- c("#5C6BF7", "#FF5959", "#5CCB5C", "#FFB111", "#AA61BB", "#FFFF5F", "#FF89EB", "#91653E", "#C1C1C1", "#5CE5D6", "#C9FF87",
 "#FFE0CC", "#AD2D5C", "#E3C4EF", "#E2D495", "#CCF1FF", "#578E52")
			
	
	if( palette[1] %in% c("hsv","rgb") ){
					col.def <- formals(rainbow)
					if(!("s" %in% labels(col.opt))){
						col.opt$s <- eval(col.def$s)
					}
					if(!("v" %in% labels(col.opt))){
						col.opt$v <- eval(col.def$v)
					}
					if(!("start" %in% labels(col.opt))){
						col.opt$start <- eval(col.def$start)
					}
					if(!("end" %in% labels(col.opt))){
						col.opt$end <- max(N-1,1)/N
					}
					##if(!("alpha" %in% labels(col.opt))){
					##	col.opt$alpha <- eval(col.def$alpha)
					##}
					colv <- rainbow(N,s = col.opt$s, v = col.opt$v, start = col.opt$start, end = col.opt$end, alpha = col.opt$alpha)
					}
				if( palette[1] == "hcl" ){
					col.def <- formals(rainbow_hcl)
					if(!("c" %in% labels(col.opt))){
						col.opt$c <- eval(col.def$c)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("start" %in% labels(col.opt))){
						col.opt$start <- eval(col.def$start)
					}
					if(!("end" %in% labels(col.opt))){
						col.opt$end <- 360 * (N - 1)/N
					}
				colv <- rainbow_hcl(N,c = col.opt$c, l = col.opt$l, start = col.opt$start, end = col.opt$end)
				}
				if( palette[1] %in% c("s","seq","sqt","sqn","sequential") ){
					col.def <- formals(sequential_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h <- eval(col.def$h)
					}
					if(!("c" %in% labels(col.opt))){
						col.opt$c <- eval(col.def$c.)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power <- eval(col.def$power)
					}
					colv <- rev(sequential_hcl(N,h = col.opt$h, c. = col.opt$c, l = col.opt$l, power = col.opt$power))
				}
				if( palette[1] %in% c("d","div","diverging","diverge") ){
					col.def <- formals(diverge_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h <- eval(col.def$h)
					}
					if(!("c" %in% labels(col.opt))){
						col.opt$c <- eval(col.def$c)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power <- eval(col.def$power)
					}
					colv <- diverge_hcl(N,h = col.opt$h, c = col.opt$c, l = col.opt$l, power = col.opt$power)
					
				}
				
				if( palette[1] %in% c("h","heat","heatcolors") ){
					col.def <- formals(heat_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h <- eval(col.def$h)
					}
					if(!("c." %in% labels(col.opt))){
						col.opt$c. <- eval(col.def$c.)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power <- eval(col.def$power)
					}
					colv <- heat_hcl(N,h = col.opt$h, c. = col.opt$c., l = col.opt$l, power = col.opt$power)
					
				}
				if( palette[1] %in% c("t","ter","terrain") ){
					col.def <- formals(terrain_hcl)
					if(!("h" %in% labels(col.opt))){
						col.opt$h <- eval(col.def$h)
					}
					if(!("c." %in% labels(col.opt))){
						col.opt$c. <- eval(col.def$c.)
					}
					if(!("l" %in% labels(col.opt))){
						col.opt$l <- eval(col.def$l)
					}
					if(!("power" %in% labels(col.opt))){
						col.opt$power <- eval(col.def$power)
					}
					colv <- terrain_hcl(N,h = col.opt$h, c. = col.opt$c., l = col.opt$l, power = col.opt$power)
					
				}


				if( palette[1] %in% c("Wijffelaars","w","wijf", "q17") ){
						colv <- wijffelaars[1:N]
				}
				if(is.null(colv)){
					colv <- rep(palette,N)[1:N]
				}
				if("alpha" %in% labels(col.opt)){
						colv <- alpha(colv, col.opt$alpha)
				}
				
if(sample) colv <- sample(colv)
	if(revert){
		return(rev(colv))
	}
	return(colv)
	
}
