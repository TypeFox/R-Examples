report.progress <-
function(current, total, percent.text = 10, percent.dot = 2){

			all.iterations <- 1:total
			percent.prog <- round(all.iterations/total, 2)*100

			rounded.percent.write <- percent.prog%/%percent.text
			rounded.percent.dot <- percent.prog%/%percent.dot
			
			current.locale <- which(all.iterations == current)
			current.percent.write <- rounded.percent.write[current.locale]
			current.percent.dot <- rounded.percent.dot[current.locale]			
			
			curr.percent.write.locale <- which(rounded.percent.write == current.percent.write)
			curr.percent.dot.locale <- which(rounded.percent.dot == current.percent.dot)
			
			
			if(current.locale == min(curr.percent.write.locale)){
				cat(rounded.percent.write[current.locale]*percent.text, "%.", sep = "")
				}else{
					if(current.locale == min(curr.percent.dot.locale)){
						cat(".")
						}					
					}
		
		}
