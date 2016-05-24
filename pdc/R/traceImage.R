trace.image <- function(img,resolution=200) 
{
	warning("Call to trace.image(...) is deprecated!");
	return(traceImage(img, resolution));
}

traceImage <- function(img,resolution=200) 
{

w <- dim(img)[1]
h <- dim(img)[2]

center.x <- w/2.0
center.y <- h/2.0

angles <- (1:resolution)/resolution*2*pi
angles.deg <- angles/(2*pi)*360
ts <- c()

# trace angle
for (angle in angles) {

	ts.value <- 0.0
	
	cur.x <- center.x
	cur.y <- center.y
	
	dx <- cos(angle)
	dy <- sin(angle)
	
	repeat
	{
		
		cur.x <- cur.x+dx
		cur.y <- cur.y+dy
		
		round.cur.x <- round(cur.x)
		round.cur.y <- round(cur.y)
		
		if ((round.cur.x >= (w-1) | round.cur.y >= (h-1) 
			| round.cur.x <= 1 | round.cur.y <= 1))  break;
	
		
		# determine black/white value			0=black
		if (img[round.cur.x,round.cur.y]<= 0) {
			a <- as.double(center.x-round.cur.x)**2
			b <- as.double(center.y-round.cur.y)**2
			ts.value <- as.double(sqrt(a+b))
		}	
		
	}
	
	
	ts <- c(ts, ts.value)
}

obj <- c()
obj$distance <- ts
obj$angles <- angles

return(obj);
}
