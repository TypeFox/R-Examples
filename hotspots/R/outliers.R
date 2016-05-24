outliers <-
function(x, p = 0.99, tail = "positive", distribution = "t", var.est = "mad", 
	center.est = "mean") {
u <- eval(parse(text = center.est))(x[!is.na(x)])
z <- hotspots(x-u, p = p, tail = tail, 
	distribution = distribution, var.est = var.est)
z$u <- u
z$tail <- tail
z$x <- z$x + u
z$data <- z$data + u
if (tail == "positive" | tail == "both")
	z$positive.cut <- z$positive.cut + u
if (tail == "negative" | tail == "both")
	z$negative.cut <- z$negative.cut + u
z$dataset_name <- deparse(substitute(x))
z }

