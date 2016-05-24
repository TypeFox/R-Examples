hotspots <-
function(x, p = 0.99, tail = "positive", distribution = "t", var.est = "mad") {
z <- list(x = x)
z$data <- x[!is.na(x)]
z$distribution <- distribution
z$var.est <- var.est
z$p <- p
z$tail <- tail
z$dataset_name <- deparse(substitute(x))
z$rrms <- (sqrt(eval(parse(text = var.est))(z$data)^2 + median(z$data)^2))
if (distribution == "normal") {
	pos <- median(z$data/z$rrms) + qnorm(p)
	neg <- median(z$data/z$rrms) - qnorm(p) } else
if (distribution == "t") {
	pos <- median(z$data/z$rrms) + qt(p, df = length(z$data)-1)
	neg <- median(z$data/z$rrms) - qt(p, df = length(z$data)-1) } else
stop("Distribution not supported")
if (tail != "positive" & tail != "negative" & tail != "both")
	stop("Tail must be \"positive\", \"negative\", or \"both\"")
if (tail == "positive" | tail == "both")
	z$positive.cut <- pos*z$rrms
if (tail == "negative" | tail == "both")
	z$negative.cut <- neg*z$rrms
class(z) <- "hotspots"
z }

