summary.hotspots <-
function(object, ...) {
z <- object
if (!inherits(z, "hotspots"))
	stop("use only with \"hotspots\" objects")
h = TRUE ; c = TRUE
if (is.null(z$positive.cut)) h = FALSE
if (is.null(z$negative.cut)) c = FALSE

if (h) {
z$num_phs <- length(z$data[z$data > z$positive.cut])
z$percent_phs <- (length(z$data[z$data > z$positive.cut])/length(z$data))*100
z$percent_phs_sum <- (sum(z$data[z$data > z$positive.cut])/sum(z$data))*100 }

if(c) {
z$num_nhs <- length(z$data[z$data < z$negative.cut])
z$percent_nhs <- (length(z$data[z$data < z$negative.cut])/length(z$data))*100
z$percent_nhs_sum <- (sum(z$data[z$data < z$negative.cut])/sum(z$data))*100 }

m <- mean(z$data)
m[2] <- median(z$data)
m[3] <- min(z$data)
m[4] <- max(z$data)
m[5] <- eval(parse(text = z$var.est))(z$data)
m[6] <- m[5]/m[2]
m <- cbind(m)
dimnames(m) <- list(c("Mean:", "Median:", "Min:", "Max:",
paste(z$var.est, ":", sep = ""), 
paste("CV (", z$var.est, "/median):", sep = "")), "")
z$m <- m

z$disprop <- disprop(z)

class(z) <- "summary.hotspots"
z }