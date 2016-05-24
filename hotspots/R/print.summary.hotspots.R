print.summary.hotspots <-
function(x, 
digits = max(3, getOption("digits") - 3), p_round = 1, top = 0, ...) {
h = TRUE ; c = TRUE
if (is.null(x$positive.cut)) h = FALSE
if (is.null(x$negative.cut)) c = FALSE
ty <- "hot spot"
if (!is.null(x$u)) ty <- "outlier"

cat("\nSource data: ", x$dataset_name, "\n", "Distribution and probability: ",
	x$distribution, ", ", x$p, "\n", sep = "")
if (h & !c) cat(paste("Tail: positive ", ty,  "s only", sep = ""))
if (c & !h) cat(paste("Tail: negative ", ty,  "s only", sep = ""))
if (h & c) cat(paste("Tail: positive and negative ", ty,  "s", sep = ""))

print.default(x$m, print.gap = 2, quote = FALSE, digits = digits)
cat("\nn = ", length(x$data))

if (h & top > 0) {
cat(paste("\n\n", top, " most disproportionate positive ", ty,  "s", sep = ""))
value <- sort(x$data, dec = TRUE)[1:top]
hn <- sort(x$disprop$positive, dec = TRUE)[1:top]
cat("\n\n")
tp <- cbind(value, hn)
dimnames(tp) <- list(1:top, c("value", "disproportionality"))
print.default(tp, print.gap = 2, quote = FALSE, digits = digits) }

if (c & top > 0) {
cat(paste("\n\n", top, " most disproportionate negative ", ty,  "s", sep = ""))
value <- sort(x$data, dec = TRUE)[1:top]
hn <- sort(x$disprop$negative, dec = TRUE)[1:top]
cat("\n\n")
tp <- cbind(value, hn)
dimnames(tp) <- list(1:top, c("value", "disproportionality"))
print.default(tp, print.gap = 2, quote = FALSE, digits = digits) }

if (h) { 
j <- cbind(cutoff = x$positive.cut, x$num_phs, x$percent_phs, x$percent_phs_sum)
dimnames(j) <- list("", c("Cutoff", paste("number positive ", ty,  "s", sep = ""), 
paste("% positive ", ty,  "s ", sep = ""), "% sum"))
cat(paste("\n\npositive ", ty,  "s:\n",sep = ""))
j[3] <- round(j[3], p_round) ; j[4] <- round(j[4], p_round)
print.default(j, print.gap = 2, quote = FALSE, digits = digits) }

if (c) {
k <- cbind(cutoff = x$negative.cut, x$num_nhs, x$percent_nhs, x$percent_nhs_sum)
dimnames(k) <- list("", c("Cutoff", paste("number negative ", ty,  "s", sep = ""),
paste("% negative", ty,  "s ", sep = ""), "% sum"))
cat(paste("\n\nnegative ", ty,  "s:\n",sep = ""))
k[3] <- round(k[3], p_round) ; k[4] <- round(k[4], p_round)
print.default(k, print.gap = 2, quote = FALSE, digits = digits) }

if (min(x$data)*max(x$data) < 0) {
cat("\n\n")
warning(call. = FALSE, "% sum may be irrelevant because both positive and negative values are present") }}
