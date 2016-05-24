local({
## Vorbereiten
require(klausuR)

## Berechne
klsr.obj <- klausur(data=klsr.data.obj)
## Drucke Ergebnisse
rk.header("klausuR: Test Evaluation")
rk.print("<h3>Global Results</h3>")
rk.print(klsr.obj@results)
rk.print("<h3>Anonymous Feedback</h3>")
rk.print(klsr.obj@anon)
rk.print("<h3>Mark definitions (effectively)</h3>")
rk.print(klsr.obj@marks.sum)
rk.print("<h3>Mean</h3>")
rk.print(klsr.obj@mean)
rk.print("<h3>Standard Deviation</h3>")
rk.print(paste("sd: ",round(klsr.obj@sd, 2), sep=""))
rk.print("<h3>Cronbach&apos;s alpha</h3>")
if(!is.na(klsr.obj@cronbach$alpha)){
  cr.alpha <- paste("alpha: ", round(klsr.obj@cronbach$alpha, 2), "<br />Confidence interval: ",
	      round(klsr.obj@cronbach$ci$LCL, 2),"-",
	      round(klsr.obj@cronbach$ci$UCL, 2)," (95%)",sep="")}
else {cr.alpha <- "Error: Cronbach's alpha is NA"}
rk.print(cr.alpha)
rk.print("<h3>Item Anlysis</h3>")
if(length(klsr.obj@item.analysis) > 1 && !is.na(klsr.obj@item.analysis)){
  item.analysis <- data.frame(Diffc=klsr.obj@item.analysis$Difficulty,
	      DiscrPwr=klsr.obj@item.analysis$Item.total,
	      PartWhole=klsr.obj@item.analysis$Item.Tot.woi,
	      Discrim=klsr.obj@item.analysis$Discrimination,
	      alphaIfDeleted=klsr.obj@item.analysis$alphaIfDeleted)
  dimnames(item.analysis)[[1]] <- dimnames(klsr.obj@item.analysis)[[1]]
} else {
  item.analysis <- "Error: Item analysis is NA"
}
rk.print(item.analysis)
assign("klsr.obj", klsr.obj, envir=globalenv())
})
