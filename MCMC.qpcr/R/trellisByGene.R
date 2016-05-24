trellisByGene=function(modelSummary,xFactor,groupFactor,nrow=2,lineWidth=0.4,whiskerWidth=0.2,pointSize=2.5, facetScales="free_y",ylab="log(abundance)",legendPos="bottom",posDodge=0.3){
	lower=NULL;upper=NULL
	pd = position_dodge(posDodge)
	gpl=ggplot( modelSummary$summary, 
	 aes_string(x = xFactor, group = groupFactor, colour = groupFactor, y = "mean")) + 
	 geom_errorbar(aes(ymin = lower, ymax = upper), lwd = lineWidth, width = whiskerWidth, position = pd) + 
	 geom_line(position = pd) + geom_point(position = pd, size = 2.5) + 
	 theme_bw() +facet_wrap(~gene,scales=facetScales,nrow=nrow) +
	 theme(strip.text.x=element_text(face="italic",size = 12, family="serif")) +
	 ylab(ylab) +
	 theme(legend.position=legendPos)
	return(gpl)}