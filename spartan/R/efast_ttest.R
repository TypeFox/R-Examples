efast_ttest <-
function(Si,rangeSi,STi,rangeSTi,OUTPUTMEASURES_TO_TTEST,NUMPARAMS,NUMCURVES,TTEST_CONF_INT)
{
	# ARRAY FOR P-VALUES FOR THE Si COMPONENT
	p_Si<-array(0,dim=c(NUMPARAMS,1,length(OUTPUTMEASURES_TO_TTEST)))
	# ARRAY FOR P-VALUES FOR THE STi COMPONENT
	p_STi<-array(0,dim=c(NUMPARAMS,1,length(OUTPUTMEASURES_TO_TTEST)))

	# NOW DO ALL THE OUTPUT MEASURES REQUESTED IN THE OUTPUTMEASURES_TO_TTEST RANGE
	# THIS SHOULD BE SOMETHING LIKE 1:3 (i.e. TO EXAMINE THE OUTPUT MEASURES 1 2 AND 3
	for(OUTMEASURE in seq(OUTPUTMEASURES_TO_TTEST))
	{
		# Compare Si or STi of parameter j with the dummy

		for(PARAM in 1:(NUMPARAMS-1))   # MINUS ONE AS ALL COMPARISONS ARE BEING DONE AGAINST THE DUMMY
		{
			# Si done first
			tTestRes<-t.test(rangeSi[PARAM,,OUTMEASURE],rangeSi[NUMPARAMS,,OUTMEASURE],alternative="greater",var.equal=FALSE,conf.level=TTEST_CONF_INT)
		
			# put in array of p values
			p_Si[PARAM,,OUTMEASURE] <- tTestRes$p.value

			# NOW DO STI
			tTestRes<-t.test(rangeSTi[PARAM,,OUTMEASURE],rangeSTi[NUMPARAMS,,OUTMEASURE],alternative="greater",var.equal=FALSE,conf.level=TTEST_CONF_INT)

			# put in array of values
			p_STi[PARAM,,OUTMEASURE] <- tTestRes$p.value

		}
		# SO THAT THE p_Si and p_STi arrays can be output with the rest, these need to be the same size.  At the moment these do not include the dummy,
		# making this difficult - so add a Null to make the lengths match
		p_Si[NUMPARAMS,,OUTMEASURE] <- 0
		p_STi[NUMPARAMS,,OUTMEASURE] <- 0
	}

	return(list(p_Si=p_Si,p_STi=p_STi))
}

