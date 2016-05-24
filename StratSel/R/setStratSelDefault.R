setStratSelDefault <-
function(){
    	temp.rob = c("($N:f)","($AIC:f)","($AIC.c:f)")
    	names(temp.rob) = c("N", "AIC", "AIC.c")
    	setSummaryTemplate(StratSel=temp.rob)
 
    	temp.lm = temp.rob[1:3]
    	setSummaryTemplate(lm =temp.lm)
	}
