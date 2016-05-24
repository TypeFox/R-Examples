compare_models <- function(out1,out2,out3=NULL,out4=NULL,out5=NULL,nested=FALSE){

	if(!nested){
		tab = rbind(model1=c("log-lik."=out1$lk, np=out1$np, BIC=out1$bic),
    	              model2=c(out2$lk, out2$np, out2$bic))
    	if(!is.null(out3)) tab = rbind(tab,model3=c(out3$lk, out3$np, out3$bic))
    	if(!is.null(out4)) tab = rbind(tab,model4=c(out4$lk, out4$np, out4$bic))
    	if(!is.null(out5)) tab = rbind(tab,model5=c(out5$lk, out5$np, out5$bic))
	}else{
		dev2 = 2*(out1$lk - out2$lk); df2 = out1$np - out2$np
		tab = rbind(model1=c("log-lik."=out1$lk, np=out1$np, BIC=out1$bic, "LR(vs. model1)"=NA, df=NA,"p-value"=NA),
    	            model2=c(out2$lk, out2$np, out2$bic, dev2, df2, 1-pchisq(dev2, df2)))
    		if(!is.null(out3)){
			dev3 = 2*(out1$lk - out3$lk); df3 = out1$np - out3$np    		
    			tab = rbind(tab,model3=c(out3$lk, out3$np, out3$bic, dev3, df3, 1-pchisq(dev3, df3)))
    		}
	    	if(!is.null(out4)){
			dev4 = 2*(out1$lk - out4$lk); df4 = out1$np - out4$np    		
    			tab = rbind(tab,model4=c(out4$lk, out4$np, out4$bic, dev4, df4, 1-pchisq(dev4, df4)))
    		}
	    	if(!is.null(out5)){
			dev5 = 2*(out1$lk - out5$lk); df5 = out1$np - out5$np
			tab = rbind(tab,model5=c(out5$lk, out5$np, out5$bic, dev5, df5, 1-pchisq(dev5, df5)))
    		}
	}
	tab

}