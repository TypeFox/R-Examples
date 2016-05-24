print.CopyDetect2<- function(x, ...){

	cat("************************************************************************","\n")
	cat("CopyDetect -- An R Package to Compute Statistical Indices for Detecting","\n")
	cat("             ","Answer Copying on Multiple-Choice Tests","\n")
	cat("","\n")
	cat("Version 1.2  2016","\n")
	cat("","\n")
	cat("Cengiz Zopluoglu","\n")
	cat("","\n")
	cat("Assistant Professor","\n")
	cat("University of Miami - Department of Educational and Psychological Studies","\n")
	cat("Research, Measurement, and Evaluation Program","\n")
	cat("","\n")
	cat("c.zopluoglu@miami.edu","\n")
	cat("*************************************************************************","\n")
	cat("","\n")
	cat("Processing Date: ",date(),"\n")
	cat("","\n")
	cat("                     Answer Key:",x$key,"\n")
	cat("","\n")
	cat("   Suspected Copier: Examinee",x$suspected.pair[1],"\n")
	cat("                Response Vector:",as.character(x$data[x$suspected.pair[1],]),"\n")
	cat("         Scored Response Vector:",sprintf(rep("%1.0f",ncol(x$data)),x$scored.data[x$suspected.pair[1],]),"\n")
	cat("","\n")
	cat("   Suspected Source: Examinee",x$suspected.pair[2],"\n")
 	cat("                Response Vector:",as.character(x$data[x$suspected.pair[2],]),"\n")
	cat("         Scored Response Vector:",sprintf(rep("%1.0f",ncol(x$data)),x$scored.data[x$suspected.pair[2],]),"\n")
	cat("","\n")
	cat("   Number-Correct Score for Examinee",x$suspected.pair[1],":",rowSums(x$scored.data[x$suspected.pair[1],]),"\n")
	cat("   Number-Correct Score for Examinee",x$suspected.pair[2],":",rowSums(x$scored.data[x$suspected.pair[2],]),"\n")
	cat("","\n")
	cat("   MLE Ability Estimate for Examinee",x$suspected.pair[1],":",round(x$theta.par[x$suspected.pair[1]],4),"\n")
	cat("   MLE Ability Estimate for Examinee",x$suspected.pair[2],":",round(x$theta.par[x$suspected.pair[2]],4),"\n")
	cat("","\n")
	incorrect.items <- which(x$scored.data[x$suspected.pair[2],]==0)	
	cat("   Number of Identical Incorrect Responses:",sprintf("%1.0f",length(which(x$data[x$suspected.pair[1],incorrect.items]==x$data[x$suspected.pair[2],incorrect.items]))),"\n")
	cat("   Number of Identical Correct Responses:",sprintf("%1.0f",length(which(x$scored.data[x$suspected.pair[1],]== 1 & x$scored.data[x$suspected.pair[2],]==1))),"\n")
	cat("   Number of Identical Responses:",sprintf("%1.0f",length(which(x$data[x$suspected.pair[1],]==x$data[x$suspected.pair[2],]))),"\n")
	cat("","\n")
	cat("","\n")

	if(is.null(x$W.index)!=TRUE) {

		if(is.na(x$W.index$p.value)!=TRUE){

			if(x$W.index$p.value<.05){  w05  <- "Flagged"} else w05  <- "Not Flagged"
			if(x$W.index$p.value<.01){  w01  <- "Flagged"} else w01  <- "Not Flagged"
			if(x$W.index$p.value<.001){ w001 <- "Flagged"} else w001 <- "Not Flagged"
		} else { w05 <- "NA"; w01 <- "NA"; w001 <- "NA"}
	}

	if(is.null(x$GBT.index)!=TRUE) {

		if(is.na(x$GBT.index$p.value)!=TRUE){

			if(x$GBT.index$p.value<.05){  gbt05  <- "Flagged"} else gbt05  <- "Not Flagged"
			if(x$GBT.index$p.value<.01){  gbt01  <- "Flagged"} else gbt01  <- "Not Flagged"
			if(x$GBT.index$p.value<.001){ gbt001 <- "Flagged"} else gbt001 <- "Not Flagged"

		} else { gbt05 <- "NA"; gbt01 <- "NA"; gbt001 <- "NA"}

	}

      if(is.null(x$K.index)!=TRUE) {

		if(is.na(x$K.index$k.index)!=TRUE){

			if(x$K.index$k.index<.05){  k05  <- "Flagged"} else k05  <- "Not Flagged"
			if(x$K.index$k.index<.01){  k01  <- "Flagged"} else k01  <- "Not Flagged"
			if(x$K.index$k.index<.001){ k001 <- "Flagged"} else k001 <- "Not Flagged"
		} else { k05 <- "NA"; k01 <- "NA"; k001 <- "NA"}
	} 

	if(is.null(x$K.index)==TRUE) { 
		k05  <- ""
		k01  <- ""		
		k001 <- ""

	} 

	if(is.null(x$K.variants)!=TRUE) {

		if(is.na(x$K.variants$K1.index)!=TRUE){
			if(x$K.variants$K1.index<.05){  k1.05  <- "Flagged"} else k1.05  <- "Not Flagged"
			if(x$K.variants$K1.index<.01){  k1.01  <- "Flagged"} else k1.01  <- "Not Flagged"
			if(x$K.variants$K1.index<.001){ k1.001 <- "Flagged"} else k1.001 <- "Not Flagged"
		} else { k1.05 <- "NA"; k1.01 <- "NA"; k1.001 <- "NA"}

		if(is.na(x$K.variants$K2.index)!=TRUE){
			if(x$K.variants$K2.index<.05){  k2.05  <- "Flagged"} else k2.05  <- "Not Flagged"
			if(x$K.variants$K2.index<.01){  k2.01  <- "Flagged"} else k2.01  <- "Not Flagged"
			if(x$K.variants$K2.index<.001){ k2.001 <- "Flagged"} else k2.001 <- "Not Flagged"
		} else { k2.05 <- "NA"; k2.01 <- "NA"; k2.001 <- "NA"}

		if(is.na(x$K.variants$S1.index)!=TRUE){
			if(x$K.variants$S1.index<.05){  s1.05  <- "Flagged"} else s1.05  <- "Not Flagged"
			if(x$K.variants$S1.index<.01){  s1.01  <- "Flagged"} else s1.01  <- "Not Flagged"
			if(x$K.variants$S1.index<.001){ s1.001 <- "Flagged"} else s1.001 <- "Not Flagged"
		} else { s1.05 <- "NA"; s1.01 <- "NA"; s1.001 <- "NA"}

		if(is.na(x$K.variants$S2.index)!=TRUE){
			if(x$K.variants$S2.index<.05){  s2.05  <- "Flagged"} else s2.05  <- "Not Flagged"
			if(x$K.variants$S2.index<.01){  s2.01  <- "Flagged"} else s2.01  <- "Not Flagged"
			if(x$K.variants$S2.index<.001){ s2.001 <- "Flagged"} else s2.001 <- "Not Flagged"
		} else { s2.05 <- "NA"; s2.01 <- "NA"; s2.001 <- "NA"}
	}

	if(is.null(x$K.variants)==TRUE) {

		k1.05  <- ""; k1.01  <- ""; k1.001 <- ""
		k2.05  <- ""; k2.01  <- ""; k2.001 <- ""
		s1.05  <- ""; s1.01  <- ""; s1.001 <- ""
		s2.05  <- ""; s2.01  <- ""; s2.001 <- ""
	}

	cat(sprintf("%50s","Alpha Level"),"\n")
	cat("","\n")
	cat(sprintf("%18s","Index"),
	    sprintf("%12s","0.05"),
          sprintf("%15s","0.01"),
          sprintf("%16s","0.001"),"\n")	
	cat("","\n")
	cat(sprintf("%18s","W"),
	    sprintf("%14s",w05),
	    sprintf("%15s",w01),
	    sprintf("%15s",w001),"\n")	

	cat(sprintf("%18s","GBT"),
	    sprintf("%14s",gbt05),
	    sprintf("%15s",gbt01),
	    sprintf("%15s",gbt001),"\n")	

	cat(sprintf("%18s","K"),
	    sprintf("%14s",k05),
	    sprintf("%15s",k01),
	    sprintf("%15s",k001),"\n")	

	cat(sprintf("%18s","K1"),
	    sprintf("%14s",k1.05),
	    sprintf("%15s",k1.01),
	    sprintf("%15s",k1.001),"\n")	

	cat(sprintf("%18s","K2"),
	    sprintf("%14s",k2.05),
	    sprintf("%15s",k2.01),
	    sprintf("%15s",k2.001),"\n")	

	cat(sprintf("%18s","S1"),
	    sprintf("%14s",s1.05),
	    sprintf("%15s",s1.01),
	    sprintf("%15s",s1.001),"\n")	

	cat(sprintf("%18s","S2"),
	    sprintf("%14s",s2.05),
	    sprintf("%15s",s2.01),
	    sprintf("%15s",s2.001),"\n")	

	cat("","\n")
	cat("","\n")
	cat("*******************************","          W Index          ","*********************************","\n")
 if(is.null(x$W.index)!=TRUE) {
      cat("","\n")
	cat("   Estimated Probabilities of Giving Response Options for Suspected Copier Examinee:","\n")
	cat("","\n")
      cat("                                Examinee",x$suspected.pair[1],"\n")
      cat("       ",sort(as.numeric(unique(x$key))),sep="          ","\n")

	for(i in 1:ncol(x$data)) { 
         oooo = c(sprintf("%10s",paste("Item",i)))
             for(jjjj in 1:length(unique(x$key))) {
                oooo = c(oooo,sprintf("%10.3f",x$GBT.index$probabilities1[i,jjjj]))
             }
         cat(oooo,"\n")
      }	
	cat("","\n")
	cat("   Expected Number of Identical Responses = ",round(x$W.index$exp.match,3),"\n")
	cat("","\n")
	cat("   Standard Deviation of The Expected Number of Identical Responses:","=",
		  round(x$W.index$sd.match,3),"\n")
	cat("","\n")
	cat("   W index value","=",round(x$W.index$W.value,4),"\n")
	cat("","\n")
	cat("   Likelihood of Agreemet (p-value)","=",round(x$W.index$p.value,5),"\n")	
	cat("","\n")
} else cat("Nominal Response Model item parameters are not provided")
	cat("","\n")
	cat("**********************","          Generalized Binomial Test          ","************************","\n")
 if(is.null(x$GBT.index)!=TRUE) {
	cat("","\n")
	cat("   Probability of Matching on Each Item:","\n")
	cat("","\n")
	for(i in 1:ncol(x$data)) { cat(sprintf("%15s",paste("Item",i)),
						 sprintf("%12.3f",x$GBT$prob.match[i]),
 				             "\n") 
                                }
	cat("","\n")
	cat("   Exact Probability Distribution for Number of Matches:","\n")
	cat("","\n")
	for(i in 1:(ncol(x$data)+1)) { cat(sprintf("%30s",paste("Probability of ",i-1," Match",sep="")),
						 sprintf("%12.3f",
                                     as.numeric(as.matrix(x$GBT$exact.prob.dist)[i,2])),
 				             "\n") 
                                }
	cat("","\n")
	cat(paste("   Probability of Observing ",
                length(which(x$data[x$suspected.pair[1],]==x$data[x$suspected.pair[2],]))," or More Matches = ",
                round(x$GBT$p.value,5),sep=""),"\n")
	cat("","\n")
} else cat("Nominal Response Model item parameters are not provided")
	cat("","\n")
	cat("***************************","          K Index          ","*******************************","\n")
	cat("","\n")
 if(is.null(x$K.index)!=TRUE) {
	cat("   Number-Incorrect Score for Examinee",x$suspected.pair[1],"(suspected copier):",
      ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],]),"\n")
	cat("","\n")
	cat("   Number-Incorrect Score for Examinee",x$suspected.pair[2],"(suspected source):",
              ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]),"\n")
	cat("","\n")
      cat("   Number of Examinees in the Subgroup of","\n") 
	cat("   Number-Incorrect Score",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],]),":",
         length(x$K.index$subgroups)
         ,"\n")
	cat("","\n")
	if(length(x$K.index$subgroups)!=0){
        cat("   Number of Identical Incorrect Responses","\n") 
        cat("   between Each Examinee in This Subgroup and Source Examinee:","\n")
        cat("","\n")

	  if(length(x$K.index$subgroups)>20) {

		r=length(x$K.index$subgroups)%/%20
      	remainder=length(x$K.index$subgroups)%%20

		for(i in 1:r) { cat("      ",sprintf(rep("%2.0f",20),
				  x$K.index$emp.agg[(20*i-19):(20*i)]),
	     		       "\n") 
      	              }
		cat("      ",sprintf(rep("%2.0f",remainder),
		    x$K.index$emp.agg[(20*r+1):(20*r+remainder)]),
	     		       "\n")
	  }
	
	  if(length(x$K.index$subgroups)<=20) {

		cat("      ",sprintf(rep("%2.0f",length(x$K.index$subgroups)),
		    x$K.index$emp.agg),"\n")
	  }
        cat("","\n")
	  cat("      Mean:",
          round(mean(x$K.index$emp.agg,na.rm=TRUE),4),"\n")
	  cat("","\n")
	  cat("   Binomial Probability of Matching on an Identical Incorrect Response","\n") 
	  cat("   with Source Examinee for This Number-Incorrect Subgroup:","\n") 
	  cat("","\n")
	  cat("      ",round(mean(x$K.index$emp.agg,na.rm=TRUE),4),
          "/",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]),"=",
	    round(round(mean(x$K.index$emp.agg,na.rm=TRUE),4)/
	    (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],])),4),"\n")
	  cat("","\n")
	  cat("   Probability of Observing",length(which(x$data[x$suspected.pair[1],incorrect.items]==x$data[x$suspected.pair[2],incorrect.items])),
            "or More Identical Incorrect","\n") 
	  cat("   Matches Using Binomial Distribution =",round(x$K.index$k.index,4),"\n") 
	  cat("","\n")
	} else cat("Binomial probability cannot be computed. K-index is not possible to compute.","\n") 

} else cat("It is not possible to compute K index for a source examinee with no incorrect response","\n") 
	cat("","\n")

	cat("***************************","    K Variants(K1,K2,S1)          ","*******************************","\n")
	cat("","\n")
 if(is.null(x$K.variants)!=TRUE) {
	cat(		 sprintf("%6s","Number"),"\n")
	cat(		 sprintf("%8s","Incorrect"),
                   sprintf("%10s","Number"),
                   sprintf("%34s","Observed Average Number of"),
			 sprintf("%31s","Predicted Average Number of"),"\n")	
	cat(		 sprintf("%11s","Score Group"),
                   sprintf("%14s","of Examinees"),
                   sprintf("%29s","Identical Incorrect Matches"),
			 sprintf("%30s","Identical Incorrect Matches"),"\n")

	cat(		 sprintf("%11s",""),
                   sprintf("%14s",""),
                   sprintf("%22s","with Source Examinee"),
			 sprintf("%30s","with Source Examinee"),"\n")
	cat("","\n")
	cat(sprintf("%12s",""),
                   sprintf("%39s",""),
                   sprintf("%12s","Linear"),
                   sprintf("%11s","Quadratic"),
                   sprintf("%10s","Loglinear"),"\n")
	cat(sprintf("%12s",""),
                   sprintf("%39s",""),
                   sprintf("%11s","Model"),
                   sprintf("%8s","Model"),
                   sprintf("%10s","Model"),"\n")
	cat("","\n")
	for(i in 1:(ncol(x$data)+1)) { 
			
		cat(sprintf("%7.0f",i-1),
                sprintf("%13.0f",length(x$K.variants$subgroups[[i]])),
		    sprintf("%23.3f",x$K.variants$mean.iden.incorrect[i]), 
                sprintf("%18.3f",x$K.variants$pred1[i]),
                sprintf("%8.3f",x$K.variants$pred2[i]),
                sprintf("%10.3f",x$K.variants$pred3[i]),"\n")
	}
	cat("","\n")
	cat("   Number-Incorrect Score for Examinee",x$suspected.pair[1],":",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],]),"\n")
	cat("","\n")
	cat("   Number-Incorrect Score for Examinee",x$suspected.pair[2],":",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]),"\n")
	cat("","\n")
	cat("   Binomial Probability of Matching on an Identical Incorrect Response","\n") 
	cat("   with Source Examinee for Number-Incorrect Subgroup",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],]),":","\n") 
      cat("","\n")

	if(x$K.variants$pred1[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]<
	   (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],])) &
	   x$K.variants$pred1[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]>0
         ){

	cat("              ","for K1 index:",
                            round(x$K.variants$pred1[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)],3),
				    "/",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]),"=",
				   round((x$K.variants$pred1[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)])/
				   (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],])),3),"\n")
	} else 
		if(x$K.variants$pred1[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]>
	         (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]))){   
                  cat("              ","for K1 index is set to .999, because predicted value was higher than 1","\n")
            } else
		   if(x$K.variants$pred1[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]<=0){
	                   cat("              ","for K1 index is set to .001, because predicted value was smaller than 0","\n")
               }

	cat("","\n")

	if(x$K.variants$pred2[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]<
	   (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],])) &
	   x$K.variants$pred2[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]>0
         ){

	    cat("              ","for K2 index:",
                            round(x$K.variants$pred2[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)],3),
				    "/",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]),"=",
				   round((x$K.variants$pred2[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)])/
				   (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],])),3),"\n")
	  } else 
		if(x$K.variants$pred2[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]>
	         (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]))){   
                  cat("              ","for K2 index is set to .999, because predicted value was higher than 1","\n")
            } else
		   if(x$K.variants$pred2[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]<=0){
	                   cat("              ","for K2 index is set to .001, because predicted value was smaller than 0","\n")
               }

	cat("","\n")

	if(x$K.variants$pred3[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]<
	   (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]))
         ){

	cat("   Predicted Average Number of Identical Incorrect Response","\n") 
	cat("   with Source Examinee for Number-Incorrect Subgroup",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],]),"\n")  
	cat("   (parameter to be used for S1 index in Poisson Distribution):          ",
      round(x$K.variants$pred3[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)],3),
      "\n")
	} else 
		if(x$K.variants$pred3[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]>
	         (ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]))){   
                  cat("   Predicted average number of identical incoccrect response is set to",
                  ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[2],]),"\n")
                  cat("   ,the maximum possible number of incorrect matches, for S1 index.","\n")
            }
	cat("","\n")
	cat("   Probability of Observing",length(which(x$data[x$suspected.pair[1],incorrect.items]==x$data[x$suspected.pair[2],incorrect.items])),
          "or More Identical Incorrect","\n") 
	cat("   Matches Using Binomial Distribution","\n") 
	cat("","\n")
	cat("     ","--- K1 index p value:",round(x$K.variants$K1.index,3),"\n")
	cat("","\n")
	cat("     ","--- K2 index p value:",round(x$K.variants$K2.index,3),"\n")
	cat("","\n")
	cat("   Probability of Observing",length(which(x$data[x$suspected.pair[1],incorrect.items]==x$data[x$suspected.pair[2],incorrect.items])),
          "or More Identical Incorrect","\n") 
	cat("   Matches Using Poisson Distribution","\n") 
	cat("","\n")
	cat("     ","--- S1 index p value:",round(x$K.variants$S1.index,3),"\n")

} else cat("It is not possible to compute K variants for a source examinee with no incorrect response","\n") 
	cat("","\n")

	cat("***************************","    S2 Index         ","*******************************","\n")
	cat("","\n")

 if(is.null(x$K.variants)!=TRUE) {

	cat(		 sprintf("%6s","Number"),"\n")
	cat(		 sprintf("%8s","Incorrect"),
                   sprintf("%10s","Number"),
                   sprintf("%34s","Observed Average Number of"),
			 sprintf("%38s","Observed Average Weighted Number of"),"\n")	
	cat(		 sprintf("%11s","Score Group"),
                   sprintf("%14s","of Examinees"),
                   sprintf("%29s","Identical Incorrect Matches"),
			 sprintf("%27s","Identical Correct Matches"),"\n")

	cat(		 sprintf("%11s",""),
                   sprintf("%14s",""),
                   sprintf("%22s","with Source Examinee"),
			 sprintf("%29s","with Source Examinee"),"\n")
	

	cat("","\n")
	for(i in 1:(ncol(x$data)+1)) { 
			
		cat(sprintf("%7.0f",i-1),
                sprintf("%13.0f",length(x$K.variants$subgroups[[i]])),
		    sprintf("%23.3f",x$K.variants$mean.iden.incorrect[i]), 
                sprintf("%29.3f",x$K.variants$weighted.iden.correct[i]),
                "\n")
	}
	cat("","\n")
	cat("Predicted Total Number of Incorrect and Weighted Correct Matches","\n")
	cat("Using Loglinear Model:","\n")
	cat("","\n")
	cat(		 sprintf("%12s","Number"),"\n")
	cat(		 sprintf("%15s","Incorrect"),
                   sprintf("%19s","Predicted"),"\n")
      cat(		 sprintf("%17s","Score Group"),
                   sprintf("%14s","Values"),"\n")
	cat("","\n")
	for(i in 1:(ncol(x$data)+1)) { 
			
		cat(sprintf("%13.0f",i-1),
                sprintf("%19.3f",x$K.variants$pred4[i]),
                "\n")
	}
	cat("","\n")
	if(x$K.variants$pred4[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]<
	   ncol(x$data)
         ){

	cat("   Predicted Average Number of Identical Response","\n") 
	cat("   with Source Examinee for Number-Incorrect Subgroup",ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],]),
	    ":",round(x$K.variants$pred4[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)],0),
	"(rounded to the nearest integer)","\n")
	} else 
		if(x$K.variants$pred4[(ncol(x$scored.data)-rowSums(x$scored.data[x$suspected.pair[1],])+1)]>
	         ncol(x$data)
               ){   
                  cat("   Predicted average number of identical response is set to",
                  ncol(x$data),"\n")
                 cat("   ,the maximum possible number of matches, for S2 index.","\n")
		 }
	cat("","\n")
	cat("S2 index p value:",round(x$K.variants$S2.index,3),"\n")
} else cat("")

}



