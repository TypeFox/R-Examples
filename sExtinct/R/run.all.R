run.all <-
function(sightingdata, alpha, test.year, data.out, plot){

##some warnings

   if(!is.logical(plot)) 
   		stop('plot must be logical', call.=FALSE)
   
   ##the functions
	O<-OLE(sightingdata, alpha)
	Strauss<-Strauss89(sightingdata, alpha, data.out)
	sol1<-Solow1993.eq2(sightingdata, alpha, test.year, data.out)
	sol2<-Solow2005.eq7(sightingdata, alpha, test.year, data.out)
	Burgmanres<-Burgman(sightingdata, alpha, test.year, data.out)
	Robson<-Robson1964(sightingdata, alpha, data.out)

	if(data.out==F){
		k=6
		results<-data.frame(Test=character(length=k),
						Estimate=numeric(length=k))
		results$Test<-c("OLE", "Strauss", "Solow1993.eq2", "Solow2005.eq7", "Robson", "Burgman")
		results$Estimate[1]<-O$Estimate
		results$Estimate[2]<-Strauss$Estimate
		results$Estimate[3]<-sol1$Estimate
		results$Estimate[4]<-sol2$Estimate
		results$Estimate[5]<-Robson[1]
		results$Estimate[6]<-Burgmanres$Estimate
		results$Estimate<-round(results$Estimate)
		for(t in 1:length(results$Estimate)){
			results$Estimate[t]<-ifelse(results$Estimate[t]>test.year, NA, results$Estimate[t])}
		}
	else{
		sol1$code<-"Solow1993.eq2"
		sol2$code<-"Solow2005.eq7"
		Burgmanres$code<-"Burgman"
		Strauss$code<-"Strauss"
		Robson$code<-"Robson"
		results<-rbind(sol1, sol2, Burgmanres, Strauss, Robson)	
		}
		if(plot==T){
		O<-OLE(sightingdata, alpha)
		Strauss<-Strauss89(sightingdata, alpha, data.out=T)
		Robson<-Robson1964(sightingdata, alpha, T)
		sol1.p<-Solow1993.eq2(sightingdata, alpha, test.year, T)
		sol2.p<-Solow2005.eq7(sightingdata, alpha, test.year, T)
		Burgmanres.p<-Burgman(sightingdata, alpha, test.year, T) 
		sol1.p$code<-"Solow1993.eq2"
		sol2.p$code<-"Solow2005.eq7"
		Burgmanres.p$code<-"Burgman"
		Strauss$code<-"Strauss"
		Robson$code<-"Robson"
		G<-rbind(sol1.p, sol2.p, Burgmanres.p, Strauss, Robson)		
		plotted<-xyplot(G$chance~G$yrs, groups=G$code, 
		type="l", lwd=2, 
		auto.key=list(columns = 2, lines=T, points=F),
		xlab="Time", ylab="Chance of persistance", ylim=c(-0.01, 1),
		panel = function(...) {
           	panel.abline(h = alpha, lty = "dotted")
           	panel.abline(v=O$Estimate, lty = "dashed", col="maroon", label="OLE")
           	panel.xyplot(...)
           	panel.text(x=O$Estimate, y=0.9, labels="OLE", col="maroon")
           	panel.text(x=max(G$yrs), y=alpha+.012, labels="alpha", col="black", cex=0.85)
         	}
           )
		}
		if(plot==T){list(plotted, results)}
		else{return(results)}
}
