
print.Cmeasures <- function(x, ...) 
{
	cl <- x$call

	if(class(x)!="Cmeasures"){
		stop("The object x must be a class Cmeasures.")
	}else{
		cat("\n")
		if(x$Frailty){
			cat("----- Concordance Measures on",x$Names.data,"dataset ----- \n")
		}else{
			cat("In Cox proportional hazards models, only one overall value\n")
			cat("of concordance probability estimation is proposed.\n")
			cat("\n")
			cat("----- Concordance Measures on",x$Names.data,"dataset ----- \n")
		}
		cat("\n")
		
		if((x$Frailty) & (x$Nboot>0)){
			neuf <- paste(rep(" ",9),collapse="")
			dix <- paste(rep(" ",10),collapse="")
			cat(paste(rep(" ",16),collapse=""),paste(neuf,"Between",neuf,neuf,"  Within",dix,neuf,"  Overall",neuf),"\n")
		}
		if(!(x$Frailty) & (x$Nboot>0)){
			neuf <- paste(rep(" ",9),collapse="")
			cat(paste(rep(" ",16),collapse=""),paste(neuf,"Overall \n"))
		}
		
		if((x$marginal==0) & !(x$Frailty)) {
			if (x$Nboot > 0) {
				CPEcond <- formatC(x$CPEcond,3,format="f")
				Cunocond <- formatC(x$Cunocond,3,format="f")
				cat("CONDITIONAL        C     Se        C.I \n")
				cat("Gonen & Heller's",	CPEcond[1,1],CPEcond[2,1],paste("[",CPEcond[3,1]," ; ",CPEcond[4,1],"]",sep=""),"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),Cunocond[1,1],Cunocond[2,1],paste("[",Cunocond[3,1]," ; ",Cunocond[4,1],"]",sep=""),"\n")
				if (x$cindex == 1) {
					cindexcond <- formatC(x$cindexcond,3,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),	cindexcond[1,1],cindexcond[2,1],paste("[",cindexcond[3,1]," ; ",cindexcond[4,1],"]",sep=""),"\n")
				}
			}else{
				CPEcond <- formatC(x$CPEcond,3,width=7,format="f")
				Cunocond <- formatC(x$Cunocond,3,width=7,format="f")
				cat("                   Overall \n")
				cat("CONDITIONAL \n")
				cat("Gonen & Heller's",CPEcond[1,1],"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),Cunocond[1,1],"\n")
				if (x$cindex == 1) {
					cindexcond <- formatC(x$cindexcond,3,width=7,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),cindexcond[1,1],"\n")
				}
			}
		}
		
		if((x$marginal==0) & (x$Frailty)) {
			if (x$Nboot > 0) {
				CPEcond <- formatC(x$CPEcond,3,format="f")
				Cunocond <- formatC(x$Cunocond,3,format="f")
				cat("CONDITIONAL        C     Se        C.I       |   C     Se        C.I       |   C     Se        C.I        \n")
				cat("Gonen & Heller's",
				CPEcond[1,1],CPEcond[2,1],paste("[",CPEcond[3,1]," ; ",CPEcond[4,1],"] |",sep=""),
				CPEcond[1,2],CPEcond[2,2],paste("[",CPEcond[3,2]," ; ",CPEcond[4,2],"] |",sep=""),
				CPEcond[1,3],CPEcond[2,3],paste("[",CPEcond[3,3]," ; ",CPEcond[4,3],"]",sep=""),"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),
				Cunocond[1,1],Cunocond[2,1],paste("[",Cunocond[3,1]," ; ",Cunocond[4,1],"] |",sep=""),
				Cunocond[1,2],Cunocond[2,2],paste("[",Cunocond[3,2]," ; ",Cunocond[4,2],"] |",sep=""),
				Cunocond[1,3],Cunocond[2,3],paste("[",Cunocond[3,3]," ; ",Cunocond[4,3],"]",sep=""),"\n")
				if (x$cindex == 1) {
					cindexcond <- formatC(x$cindexcond,3,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),
					cindexcond[1,1],cindexcond[2,1],paste("[",cindexcond[3,1]," ; ",cindexcond[4,1],"] |",sep=""),
					cindexcond[1,2],cindexcond[2,2],paste("[",cindexcond[3,2]," ; ",cindexcond[4,2],"] |",sep=""),
					cindexcond[1,3],cindexcond[2,3],paste("[",cindexcond[3,3]," ; ",cindexcond[4,3],"]",sep=""),"\n")
				}
			}else{
				CPEcond <- formatC(x$CPEcond,3,width=7,format="f")
				Cunocond <- formatC(x$Cunocond,3,width=7,format="f")
				cat("                   Between Within  Overall \n")
				cat("CONDITIONAL \n")
				cat("Gonen & Heller's",CPEcond[1,1],CPEcond[1,2],CPEcond[1,3],"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),Cunocond[1,1],Cunocond[1,2],Cunocond[1,3],"\n")
				if (x$cindex == 1) {
					cindexcond <- formatC(x$cindexcond,3,width=7,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),cindexcond[1,1],cindexcond[1,2],cindexcond[1,3],"\n")
				}
			}
		}
		
		if((x$marginal==1) & !(x$Frailty)) {
			if (x$Nboot > 0) {
				CPEmarg <- formatC(x$CPEmarg,3,format="f")
				Cunomarg <- formatC(x$Cunomarg,3,format="f")
				cat("MARGINAL           C     Se        C.I \n")
				cat("Gonen & Heller's",CPEmarg[1,1],CPEmarg[2,1],paste("[",CPEmarg[3,1]," ; ",CPEmarg[4,1],"]",sep=""),"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),Cunomarg[1,1],Cunomarg[2,1],paste("[",Cunomarg[3,1]," ; ",Cunomarg[4,1],"]",sep=""),"\n")
				if (x$cindex == 1) {
					cindexmarg <- formatC(x$cindexmarg,3,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),cindexmarg[1,1],cindexmarg[2,1],paste("[",cindexmarg[3,1]," ; ",cindexmarg[4,1],"]",sep=""),"\n")
				}
			}else{
				CPEmarg <- formatC(x$CPEmarg,3,width=7,format="f")
				Cunomarg <- formatC(x$Cunomarg,3,width=7,format="f")
				cat("                   Overall \n")
				cat("MARGINAL \n") 
				cat("Gonen & Heller's",CPEmarg[1,1],"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),Cunomarg[1,1],"\n")
				if (x$cindex == 1) {
					cindexmarg <- formatC(x$cindexmarg,3,width=7,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),cindexmarg[1,1],"\n")
				}
			}
		}
		
		if((x$marginal==1) & (x$Frailty)) {
			if (x$Nboot > 0) {
				CPEmarg <- formatC(x$CPEmarg,3,format="f")
				Cunomarg <- formatC(x$Cunomarg,3,format="f")
				cat("MARGINAL           C     Se        C.I       |   C     Se        C.I       |   C     Se        C.I        \n")
				cat("Gonen & Heller's",
				    CPEmarg[1,1],CPEmarg[2,1],paste("[",CPEmarg[3,1]," ; ",CPEmarg[4,1],"] |",sep=""),
				    CPEmarg[1,2],CPEmarg[2,2],paste("[",CPEmarg[3,2]," ; ",CPEmarg[4,2],"] |",sep=""),
				    CPEmarg[1,3],CPEmarg[2,3],paste("[",CPEmarg[3,3]," ; ",CPEmarg[4,3],"]",sep=""),"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),
				    Cunomarg[1,1],Cunomarg[2,1],paste("[",Cunomarg[3,1]," ; ",Cunomarg[4,1],"] |",sep=""),
				    Cunomarg[1,2],Cunomarg[2,2],paste("[",Cunomarg[3,2]," ; ",Cunomarg[4,2],"] |",sep=""),
				    Cunomarg[1,3],Cunomarg[2,3],paste("[",Cunomarg[3,3]," ; ",Cunomarg[4,3],"]",sep=""),"\n")
				if (x$cindex == 1) {
					cindexmarg <- formatC(x$cindexmarg,3,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),
					    cindexmarg[1,1],cindexmarg[2,1],paste("[",cindexmarg[3,1]," ; ",cindexmarg[4,1],"] |",sep=""),
					    cindexmarg[1,2],cindexmarg[2,2],paste("[",cindexmarg[3,2]," ; ",cindexmarg[4,2],"] |",sep=""),
					    cindexmarg[1,3],cindexmarg[2,3],paste("[",cindexmarg[3,3]," ; ",cindexmarg[4,3],"]",sep=""),"\n")
				}
			}else{
				CPEmarg <- formatC(x$CPEmarg,3,width=7,format="f")
				Cunomarg <- formatC(x$Cunomarg,3,width=7,format="f")
				cat("                   Between Within  Overall \n")
				cat("MARGINAL \n") 
				cat("Gonen & Heller's",CPEmarg[1,1],CPEmarg[1,2],CPEmarg[1,3],"\n")
				cat("Uno's",paste(rep(" ",10),collapse=""),Cunomarg[1,1],Cunomarg[1,2],Cunomarg[1,3],"\n")
				if (x$cindex == 1) {
					cindexmarg <- formatC(x$cindexmarg,3,width=7,format="f")
					cat("Harrell's",paste(rep(" ",6),collapse=""),cindexmarg[1,1],cindexmarg[1,2],cindexmarg[1,3],"\n")
				}
			}
		}
		
		cat("\n")
		cat("----- Information -----\n")
		cat("\n")
		print(x$frequencies,row.names=F,digits=3)
		cat("\n")
		print(x$Npairs,row.names=T)
		cat("\n")
	}
}

