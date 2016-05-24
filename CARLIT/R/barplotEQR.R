barplotEQR <-
function(EQR.df){cols <- rep(NA, length(unique(EQR.df[,1])))

for(n in 1:length(cols)){
	
	if(as.numeric(paste(unique(EQR.df[n,5])))<=1) {
		if(as.numeric(paste(unique(EQR.df[,5])))[n]>0.75){cols[n]=c("blue")}
		else {if(as.numeric(paste(unique(EQR.df[,5])))[n]>0.60){cols[n]=c("green")}
			else{if(as.numeric(paste(unique(EQR.df[,5])))[n]>0.40){cols[n]=c("yellow")}
				else{if(as.numeric(paste(unique(EQR.df[,5])))[n]>0.25){cols[n]=c("orange")}
					else{cols[n]=c("red")}}
				}
			}
		}
}	


barplot(unique(EQR.df[,5]), ylim=c(0, 1.1), ylab="", xlim=c(0, (length(unique(EQR.df[,1]))+2)), col=cols)

h <- c(0.25, 0.40, 0.60, 0.75, 1)

for (i in 1:length(h)){
	abline(h=h[i], lty=2)
}

par(new=T)

barplot(unique(EQR.df[,5]), names.arg=unique(EQR.df[,1]), ylim=c(0, 1.1), ylab="EQR", xlim=c(0, (length(unique(EQR.df[,1]))+2)), col=cols)

ES <- c("Bad", "Poor", "Moderate", "Good", "High")

for (i in 1:length(h)){
	text(((length(unique(EQR.df[,1])))+(((length(unique(EQR.df[,1]))+2)-(length(unique(EQR.df[,1]))))/1.5)), h[i]-0.1, ES[i], pos=4)
}
}
