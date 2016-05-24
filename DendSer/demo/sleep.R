# Sleep example as in Section 4.3 of Advances in dendrogram seriation for application to visualization, D. Earle and C. Hurley



library(DendSer)


for (pk in c("scagnostics","alr3","RColorBrewer")){
	if (!(pk %in% rownames(installed.packages()))) install.packages(pk)

}

library(scagnostics); library(alr3); 
library(RColorBrewer)



data(sleep1)
dat <- na.omit(sleep1)
colnames(dat) <- c("SW","PS" ,"TS" ,"Bd", "Br","L","GP","P" ,"SE" , "D"  )
#--------------
dat<- na.omit(dat)
sc <- scagnostics(dat)
g <- as.matrix(scagnosticsGrid(sc));

which.sc <- "Outlying"

d <- diag(ncol(dat))
d[g] <- sc[which.sc,]# outliers of plots
d <- 1-as.dist(t(d))



brks <- quantile(sc[which.sc,],c(0,.5,.8,1))
d.color <- dmat.color(1-d,brewer.pal(9,"Greens")[c(1,3,5)],breaks=brks)
#d.color <- dmat.color(1-d)


dm <- as.matrix(d)
o<- dser(d,cost=costLPL)

dev.new(width=5.5, height=5.25)  # Figure 13
cpairs(dat, border.col=1,order=o, oma=c(1,1,1,1), panel.colors=d.color,yaxt="n",xaxt="n", gap=0.2, cex.main=1.5,pch=20,col="grey40")


