#######################################
### EcoVirtual -Internal Functions ####
#######################################


##############################################
### Island Biogeography and Neutral Theory ###
##############################################

## fuction rich used in 'simHub1' 'simHub2' simHub3'
rich <- function(x)length(unique(x))

## function animaIls used in 'archip' function
animaIsl=function(riq.tempo, ar.isl, locxy, sprain, col_riq=col_riq, S=S)
{
          Nspp=S
          nsppobs=max(riq.tempo)
          maxt=dim(riq.tempo)[1]
          nIsl<-length(ar.isl)
          comp.max<-max(ar.isl)
          tempo=length(riq.tempo)
          col_spp=rainbow(max(riq.tempo))
          col_func1=colorRamp(c("white", "green4"))
          col_func2=colorRamp(c("green4", "royalblue"))
          col_riq1=rgb(col_func1(seq(0,1, length.out=nsppobs)), maxColorValue=255) 
          col_riq2=rgb(col_func2(seq(0,1, length.out=Nspp-nsppobs)), maxColorValue=255) 
          col_riq=c(col_riq1,col_riq2)
          ## aqui inicia o grafico
          layout(matrix(data=c(2,1), nrow=2, ncol=1), widths=c(1,1), heights=c(5,1))
          old<-par(mar=c(2,2,1,3))
          image(x=1:Nspp, y=1, matrix(data=1:Nspp, nrow=Nspp,ncol=1),col=col_riq, ylab="",xlab=paste("cycle", 1:length(maxt)), xaxt="n", yaxt="n", main="Richness")
          axis(3, at=c(1.5,Nspp),tick=FALSE, labels=c("0", Nspp), mgp=c(0,0,0))
          polygon(x=c(1.5,1.5,Nspp,Nspp), y=c(0.6,1.4,1.4,0.6), lwd=2)
          plot(0:comp.max, 0:comp.max, usr=c(0,comp.max,0,comp.max), type="n", yaxt="n", xaxt="n", xlab="", ylab="", bty="n", main="Passive Sampling and Area ",mar=c(0,2,3,2), oma=c(0,0,0,0))
          segments(x0=c(0,0,comp.max,0), y0=c(0,0,0,comp.max), x1=c(0,rep(comp.max,3)), y1=c(comp.max,0,comp.max,comp.max))
          segments(x0=c(rep(0,nIsl), ar.isl), y0=c(ar.isl,rep(0,nIsl)), x1=c(ar.isl,ar.isl), y1=c(ar.isl,ar.isl))
          for (i in 2:maxt)
          {
                    lxy=locxy[[i]]
                    nspp=riq.tempo[i,]
                    for(f in nIsl:1)
                    {
                              vert=ar.isl[f]
                              polygon(x=c(0,vert, vert,0),y=c(0,0,vert,vert), col=col_riq[nspp[f]] )
                    }
                    points(lxy[,1],lxy[,2], col=col_spp[sprain[[i]]], pch=16)
                    Sys.sleep(.1)
          }
          par(old)
}


## grColExt used in 'animaColExt' and 'bioGeoIsl'
grColExt=function(E , I , P, area)
{
	# Following code avoids spurious NOTE by R CMD check:
	x <- NULL; rm(x);

	S = I*P/(I+E) ; T = I*E/(I+E)
	nIsl=length(E)
	corIsl=rainbow(nIsl)
	curve(I[1]-I[1]*x/P[1],0,P[1],bty="n",xlab="Number of Species", ylab="Rate",xaxt="n",yaxt="n", font.lab=2,lwd=2,ylim=c(0,1),  main="Island Biogeography", col=corIsl[1])
	curve((E[1]/P[1])*x,0,P,lwd=2,add=TRUE, col=corIsl[1], lty=2)
	legend("top", legend=c("Colonization", "Extinction"),  bty="n",lty=c(1,2))
	abline(v=0)
	abline(h=0)
	mtext("St",side=1,at=P,font=2, line=1)
	linhas=seq(0,1.5, length.out=nIsl)
	for(i in 1:nIsl)	
	{
		curve(I[i]-I[i]*x/P,0,P,lwd=2,add=TRUE, col=corIsl[i], lty=1)
		curve((E[i]/P)*x,0,P,lwd=2,add=TRUE, col=corIsl[i], lty=2)
		mtext(paste("S", i, sep=""),side=1,at=S[i], cex=0.8,font=2,col=corIsl[i], line=linhas[i])
		mtext(paste("T", i, sep=""),side=2,at=T[i],cex=0.8,font=2,las=1,col=corIsl[i], line=linhas[i])
		points(S[i],T[i],col=corIsl[i],pch=16,cex=1)
		if(length(unique(area))>1)
		{
			siz_ar=2 +(area/max(area))
			points(S[i],T[i],col=corIsl[i],cex=siz_ar[i])
		}
		segments(S[i],T[i],S[i],0,lty=3,col=corIsl[i])
		segments(S[i],T[i],0,T[i],lty=3,col=corIsl[i])
		Sys.sleep(0.1)
	}	
}

#grColExt(E = .5 , I = .5 , P = 100, area=1:10)


## animaRandWalk used in 'randWalk'
animaRandWalk = function(rwData, time=2, sleep=0.1)
{
          xplus=max(time)*0.1
          ymax=max(apply(rwData, 2, max))[1]
          plot(time, rwData[,which.max(apply(rwData, 2, max))[1]], xlab="Steps", ylab="Distance from the edge",cex.axis=1.2, cex.lab=1.2,ylim=c(-.1* ymax,ymax), main="Random Walk", cex.main=1.5, type="n", xlim=c(0,max(time)))
          polygon(x=c(-xplus, -xplus, max(time)+xplus, max(time)+xplus), y=c(ymax*-0.15,0,0,ymax*-0.15), col="gray")
          text(max(time)/2, -0.05* ymax, labels="Absorption Surface", col="red", cex=1.5)
          n=dim(rwData)[2]
          ncolors= rainbow(n)
          for(i in 2:length(time))
          {
                    for(j in 1:n)
                    {
                              lines(time[1:i], rwData[1:i,j], col=ncolors[j], lty=j )
                    }
                    Sys.sleep(sleep)
          }
}


## animaGame used in 'extGame'
animaGame = function(xGame, total, sleep=0.01)
{
          xmax=length(xGame)
          xseq=1:xmax
          if(xmax>1e3){sleep=0}
          if(xmax>1e4)
          {
                    indx=ceiling(seq(1,xmax, len=1000)) 
                    xGame=xGame[indx]
                    xseq=xseq[indx]
          }
          plot(0:xmax, seq(0,total, len=xmax+1), xlab="Cycle", ylab="Money amount",cex.axis=1.2, cex.lab=1.2, ylim=c(-.1* total,total+total*0.1), main="Zero Sum Game", cex.main=1.5, type="n", cex.sub=0.9)
          abline(h=total/2, lty=2, col="red")
          cores= c("blue","black")
          for(i in 2:xmax)
          {
                    lines(xseq[1:i], xGame[1:i], col=cores[1], lty=2)
                    lines(xseq[1:i], total - xGame[1:i], col=cores[2], lty=3)
                    Sys.sleep(sleep)
          }
          polygon(x=c(-.2* xmax, -.2* xmax, xmax+ 0.1*xmax, xmax+ 0.1*xmax), y=c(-.2*total,0,0,-.2* total), col="gray")
          polygon(x=c(-.2*xmax, -.2*xmax, xmax+ 0.1*xmax, xmax+ 0.1*xmax), y=c(total,total+total*.5,total +total*.5,total), col="gray")
          text(xmax/2, - 0.05* total, labels="Loser", col="red", cex=1.5)
          text(xmax/2, total + 0.05* total, labels="Winner", col="green", cex=1.5)
}


## animaHub used in 'simHub1', 'simHub2', 'simHub3'
animaHub=function(dadoHub, sleep=0.1)
{
          maxsp=max(dadoHub)[1]
          uniqsp=unique(as.numeric(dadoHub))
          nind=dim(dadoHub)[1]
          nsim=dim(dadoHub)[2]
          ciclo=as.numeric(colnames(dadoHub))
          pb = tkProgressBar(title = "Simulation Progress", max = nsim)
          riq=apply(dadoHub, 2, rich)
          ## definindo o tamanho do retangulo
          lado<-round(sqrt(nind))
          lado2<-ceiling(nind/lado)
          lastLine=lado*lado2 - nind
          cormix=sample(rainbow(maxsp+10))
          cor=c("#FFFFFF", cormix)
          mcor<-c("#FFFFFF00","#000000")
          spcol<-c(rep(0, lastLine),dadoHub[,1])
          ############ escala das especies da metapopulacao ########
          layout(matrix(data=c(2,1), nrow=2, ncol=1), widths=c(1,1), heights=c(5,1))
          old<-par(mar=c(2,2,1,2))
          image(x=1:maxsp, y=1, matrix(data=1:maxsp, nrow=maxsp,ncol=1),col=rainbow(maxsp), ylab="",xlab="", xaxt="n", yaxt="n", main="Metacommunity Species colors", cex.main=0.8)
          axis(3, at = c(1,maxsp), labels = c(1, maxsp), tick = FALSE, mgp=c(1,0,0), cex.axis=0.8)
          hmat=matrix(spcol,ncol=lado, nrow=lado2)
          par(mar=c(2,2,2,2))
          image(hmat, col=cor[sort(unique(as.numeric(hmat)))], xaxt="n", yaxt="n")
          grid(nx=lado2, ny=lado)
          for (i in 2:nsim) 
          {
                    mvf=dadoHub[,i-1]!=dadoHub[,i]
                    matm<-matrix(c(rep(FALSE, lastLine),mvf ),ncol=lado, nrow=lado2)
                    image(matm,col=mcor, add=TRUE)
                    Sys.sleep(sleep)
                    spcol<-c(rep(0, lastLine),dadoHub[,i] )
                    cores=cor[sort(unique(spcol)+1)]
                    scol<-sort(unique(spcol))
                    lcol<-length(scol)
                    mcol<-match(spcol, scol)
                    hmat=(matrix(mcol,ncol=lado, nrow=lado2))
                    image(hmat, col=cores, add=TRUE)
                    grid(nx=lado2, ny=lado)
                    setTkProgressBar(pb, value = i, label = paste("Simulation #", ciclo[i], sep="")) 
          }
          close(pb)
}



###############################
### Two Species Competition ###
###############################

### Meta competition 'animaMetaComp' used in 'metaComp'
animaMetaComp=function(dados)
{
nsim=dim(dados)[3]
ln=dim(dados)[1]
cl=dim(dados)[2]
op=par(mar=c(1,2,2,2))
layout(matrix(c(2,1), ncol=1, nrow=2), heights=c(5,1),widths=c(1,1))
plot(1:10,1:10,xaxt="n", yaxt="n", xlab="", ylab="", cex=0.8,type="n", , bty="n")
legend(0.5,11.8,ncol=4, legend=c("not available", "empty", "sup. competitor", "inf. competitor"), pch=c(15,22,15,15), title="Patches legend", col=c("red","black", "blue", "green"),bty="n")
image(0:ln, 0:cl, dados[,,1], col=c("red", "white","blue" ,"green") , breaks=c(-0.9,-0.001,0.1,1.5,2.9),main="Metapopulations Competition",  xlab="", ylab="")
grid(ln,cl)
Sys.sleep(.5)
	for(i in 2:nsim)
	{
	par(new=TRUE)
image(0:ln, 0:cl, dados[,,i], col=c("red", "white","blue" ,"green") , breaks=c(-0.1,-0.001,0.1,1.9,2.9), xlab="", ylab="")
grid(ln,cl)
	Sys.sleep(.1)
	}
}



#######################
### Metapopulations ###
#######################

## animaMeta2 used in 'metaPop', 'metaEr', 'metaCi', 'metaCier'
animaMeta2=function(dados)
{
nsim=dim(dados)[3]
ln=dim(dados)[1]
cl=dim(dados)[2]
image(0:ln, 0:cl, dados[,,1], col=c("white", "green") , breaks=c(0,0.99,5),main="Metapopulation Dynamics", sub=paste("Initial configuration from", nsim," simulations",  sep=""), xlab="", ylab="")	
grid(ln,cl)
Sys.sleep(.5)
conta12=dados[,,1]+ (2*dados[,,2])
image(0:ln, 0:cl, conta12, col=c("white","red","lightgreen", "darkgreen") , breaks=c(0,0.9,1.9,2.9,3.9),main="Metapopulation Dynamics", sub=paste("red= extinction; light green= colonization; dark green = permanence \n maximum time = ", nsim, sep=""), xlab="", ylab="")
	for(i in 3:nsim)
	{
	conta12=dados[,,(i-1)]+ (2*dados[,,i])
	image(0:ln, 0:cl, conta12, col=c("white","red","lightgreen", "darkgreen") , breaks=c(0,0.9,1.9,2.9,3.9), xlab="", ylab="", add=TRUE)
	Sys.sleep(.1)
	}
}


## grFim used in 'metaPop', 'metaEr', 'metaCi', 'metaCier'
grFim=function(dados)
{
op=par(mfrow=c(2,2))
nsim=dim(dados)[3]
ln=dim(dados)[1]
cl=dim(dados)[2]
image(0:ln, 0:cl, dados[,,1], col=c("white", "green") , breaks=c(0,0.99,5),main="Metapopulation Dynamics", sub=paste("time = 1/", nsim, sep=""), xlab="", ylab="")	
grid(ln,cl)
	for(ts in c(4,2,1))
	{
	sim=round(nsim/ts)
	conta12=dados[,,(sim-1)]+ (2*dados[,,sim])
	image(0:ln, 0:cl, conta12, col=c("white","red","lightgreen", "darkgreen") , breaks=c(0,0.9,1.9,2.9,3.9),main="Metapopulation Dynamics", sub=paste("red= extinction; light green= colonization;\n dark green = permanence \t time = ", sim, "/", nsim, sep=""), xlab="", ylab="")
	}
par(op)
}
##############################
### Multispecies Functions ###
##############################
###############################
#Trade-off Multispecies Graphic
### 
gr.toff=function(rq, fsp1,pe,add=FALSE,...)
{
#	rq <- as.numeric(tclvalue(rqVar))
#	fsp1 <- as.numeric(tclvalue(fsp1Var))
#	pe <- as.numeric(tclvalue(peVar))
	rank=1:rq
	ci= pe/(1-fsp1)^(2*rank-1)
	px= fsp1*(1-fsp1)^(rank-1)
		if(add==FALSE)
		{
		toff<-dev.new( width=5, height=5)
		}
	old<-par(mar=c(3,3,3,3))
	plot(ci~rank,col="red",ylim=c(0,max(ci)*1.1), type="b", ann=FALSE, axes=FALSE)
	axis(4, cex.axis=0.8)#, yaxp=c(0,3,3))
	par(new=TRUE)
	plot(px~rank, ylim=c(0,fsp1),type="b", bty="n",  ann=FALSE, cex.axis=0.8)#yaxt="n", xaxp=c(0,10,5))
	#axis(2, cex.axis=0.8)#, yaxp=c(0,0.2,4))
	mtext("Species competitive rank", 1, 2, cex=0.9)
	mtext("Abundance", 2, 2, cex=0.9)
	mtext("Colonization rate", 4, 2, cex=0.9)
	mtext("Trade-off Species Rank ", 3, 0, cex=1.2)
	par(old)
}
###############################
## animaCena used in 'regNicho'
animaCena=function(dados)
{
nt=dim(dados)[3]
dev.new()
op=par(mfrow=c(5,5),  mar=c(0.1,0.1,0.1,0.1))
	for(i in 1:nt)
	{
	image(dados[,,i], main="",  bty="n",xaxt='n',yaxt='n', col=c("white", "yellow", "orange", "blue", "green"))
	grid(dim(dados)[2],dim(dados)[1])
	}
dev.new()
par(mfrow=c(2,2))
image(dados[,,1], main= paste("Patch occupancy\n \t time=", 1 ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
image(dados[,,round(nt/3)], main= paste("Patch occupancy\n \t time=", round(nt/3) ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
image(dados[,,round(2*nt/3)], main= paste("Patch occupancy\n \ttime=", round(2*nt/3) ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
image(dados[,,nt], main= paste("Patch occupancy\n \t time=", nt ),  bty="n",xaxt='n',yaxt='n',col=c("white", "yellow", "orange", "blue", "green"))
grid(dim(dados)[2],dim(dados)[1])
par(op)
}
