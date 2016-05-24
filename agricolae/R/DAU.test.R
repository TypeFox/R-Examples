`DAU.test` <-
function (block, trt, y, method = c("lsd","tukey"),alpha=0.05,group=TRUE,console=FALSE)
{
    method<-match.arg(method)
	if(method =="lsd") snk=5
	if(method =="tukey") snk=6
    block.unadj <- as.factor(block)
    trt.adj <- as.factor(trt)
    block.adj <- as.factor(block)
    trt.unadj <- as.factor(trt)
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    modelo1 <- formula(paste(name.y,"~ block.unadj + trt.adj"))
    model1 <- lm(modelo1)
    glerror <- df.residual(model1)
    MSerror <- deviance(model1)/glerror
    modelo2 <- formula(paste(name.y,"~ trt.unadj + block.adj"))
    model2 <- lm(modelo2)
    r <- unique(table(trt.adj))
    b <- nlevels(block.unadj)
    ntr <- nlevels(trt.adj)  # nro trt
    Means<-tapply.stat(y,trt, function(x) mean(x,na.rm=TRUE))
	mi <- tapply.stat(y, trt, function(x) min(x,na.rm=TRUE))
	ma <- tapply.stat(y, trt, function(x) max(x,na.rm=TRUE))
	n.rep <- tapply.stat(y, trt, function(x) length(na.omit(x)))
	sds<- tapply.stat(y, trt, function(x) sd(x,na.rm=TRUE))
	std.err<-sds[,2]/sqrt(n.rep[,2])
	Means<-data.frame(Means,std=sds[,2],r=n.rep[,2],Min=mi[,2],Max=ma[,2])
	names(Means)[1:2] <- c(name.t, name.y)
	mean.y<- mean( y, na.rm=TRUE )
	nameTrt<-Means[,1] 
    #mean.block<- data.frame(mean.block,ri=mean.block[,2]-mean.y )
    r.trt<-tapply.stat( y, trt, length)
    names(r.trt)[2]="N"
    r.trt<-data.frame(r.trt,control=r.trt[,2]==b,means=Means[,2],mean.adj =Means[,2], 
			block="",std.err=sqrt(MSerror/r.trt[,2]))
    r.trt[,6]<-as.character(r.trt[,6])
    estado <- NULL 
    n<-length(y)
	for (i in 1:n) {
    for (j in 1:ntr) {
    if (trt[i] == r.trt[j,1]) {
    estado[i]<-r.trt[j,3];
    break;
    }
    }}
    datos <- data.frame(block, trt, y, control=estado)
    mean.block<-tapply.stat( datos[datos$control,3],datos[datos$control,1], function(x) mean(x,na.rm=TRUE))
    mean.y<- mean( datos[datos$control,3], na.rm=TRUE )
    mean.yy<-mean( datos[,3], na.rm=TRUE )
    mean.block<- data.frame(mean.block,ri=mean.block[,2]-mean.y )
    tc<-sum(datos[,3])^2/n
    # treatment unadjusted
    # comunes vs aumentando
    sc1<-sum(datos[datos$control,3])^2/length(datos[datos$control,3]) + 
    sum(datos[!datos$control,3])^2/length(datos[!datos$control,3]) - tc
    # entre comunes
    unadj <- datos[datos$control,1:3]
    sc2 <- anova(lm(unadj[,3] ~ unadj[,2]))[1,2]
    # entre aumentados
    sc3 <- anova(model2)[1,2] - sc1 - sc2
    # treatment adjusted    
    # Comunes + comunes vs aumentados
    A1 <- anova(model1)
    sc4<- A1[2,2] - sc2 
    A1 <- rbind(A1,A1[3,],A1[3,])
    rownames(A1)[3:5]<-c("Control","Control + control.VS.aug.","Residuals")
    A1[3,1]<-length(r.trt[r.trt$control,3]) - 1
    A1[4,1]<-ntr-length(r.trt[r.trt$control,3])
    A1[3,2]<-sc2
    A1[4,2]<-sc4
    A1[3,3]<-A1[3,2]/A1[3,1]
    A1[4,3]<-A1[4,2]/A1[4,1]
    A1[3,4]<-A1[3,3]/A1[5,3]
    A1[4,4]<-A1[4,3]/A1[5,3]
    A1[3,5]<-1-pf(A1[3,4],A1[3,1],A1[5,1])
    A1[4,5]<-1-pf(A1[4,4],A1[4,1],A1[5,1])
    A1[1,4]<-NA
    A1[1,5]<-NA    
    A2<- anova(model2)
    A2<- rbind(A2,A2[3,],A2[3,],A2[3,])
    rownames(A2)[3:6]<-c("Control","Augmented","Control vs augmented","Residuals")
    A2[3,1]<-length(r.trt[r.trt$control,3]) - 1
    A2[4,1]<-length(r.trt[!r.trt$control,3]) - 1
    A2[5,1]<-1
    A2[3,2]<-sc2
    A2[4,2]<-sc3
    A2[5,2]<-sc1
    A2[3,3]<-A2[3,2]/A2[3,1]
    A2[4,3]<-A2[4,2]/A2[4,1]
    A2[5,3]<-A2[5,2]/A2[5,1]
    A2[3,4]<-A2[3,3]/A2[6,3]
    A2[4,4]<-A2[4,3]/A2[6,3]
    A2[5,4]<-A2[5,3]/A2[6,3]
    A2[3,5]<-1-pf(A2[3,4],A2[3,1],A2[6,1])
    A2[4,5]<-1-pf(A2[4,4],A2[4,1],A2[6,1])        
    A2[5,5]<-1-pf(A2[5,4],A2[5,1],A2[6,1])
    A2[1,4]<-NA
    A2[1,5]<-NA
# Calcula mean adj y ubica el bloque
    for(i in 1:ntr) {
    if(!r.trt[i,3]) {
    for(j in 1:n) {
    if(datos[j,2]== r.trt[i,1] ) {
    marca<-as.character(datos[j,1])
    break
    }
    }
    for(l in 1:b) {
    if(mean.block[l,1]==marca) {
    r.trt[i,6] <- marca
    r.trt[i,5] <- r.trt[i,5]- mean.block[l,3]
    break
    }
    }
    }
    }
	n.rep
# matriz variancia
comunes<-sum(r.trt[,3])
V<-diag(0,ntr)
rownames(V)<-r.trt[,1]
colnames(V)<-r.trt[,1]  
for(i in 1:ntr){
for(j in 1:ntr) {
if(i!=j) {
if (r.trt[i,3]& r.trt[j,3])   V[i,j]<- 2/b
if (!r.trt[i,3]&!r.trt[j,3])  {
    if(r.trt[i,6]==r.trt[j,6] ) V[i,j]<- 2
    else  V[i,j]<- 2*(1+1/comunes)
    }
if (r.trt[i,3] != r.trt[j,3]) V[i,j]<- 1+1/b+1/comunes-1/(b*comunes)

}}}
V<-MSerror*V
#
CV<-round(cv.model(model1), 1)
if(console){
cat("\nANALYSIS DAU: ", name.y, "\nClass level information\n")
    cat("\nBlock: ", unique(as.character(block)))
    cat("\nTrt  : ", as.character(r.trt[,1]))
    cat("\n\nNumber of observations: ", length(y), "\n")
    cat("\nANOVA, Treatment Adjusted\n")
    print(A1)
    cat("\nANOVA, Block Adjusted\n")
    print(A2)
    cat("\ncoefficient of variation:", CV, "%\n")
    cat(name.y, "Means:", mean.yy, "\n")
}
 SE.dif<-rbind(sqrt(2*MSerror/b),sqrt(2*MSerror),sqrt(2*MSerror*(1+1/comunes)),
 sqrt(MSerror*(1+1/b+1/comunes-1/(b*comunes))))
rownames(SE.dif)<- c("Two Control Treatments","Two Augmented Treatments (Same Block)",
"Two Augmented Treatments(Different Blocks)","A Augmented Treatment and A Control Treatment")
colnames(SE.dif)<-"Std Error Diff."
#    cat("\nCritical Differences (Between)         Std Error Diff.\n")
if(console){cat("\nCritical Differences (Between)\n")
print(SE.dif)
cat("\n")}   
if(!group){
#	Omeans<-order(mean.adj,decreasing = TRUE)
#	Ordindex<-order(Omeans)
	comb <- utils::combn(ntr, 2)
	nn <- ncol(comb)
	dif <- rep(0, nn)
	sig <- rep(" ",nn)
	pvalue <- dif
	odif<-dif
	for (k in 1:nn) {
		i <- comb[1, k]
		j <- comb[2, k]
		dif[k] <- r.trt[i, 5] - r.trt[j, 5]
		tc <- abs(dif[k])/sqrt(V[i,j])
    if (method == "lsd")
    pvalue[k] <- 2 * round(1 - pt(tc, glerror), 4)
    if (method == "tukey")
    pvalue[k] <- round(1 - ptukey(tc*sqrt(2), ntr, glerror), 4)
	sig[k]<-" "
	if (pvalue[k] <= 0.001) sig[k]<-"***"
	else  if (pvalue[k] <= 0.01) sig[k]<-"**"
	else  if (pvalue[k] <= 0.05) sig[k]<-"*"
	else  if (pvalue[k] <= 0.1) sig[k]<-"."
    }
	groups <- NULL
	tr.i <- nameTrt[comb[1, ]]
	tr.j <- nameTrt[comb[2, ]]
	comparison<-data.frame("Difference" = dif, pvalue=pvalue,"sig."=sig)
	rownames(comparison)<-paste(tr.i,tr.j,sep=" - ")
    }
	if(group){
	if(console){
		cat("\n\nMeans with the same letter are not significantly different.")
		cat("\n\nGroups, Treatments and means\n")}
		groups <- order.group(trt=r.trt[,1], r.trt[,5], r.trt[,2], MSerror=NULL, Tprob=NULL, 
		std.err=r.trt[,7],parameter=1,snk, DFerror=glerror,alpha,sdtdif=1, vartau=V,console=console)
		names(groups)[2]<-"mean.adj"
		groups<-groups[,1:3]
		comparison=NULL
	}
    if(console){cat("\nComparison between treatments means\n")
    cat("\n<<< to see the objects: comparison and means  >>>\n\n")}
    #pvalue <- as.dist(pvalue)
    #    N = r, std.err = sqrt(diag(vartau)))
    means<-data.frame(Means[,1:6],mean.adj=r.trt[,5],SE=r.trt[,7],block=r.trt[,6])
    rownames(means)<- means[,1]
    means<-means[,-1]
	parameters<-data.frame(treatments=ntr,Controls=comunes,
	Augmented=ntr-comunes,blocks=b,alpha=alpha,test="DAU")
	statistics<-data.frame(Mean=mean.yy,CV=CV)
	rownames(parameters)<-" "
	rownames(statistics)<-" "
	output<-list(means = means, parameters=parameters, statistics=statistics,
	comparison=comparison,groups=groups,SE.difference=SE.dif,vartau = V)
	invisible(output)
}


