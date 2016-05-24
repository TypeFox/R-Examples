summary.sparsereg<-function(object,...) {
summary.inner.sparsereg <- function(object,interval=.9,ci="quantile",order="sparse",
normal=TRUE,select="mode",printit=TRUE,stage=NULL,... ){

if(object$modeltype=="onestage") stage<-1
if(object$modeltype=="twostage"&length(stage)==0) stop("Please specify the stage for which you want results")
		if(stage==2){
			object$beta.ci<-object$beta.ci_2
			object$beta.mode<-object$beta.mode_2
			object$beta.mean<-object$beta.mean_2			
		}

if(normal) ests<-object$beta.ci else ests<-object$beta.mode
if(select=="mode") est.obj<-object$beta.mode else est.obj<-object$beta.mean
point.ests<-apply(object$beta.mode,2,median)
int.calc<-sort(c((1-interval)/2,1-(1-interval)/2))
if(ci=="quantile") quants.out<-cbind(apply(ests,2,quantile,min(int.calc)),
						apply(ests,2,quantile,max(int.calc)) )
if(ci=="HPD") quants.out<-HPDinterval(ests,interval)
p.vals<-colMeans(object$beta.mode!=0)

sum.mat<-cbind(colSums(ests<0),colSums(ests==0),colSums(ests>0))
prob.mat<-t(apply(sum.mat,1,FUN=function(x) x/sum(x)))
if(select=="mode") rows.keep<-(point.ests!=0)
if(select=="all")  rows.keep<-(1:length(point.ests))
if(is.numeric(select)) rows.keep<-(p.vals>select)

table.init<-as.matrix(cbind(point.ests,quants.out,prob.mat))
table.out<-as.matrix(table.init[rows.keep,],ncol=6)
if(sum(rows.keep)==1) table.out<-t(table.out)
colnames(table.init)<-colnames(table.out)<-c("Posterior Median",paste(round(int.calc[1]*100,2),"%",sep=""),paste(round(int.calc[2]*100,2),"%",sep=""),"Pr(b<0)","Pr(b=0)","Pr(b>0)")

if(order=="all") order.rows<-1:nrow(table.init)
if(order=="sparse") order.rows<-which(point.ests!=0)
if(order=="magnitude") {
	order.rows<-sort(abs(point.ests),index.return=TRUE,decreasing=TRUE)$ix
	order.rows<-order.rows[abs(point.ests[order.rows])>0]
	}
if(order=="alphabetical") {
	order.rows<-sort(rownames(table.init),index.return=TRUE)$ix
	order.rows<-order.rows[abs(point.ests[order.rows])>0]
	}	
if(order=="alphabeticalall") {
	order.rows<-sort(rownames(table.init),index.return=TRUE)$ix
	}	
if(!printit){
	out<-list("table"=table.init)
return(out)

	
}


if(printit){
cat('LASSOplus results:')
cat('\n----------------------------')
cat('\nVariable Selection: ')
cat(paste('\n   Original variables:', ncol(object$X)))
cat(paste('\n   Selected variables:', sum(point.ests!=0)))
cat('\n----------------------------')
cat('\nCoefficents:\n ')
print(signif(table.init[order.rows,],3))
cat('----------------------------')
cat('\n Posterior intervals using ')
cat(ifelse(ci=="quantile","quantiles of ","HPD interval of "))
cat(ifelse(normal,"the approximate confidence interval","the posterior density"))
out<-list("table"=table.init)
invisible(out)

}

}




summary.inner.sparsereg(object,...)
}


plot.sparsereg<-function(x,...){
plot.sparsereg.inner<-function(x, main1="Main Effects", main2="Interaction Effects", main3="Three-way Interaction Effects",xlabel="Effect", plot.one=FALSE,stage,...) {
if(x$modeltype=="onestage") stage<-1
if(x$modeltype=="twostage"&stage!=1&stage!=2) stop("Please specify the stage for which you want results")
		if(stage==2){
			x$beta.ci<-x$beta.ci_2
			x$beta.mode<-x$beta.mode_2
			x$beta.mean<-x$beta.mean_2			
		}
x$modeltype<-"onestage"

names.est<-colnames(x$beta.mode)
sum1<-summary(x,printit=FALSE)$table
inter.count<-unlist(lapply(sapply(names.est,strsplit," x "),length))

length.names<-sapply(names.est,nchar)
#for(i in 1:length(names.est)) for(j in length.names[i]:max(length.names)) names.est[i]<-paste(" ",names.est[i],sep="")
plot.mat<-as.matrix(x$beta.mode)
num.mat<-name.mat<-plot.mat

for(i in 1:ncol(plot.mat)) {
	name.mat[,i]<-colnames(plot.mat)[i]
	num.mat[,i]<-i
	}
plot.vec<-as.vector(plot.mat)
num.vec<-as.vector(num.mat)
name.vec<-as.vector(name.mat)
plotobj<-data.frame(names.est,sum1[,1:3])
names(plotobj)<-c("names","median50","lo","hi")
drop.zero<-(plotobj[,2]!=0)


d.nz<-plotobj[drop.zero,]
inter.nz<-inter.count[drop.zero]
lo<-hi<-median50<-NULL

if(sum(inter.nz==1)>0){
g1<-ggplot(d.nz[inter.nz==1,], aes(x=names, y=median50))+
	geom_errorbar(aes(ymin=lo,ymax=hi),width=.1)+coord_flip()+ theme(legend.position="none")+ylab(xlabel)+ theme(axis.title.y = element_blank())+geom_point(shape="|")+
ggtitle(main1)+ylim(range(d.nz[,-1]))
}


if(sum(inter.nz==2)>0){
g2<-ggplot(d.nz[inter.nz==2,], aes(x=names, y=median50))+
geom_errorbar(aes(ymin=lo,ymax=hi),width=.1)+coord_flip()+ theme(legend.position="none")+ylab(xlabel)+ theme(axis.title.y = element_blank())+geom_point(shape="|")+
ggtitle(main2)+ylim(range(d.nz[,-1]))
}


if(sum(inter.nz==3)>0){
g3<-ggplot(d.nz[inter.nz==3,], aes(x=names, y=median50))+
geom_errorbar(aes(ymin=lo,ymax=hi),width=.1)+coord_flip()+ theme(legend.position="none")+ylab(xlabel)+ ggtitle(main3)+theme(axis.title.y = element_blank())+geom_point(shape="|")+ylim(range(d.nz[,-1]))
}

if(plot.one==1 & sum(inter.nz==1)==0){ stop("You do not have any main effects to plot.  Check the value of plot.one.")}
if(plot.one==2 & sum(inter.nz==2)==0){ stop("You do not have any interaction effects to plot. Check the value of plot.one.")}
if(plot.one==3 & sum(inter.nz==3)==0){ stop("You do not have any two-way interaction effects to plot. Check the value of plot.one.")}


plot.two<-function(g1,g2){
gp1<- ggplot_gtable(ggplot_build(g1))##from grid.table
gp2<- ggplot_gtable(ggplot_build(g2))
maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3])
gp1$heights<-gp1$heights
gp2$heights<-gp2$heights

gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth
grid.arrange(gp1,gp2, heights=c(sum(inter.nz==1),sum(inter.nz==2)))
}

plot.three<-function(g1,g2,g3){
gp1<- ggplot_gtable(ggplot_build(g1))##from grid.table
gp2<- ggplot_gtable(ggplot_build(g2))
gp3<- ggplot_gtable(ggplot_build(g3))
maxWidth = unit.pmax(gp1$widths[2:3], gp2$widths[2:3], gp3$widths[2:3])
gp1$widths[2:3] <- maxWidth
gp2$widths[2:3] <- maxWidth
gp3$widths[2:3] <- maxWidth
grid.arrange(gp1,gp2, heights=c(sum(inter.nz==1),sum(inter.nz==2),sum(inter.nz==3)))##Uses gridExtra
}

if(plot.one==FALSE){
if(mean(inter.nz==1)==1) plot(g1)
if(mean(inter.nz==2)==1) plot(g2)
if(mean(inter.nz==3)==1) plot(g3)
if(length(unique(inter.nz))==2){
	if(mean(inter.nz%in%c(1,2))==1) plot.two(g1,g2)
	if(mean(inter.nz%in%c(1,3))==1) plot.two(g1,g3)
	if(mean(inter.nz%in%c(2,3))==1) plot.two(g2,g3)
}
if(length(unique(inter.nz))==3) plot.three(g1,g2,g3)
}

if(plot.one==1) plot(g1)
if(plot.one==2) plot(g2)
if(plot.one==3) plot(g3)


}



plot.sparsereg.inner(x,...)
}
#plots the densities of the posterior

violinplot<-function(x,columns=NULL,newlabels=NULL,type="mode",stage=NULL){
	if(x$modeltype=="onestage") stage<-1
		if(stage!=1 & stage!=2 & x$modeltype=="twostage") stop("For a two stage object, stage must be set to 1 or 2")

		if(stage==2){
			x$beta.ci<-x$beta.ci_2
			x$beta.mode<-x$beta.mode_2
			x$beta.mean<-x$beta.mean_2			
		}

    input<-x
	if(type=="mode") m<-input$beta.mode
	if(type=="mean") m<-input$beta.mean
	if(type=="ci") m<-input$beta.ci

	if(length(columns)==0) columns<-which(apply(input$beta.mode,2,median)!=0)
	if(length(columns)==0) stop("No covariates were selected.")
	num.plot<-name.plot<-m.plot<-as.matrix(m[,columns])
	if(is.numeric(columns))	colnames.fill<-colnames(m)[columns]
	if(is.character(columns))	colnames.fill<-colnames(m)[colnames(m)%in%columns]

	if(length(newlabels)>0) colnames.fill<-newlabels
	for(i in 1:length(columns)) {
		num.plot[,i]<-i
		name.plot[,i]<-colnames.fill[i]
	}
	data.plot<-data.frame(
	as.vector(name.plot),as.vector(num.plot),as.vector(m.plot)
	)
	names(data.plot)<-c("names","number","value")
	..density..<-value<-data<-number<-names<-NULL
	
		ylim.plot<-range(data.plot$value)

ggplot(data.plot,aes(x=value)) +
  stat_density(aes(ymax = ..density..,  ymin = -..density..),
    fill = "grey50", colour = "grey50",
    geom = "ribbon", position = "identity") +
  facet_grid(. ~ names) +
  coord_flip()+xlim(ylim.plot)

}


difference<-function(x,type="mode",var1=NULL,var2=NULL,plot.it=TRUE, main="Difference",xlabel="Effect",
ylabel="Density" ){
  input<-x
	if(type=="mode") out<-input$beta.mode
	if(type=="mean") out<-input$beta.mean
	if(type=="ci")  out<-input$beta.ci
  	out<-as.matrix(out)
  #rename the variables names some for better formatting and synching with summary.sparsereg etc..  
  colnames(out)<-gsub("_",": ",colnames(out))

  #Select the variables
  if(class(var1)=="character") t1<-which(colnames(out)==var1)
  if(class(var2)=="character")   t2<-which(colnames(out)==var2)
  p1<-out[,var1]
  p2<-out[,var2]
  #take the different t1-t2
  diffpt<-p1-p2
  m <- ggplot(data.frame(diffpt), aes(x = diffpt))+ geom_density()+ggtitle(main)+xlab(xlabel)+ylab(ylabel)
	if(plot.it) plot(m)
  #display quantiles
  print(quantile(diffpt))
}

print.sparsereg<-function(x,...) {
	print.sparsereg.inner<-function(x,stage=NULL,...){
	if(x$modeltype=="onestage") stage<-1
	if(x$modeltype=="twostage"&length(stage)==0) stop("Please specify the stage for which you want results")
		if(stage==2){
			x$beta.ci<-x$beta.ci_2
			x$beta.mode<-x$beta.mode_2
			x$beta.mean<-x$beta.mean_2			
		}
	
	x$beta.mode<-as.mcmc(x$beta.mode)
	summary(x$beta.mode)
	}
	return(print.sparsereg.inner(x,...))

}

