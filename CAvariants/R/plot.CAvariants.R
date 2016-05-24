plot.CAvariants<-function(x,  firstaxis=1, lastaxis=2, cex=0.8,cex.lab=0.8,prop=1,
plottype="biplot", biptype = "row",M=x$r,scaleplot=1,
posleg="topleft",pos=2,ell=FALSE,alpha=0.05,...) {

## internal function to plot a  single picture
##

if ((firstaxis<1)|(firstaxis>x$maxaxes-1)) stop(paste("incorrect first axis =", firstaxis, "\n\n")) 
if (lastaxis>x$maxaxes) stop(paste("incorrect last axis =", lastaxis, "\n\n")) 
if (firstaxis>=lastaxis) stop(paste("last axis must be greater than first axis\n\n"))
#if (!any(plottype==c("classic","classical","c","biplot","bip","b"))) stop(paste("Must be specified the kind of graph: classic, or biplot"))

# Groups file must have no blank line at start and only one between sections
# group number   group name   symbol   colour   plot ellipse? 
n<-sum(x$DataMatrix)
I<-nrow(x$DataMatrix)
rowgroup <- list(1:x$rows,rep(1,x$rows))
rowgrlab <- list(1,"","*","red","T")
colgroup <- list(1:x$cols,rep(1,x$cols))
colgrlab <- list(1,"","+","blue","T")
######################################################
# Plot row and col coordinates
#########################################
if ((plottype=="Classical")|(plottype=="classical")|(plottype=="classic")|(plottype=="c")) {
nthings<-x$cols
nvars<-x$rows
cord1<- x$Cprinccoord 
cord2<-x$Rprinccoord
dmu=diag(x$inertias[,1])
inertiapc=x$inertias[,2] #inertia in percentage of row axes 
dimnames(cord1)[1]<-dimnames(x$DataMatrix)[2]
dimnames(cord2)[1]<-dimnames(x$DataMatrix)[1]
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
thinglabels<-x$collabels
varlabels<-x$rowlabels
main="Classical plot"
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
cat("\n ERROR: NO CLASSICAL PLOT for ordered analysis. ONLY A BIPLOT can be constructed  (Please change 'plottype' and specify 'biptype')\n")
stop()
#cord1<-x$Cstdcoord/scaleplot
#cord2<-x$Rprinccoord*scaleplot
#inertiapc=x$inertias2[,2] #inertia in percentage of row axes 
#ell=ell
#--------news
}
}#end classical plot
         if ((plottype=="Biplot")|(plottype=="biplot")|(plottype=="bip")|(plottype=="b")){
if ((biptype=="rows")|(biptype=="Rows")|(biptype=="row")|(biptype=="r")) 
{
plottype<-"biplot"
#biptype<-"row"
#if (ell==TRUE) { scaleplot<-1
#}
cord1<-x$Rprinccoord*scaleplot
cord2<-x$Cstdcoord/scaleplot
nthings<-x$rows
nvars<-x$cols
thinglabels<-x$rowlabels
varlabels<-x$collabels
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
main<-"Row Isometric Biplot"
inertiapc=x$inertias[,2] #inertia of column poly
dmu=diag(x$inertias[,1])
dimnames(cord2)[1]<-dimnames(x$DataMatrix)[2]
dimnames(cord1)[1]<-dimnames(x$DataMatrix)[1]
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
#if (ell==TRUE) { scaleplot<-1
#}
cord2<-x$Rprinccoord*scaleplot
cord1<-x$Cstdcoord/scaleplot
nthings<-x$cols
nvars<-x$rows
thinglabels<-x$collabels
varlabels<-x$rowlabels
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-rowgroup
vargrlab<-rowgrlab
inertiapc=x$inertias2[,2] #inertia of column poly
dmu=diag(x$inertias2[,1])
dimnames(cord2)[1]<-dimnames(x$DataMatrix)[1]
dimnames(cord1)[1]<-dimnames(x$DataMatrix)[2]
}#end catype

} #end bip row
else{
plottype<-"biplot"
#biptype<-"columns"
if (ell==TRUE) { scaleplot<-1
}
if ((x$catype=="CA")|(x$catype=="NSCA")){
cord1<- x$Cprinccoord*scaleplot
cord2<-x$Rstdcoord/scaleplot
nthings<-x$cols
nvars<-x$rows
thinggroup<-colgroup
thinggrlab<-colgrlab
vargroup<-colgroup
vargrlab<-rowgrlab
thinglabels<-x$collabels
varlabels<-x$rowlabels
main<-"Column Isometric Biplot"
inertiapc=x$inertias[,2] #inertia of row 
dmu=diag(x$inertias[,1])
dimnames(cord1)[1]<-dimnames(x$DataMatrix)[2]
dimnames(cord2)[1]<-dimnames(x$DataMatrix)[1]
}
if ((x$catype=="DONSCA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA"))
{
if (ell==TRUE) {
scaleplot<-1}
cord2<- x$Cprinccoord*scaleplot
cord1<-x$Rstdcoord/scaleplot
nthings<-x$rows
nvars<-x$cols
thinggroup<-rowgroup
thinggrlab<-rowgrlab
vargroup<-colgroup
vargrlab<-colgrlab
thinglabels<-x$rowlabels
varlabels<-x$collabels
inertiapc=x$inertias[,2] #inertia of row poly
dmu=diag(x$inertias[,1])
dimnames(cord1)[1]<-dimnames(x$DataMatrix)[1]
dimnames(cord2)[1]<-dimnames(x$DataMatrix)[2]
}#end catype
}#end bip column
}


###################################################################################ok without choice plottype
#repeat{

if ((x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA"))
{
 cat("\n Looking at the Trends of rows and columns\n")
#browser()
#trendplot(x@mj,(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
#windows()
#trendplot(x@mi,t(x@pcc), position="topleft",main="Original Trends of Centered #Column Profile",xlab="ordered scores")
############################################################## reconstructed TREND
#windows()
trendplot(x$mj,(x$Trend), posleg=posleg,main="Reconstructed rows of the centred column profile",xlab="ordered scores",prop=prop)
plot.new()
trendplot(x$mi,t(x$Trend), posleg=posleg,main="Reconstructed columns of the centred column profile",xlab="ordered scores",prop=prop)
#browser()
}

##############################################
picsize1<-c(range(cord1[,c(firstaxis,lastaxis)], cord2[,c(firstaxis,lastaxis)])/prop)

if (picsize1[1]>=picsize1[2]) stop(paste("incorrect axis scale picsize =", picsize1[1], picsize1[2], "\n\n"))
########################################################################################## 
if ((x$catype=="DONSCA")||(x$catype=="DOCA"))
{
plotone (firstaxis,lastaxis,plottype=plottype,things=x$catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,type="b",
catype=x$catype,pos=pos) 
}
######################
if ((x$catype=="SOCA")||(x$catype=="SONSCA"))
{
if (biptype=="row"){type="b"}
else {type="p"}
plotone (firstaxis,lastaxis,plottype=plottype,things=x$catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,catype=x$catype,type=type,pos=pos) 
}

###############################
if ((x$catype=="CA")||(x$catype=="NSCA"))
{
plotone (firstaxis,lastaxis,plottype=plottype,things=x$catype,nthings,nvars,cord1,cord2,
inertiapc=inertiapc, thinggroup,thinggrlab,vargroup,vargrlab, thinglabels, varlabels,picsize=picsize1,cex=cex,type="p",catype=x$catype,pos=pos) 
}
################################################################################
#cat("\nIncluding Beh's Confidence Ellipses\n")
################################################################################
if (ell==TRUE) {

if (((x$catype=="NSCA")|(x$catype=="CA")|(x$catype=="DOCA")|(x$catype=="SOCA")|(x$catype=="SONSCA")|(x$catype=="DONSCA")) & (plottype=="biplot")&(biptype=="row")|(biptype=="r")|(biptype=="rows")){
cordr<-cord2
cordc<-cord1
cord1<-cordr
cord2<-cordc
}

plot.new()

switch(x$catype, "CA"=caellipse(N=x$DataMatrix,a1=firstaxis,a2=lastaxis,alpha=alpha,M=M,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,b=x$Cstdcoord,g=cord1,fr=cord2,dmu=dmu,inertiapc=inertiapc,
plottype=plottype,biptype=biptype,pos=pos,arrow=T,length=0,graphy=T), 

"SOCA"=caellipse(N=x$DataMatrix,a1=firstaxis,a2=lastaxis,alpha=alpha,M=M,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=solve(x$Imass^0.5)%*%x$Rstdcoord,b=solve(x$Jmass^0.5)%*%x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,inertiapc=inertiapc,
plottype=plottype,biptype=biptype,pos=pos,arrow=T,length=0,graphy=T),

"DOCA"=caellipse(N=x$DataMatrix,a1=firstaxis,a2=lastaxis,alpha=alpha,M=M,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=solve(x$Imass^0.5)%*%x$Rstdcoord,
b=solve(x$Jmass^0.5)%*%x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,inertiapc=inertiapc,
plottype=plottype,biptype=biptype,pos=pos,arrow=T,length=0,graphy=T), 

"NSCA"=nscaellipse(x$DataMatrix,a1=firstaxis,a2=lastaxis,alpha=alpha,M=M,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord1,fr=cord2,dmu=dmu, tauden=x$tauden,
inertiapc=inertiapc,
plottype=plottype,biptype=biptype,pos=pos,arrow=T,length=0,graphy=T),

"SONSCA"=nscaellipse(x$DataMatrix,a1=firstaxis,a2=lastaxis,alpha=alpha,M= M,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,tauden=x$tauden,inertiapc=inertiapc,
plottype=plottype,biptype=biptype,pos=pos,arrow=T,length=0,graphy=T), 

"DONSCA"=nscaellipse(x$DataMatrix,a1=firstaxis,a2=lastaxis,alpha=alpha,M=M,cex=cex,cex.lab=cex.lab,prop=prop,
Imass=x$Imass,Jmass=x$Jmass,a=x$Rstdcoord,
b=x$Cstdcoord,g=cord2,fr=cord1,dmu=dmu,tauden=x$tauden,inertiapc=inertiapc,
plottype=plottype,biptype=biptype,pos=pos,arrow=T,length=0,graphy=T))

}#end if ellipse

}
