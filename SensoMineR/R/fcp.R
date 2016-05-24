fcp <- function(X,group,scale=TRUE, ncp = NULL, axes=c(1,2), name.group = NULL, level.conf = 0.95, nbsim=500, nbchoix=NULL, cex=1, color=NULL, title=NULL, new.plot=TRUE, graph=c("ind","var","ellipse")){

if (scale) type = rep("s",length(group))
if (!scale) type = rep("c",length(group))
res <- MFA(X,group=group,type=type,ncp=Inf,name.group=name.group,graph=FALSE)
if (new.plot) dev.new()
if ("ind"%in%graph) plot(res,cex=1,title=title,color=color,new.plot=FALSE)
if ("var"%in%graph){
  if ("ind"%in%graph) dev.new()
  plot(res,choix="var",cex=1,title=title,color=color,new.plot=FALSE,habillage="group")
}
result <- list()
result$mfa <- res
if ("ellipse"%in%graph){
  if (("ind"%in%graph)|("var"%in%graph)) dev.new()
  result$ellipse <- boot(X, method = "freechoice", axes = axes, scale = scale, ncp = ncp, group = group, nbsim = nbsim, level.conf = level.conf,
    nbchoix = nbchoix,color = color,cex = cex, title = title, new.plot = FALSE)
}
return(result) 
}
