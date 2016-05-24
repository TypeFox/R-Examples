plotBVS = function(results, num.models=100, num.snps=20, num.regions=20, plot.coef=FALSE, true.coef=NULL,main=NULL, regions=NULL, type="s", prop.cases=NULL,...) {
  
  which <- data.frame(results$Which[,1:(dim(results$Which)[2]-2)])
  #for(i in 1:dim(which)[2]){which[,i] = as.numeric(which[,i])}
  which.r <- results$Which.r
  u.regions = colnames(which.r)
  snps <- colnames(results$Which)[-c((dim(which)[2]+1):(dim(which)[2]+3))]
  p = length(snps)
  num.snps = min(p,num.snps)
  num.models = min((dim(which)[1]-1),num.models)
  postprob <- as.numeric(results$Which[,dim(results$Which)[2]])
  bf <- as.numeric(results$MargBF)
  bf.r <- as.numeric(results$Marg.RBF)
  if(plot.coef==TRUE){
  	coef = results$Coef
  }

  null.ind = c(1:dim(which)[1])[apply(which,1,paste,collapse="")==paste(rep(0,p),collapse="")]
  null.post <- postprob[null.ind]
  which = as.matrix(which[-null.ind,])
  if(length(which.r)>0){
  	  which.r = as.matrix(which.r[-null.ind,])}  
  postprob = postprob[-null.ind]
  model.ord <- order(-postprob)
  which <- as.matrix(which[model.ord, ])
  if(length(which.r)>0){
      which.r <- as.matrix(which.r[model.ord, ])}  
  postprob <- postprob[model.ord]
  if(plot.coef==TRUE){
  	coef = coef[-null.ind]
  	coef = coef[model.ord]
  }	
  

  ## Graphic Parameters
  if(type=="s"){
  nmodel <- num.models
    nvar <- num.snps
    ordr <- order(-bf)
    snps = snps[ordr][1:nvar]
    regions = regions[ordr][1:nvar]
    rownms <- paste(snps,regions,sep="\n")
    clr <- c("#FFFFFF", "#A020F0", "#0000CD")
    ordr <- ordr[1:nvar]
    color.matrix = which[1:(nmodel),ordr[1:nvar]]
    for(i in 1:dim(color.matrix)[2]){color.matrix[,i] = as.numeric(color.matrix[,i])}
    color.matrix <- color.matrix + 2
    if(length(true.coef)>0){
        prob.labels <- paste("Marg BF:",as.numeric(round(bf[ordr[1:nvar]], 2)),
                             " \nTrue OR:",round(true.coef[ordr[1:nvar]],2),sep="")}
    if(length(prop.cases)>0){
    	prob.labels <- paste("Marg BF:",as.numeric(round(bf[ordr[1:nvar]], 2)),
                             " \nCases: ",prop.cases[ordr[1:nvar],1]," Controls: ",prop.cases[ordr[1:nvar],2],sep="")}
    if(length(true.coef)==0 & length(prop.cases)==0){
    	prob.labels <- paste("Marg BF:",as.numeric(round(bf[ordr[1:nvar]], 2)))}
    

    
  ## matrix of colors white, purple(la), blue(dom), red(rec)    
  keep.mar <- par(mar=c(5, 6, 4, 2) + 0.1)
  par(las=1, mar=c(8, 10, 5, 10), ps=10, font=2)
  maintitle = main
  if(length(main)==0){
     maintitle=paste("SNP Inclusions of Top Models \nGlobal BF=",round(results$Global,1))}
  prob.axis <- postprob[1:(nmodel)]
  prob.axis <- prob.axis/sum(prob.axis) 
  image(c(0, cumsum(prob.axis)), 1:nvar, as.matrix(color.matrix), col=clr, 
        xlab="Models Ranked by Post. Prob.", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1),
        main=maintitle)

  xat <- (cumsum(prob.axis) + c(0, cumsum(prob.axis[-nmodel]))) / 2
  
  if(plot.coef==FALSE){
      axis(1, at=xat, labels=c(1:num.models))}
  if(plot.coef==TRUE){
  	  beta.labels <- round(coef[1:(nmodel)],1)
  	  if(nmodel>5){
         beta.labels[5:nmodel] = NA}
  	  axis(1, at=xat,labels=beta.labels)
  }    
  axis(2, at=1:nvar, labels=rownms)
  axis(4, at=1:nvar, labels=prob.labels)
  par(mar=keep.mar)}
  
  if(type=="r"){
    nmodel <- num.models
    nvar <- min(length(u.regions),num.regions)
    ordr <- order(-bf.r)
    u.regions = u.regions[ordr]
    rownms = u.regions[1:nvar]
    clr <- c("#FFFFFF", "#A020F0", "#0000CD")
    ordr <- ordr[1:nvar]
    which.r[which.r>1] = 1
    color.matrix = which.r[1:(nmodel),ordr[1:nvar]]
    for(i in 1:dim(color.matrix)[2]){color.matrix[,i] = as.numeric(color.matrix[,i])}
    color.matrix <- color.matrix + 2
    if(length(true.coef)>0){
        prob.labels <- paste("Region BF:",as.numeric(round(bf.r[ordr[1:nvar]], 2)),
                             " \nTrue OR:",round(true.coef[ordr[1:nvar]],2),sep="")}
    if(length(true.coef)==0){
    	prob.labels <- paste("Region BF:",as.numeric(round(bf.r[ordr[1:nvar]], 2)))}
    

    
  ## matrix of colors white, purple(la), blue(dom), red(rec)    
  keep.mar <- par(mar=c(5, 6, 4, 2) + 0.1)
  par(las=1, mar=c(8, 12, 5, 12), ps=10, font=2)
  maintitle = main
  if(length(main)==0){
     maintitle=paste("Region Inclusions of Top Models \nGlobal BF=",round(results$Global,1))}
  prob.axis <- postprob[1:(nmodel)]
  prob.axis <- prob.axis/sum(prob.axis) 
  image(c(0, cumsum(prob.axis)), 1:nvar, as.matrix(color.matrix), col=clr, 
        xlab="Models Ranked by Post. Prob.", ylab="", xaxt="n", yaxt="n", xlim=c(0, 1),
        main=maintitle)

  xat <- (cumsum(prob.axis) + c(0, cumsum(prob.axis[-nmodel]))) / 2
  if(plot.coef==FALSE){
      axis(1, at=xat, labels=c(1:num.models))}
  if(plot.coef==TRUE){
  	  beta.labels <- round(coef[1:(nmodel)],1)
  	  if(nmodel>5){
         beta.labels[5:nmodel] = NA}
  	  axis(1, at=xat,labels=beta.labels)
  }  
  axis(2, at=1:nvar, labels=rownms)
  axis(4, at=1:nvar, labels=prob.labels)
  par(mar=keep.mar)}
}


	

	

