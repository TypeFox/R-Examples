MAT <-
function(ipar,
              resp,
              cors,
              target.content.dist=NULL,
              content.cat=NULL,
              ncc=1,
              content.order=NULL,
              p=stop("p is required"),
              selectionMethod=c("D","A","C","R"),
              selectionType=c("FISHER","BAYESIAN"),
              c.weights=NA,
              stoppingCriterion=c("CONJUNCTIVE","COMPENSATORY"),
              topN=1,
              minNI=10,
              maxNI=30,
              minSE=0.3,
              D=1.0,
              maxIter=30,
              conv=0.001,
              minTheta=-4,
              maxTheta=4,
              plot.audit.trail=T,
              theta.labels=NULL,
              easiness=T) {

call<-match.call()

if (is.data.frame(ipar)) ipar<-cbind(as.matrix(ipar[paste("a",1:p,sep="")]),as.matrix(ipar["d"]),as.matrix(ipar["c"]))

if (!easiness) ipar[,p+1]<- -ipar[,p+1]

if (is.data.frame(resp)) resp<-as.matrix(resp)

selectionMethod<-toupper(selectionMethod)
selectionMethod<-match.arg(selectionMethod)

stoppingCriterion<-toupper(stoppingCriterion)
stoppingCriterion<-match.arg(stoppingCriterion)

selectionType<-toupper(selectionType)
selectionType<-match.arg(selectionType)

if (p==1) sigma<-1 
      else if (p>1) {
sigma<-as.matrix(cors)
sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
if (dim(sigma)[1]!=dim(sigma)[2] || p!= dim(sigma)[1]) stop("ERROR: p and cors non-conforming")
} 
      else stop("ERROR: p and cors non-conforming")

if (selectionMethod=="C") {
if (sum(c.weights)!=1 || length(c.weights) !=p) {
selectionMethod="D"
c.weights<-rep(1/p,p)
warning("ERROR: c.weights do not sum to 1.0 or have unexpected values - D optimality will be used")
}
}

ni<-nrow(ipar)
nj<-nrow(resp)

available.items<-rep(TRUE,ni)

content.balancing<-F

if (ncc>1 && !is.null(target.content.dist) & !is.null(content.cat)){

if (abs(sum(target.content.dist)-1)>.1) warning("ERROR: target content proportions does not add up to 1.0\n:content balancing will not be used")
else if (length(target.content.dist)!=ncc) warning("ERROR: ncc does not match the number of target proportions\n:content balancing will not be used")
else if (length(content.cat)!=ni) warning("ERROR: number of rows in the content control file does not match the number of items in the bank\n:content balancing will not be used")
else {
overall.content.freq<-numeric(ncc)
content.balancing<-T

if (any(target.content.dist<=0)) {
available.items[which(content.cat %in% which(target.content.dist<=0))]<-FALSE
}
if (any(target.content.dist<0)) {
  warning("ERROR: values in target.cont.dist cannot be negative")
}
}
}
###add warning messages for content ordering
if (!is.null(content.order)) {
  if (anyDuplicated(content.order)) { content.order <- NULL
    warning("ERROR: content categories are duplicated in content.order\n:content ordering will not be used")
  } else if (any(content.order<1) | any(content.order>p)) { content.order <- NULL
    warning("ERROR: content categories are out of range in content.order\n: must specify content categories as integers
            from 1 to ncc, where ncc is the number of content categories\n:content ordering will not be used")
  } else if (length(content.order)!= ncc) { content.order <- NULL
    warning("ERROR: number of content categories specified in content.order is not equal to ncc, where ncc is the number of
            content categories\n:content ordering will not be used")
  } else if (!is.null(content.order) & (minNI!=maxNI)) { content.order <- NULL
    warning("ERROR: minNI must equal maxNI (i.e., fixed-length test) when content.order is specified\n:content ordering will not be used.") 
  } else if (is.null(target.content.dist)) { content.order <- NULL
    warning("ERROR: content distributions must be specified to use content ordering\n:content ordering will not be used.") 
  }
}

next.content<-function() {
  if (!is.null(content.order)) { 
    completed.content <- which((target.content.dist-current.content.dist)<=0)
    available.content <- 1:ncc
    available.content[completed.content] <- NA
    available.content <- available.content[content.order]
    available.content <- available.content[!is.na(available.content)]
    idx <- available.content[1]
    return(idx)
  }
  else {
    available.content<-which(target.content.dist>0)
    idx<-which.max(target.content.dist[available.content]-current.content.dist[available.content])
    return(available.content[idx])
  }
}

update.content.dist<-function() {
  idx<-content.cat[item.selected]
  current.content.freq[idx]<<-current.content.freq[idx]+1
  overall.content.freq[idx]<<-overall.content.freq[idx]+1
  if (!is.null(content.order)) {
    content.dist.denom <- c(minNI,ni.given)
    current.content.dist<<-current.content.freq/(content.dist.denom[which.max(content.dist.denom)])
  } else {
    current.content.dist<<-current.content.freq/ni.given
  }
}

estimates.full<-SCORE_cpp(ipar,resp,p,sigma,maxIter=maxIter,conv=conv,D=D,Fisher=T)

estimates.full$theta <- round(estimates.full$theta,4)
estimates.full$SE <- round(estimates.full$SE,4)

TH<-matrix(nrow=nj,ncol=p)
SE<-matrix(nrow=nj,ncol=p)

start.theta<-rep(0,p)

items.used<-matrix(NA,nj,maxNI)
selected.item.resp<-matrix(NA,nj,maxNI)
ni.administered<-numeric(nj)
theta.CAT<-matrix(NA,nj,p)
se.CAT<-matrix(NA,nj,p)
theta.history<-array(NA,c(nj,maxNI,p))
se.history<-array(NA,c(nj,maxNI,p))

if (plot.audit.trail) dev.new(record=T,width=10,height=6.5)
if (is.null(theta.labels)) theta.labels<-paste("Theta",1:p,sep=" ")

for (j in 1:nj) {
  se.met<-logical(p)
  theta.current<-start.theta
  rs<-resp[j,]
  
  W<-matrix(0,nrow=p,ncol=p) #initialize
  
  crit.met<-FALSE
  items.available<-available.items
  
  for (i in 1:ni) {
    if (is.na(rs[i]) || (rs[i] !=0 && rs[i] !=1)) items.available[i]<-FALSE
  }
  items.given<-rep(FALSE,ni)
  
  max.to.administer<-ifelse(sum(items.available)<=maxNI,sum(items.available),maxNI)
  ni.given<-0
  
  if (content.balancing) {
    current.content.dist<-numeric(ncc)
    current.content.freq<-numeric(ncc)
  }
  
  
  while (crit.met==FALSE && ni.given<max.to.administer) {
    
    available <- items.available
    
    if (content.balancing) {
      content.available<-available & (content.cat==next.content())
      
      if (sum(content.available)>=1) available<-content.available
    }
    
    info <- selectItem_cpp(ipar,available,items.given,theta.current,p,sigma,D=D,method=selectionMethod,selectionType=selectionType,c_weights=c.weights,content_balancing=content.balancing,topN=topN)
    
    info[!available] <- NA
    
    info.index<-rev(order(info))
    if (topN > 1 && sum(items.available) > topN) { item.selected<-info.index[sample(topN,1)]
    } else { item.selected=which(info==max(info,na.rm =T))[1] }
    
    ni.given<-ni.given+1
    items.used[j,ni.given]<-item.selected
    
    
    if (content.balancing) { update.content.dist() }
    
    items.available[item.selected]<-FALSE
    items.given[item.selected]<-TRUE
    
    selected.item.resp[j,ni.given]<-resp[j,item.selected]
    
    estimates<-score_cpp(ipar[items.used[j,1:ni.given],,drop=F],rs[items.used[j,1:ni.given]],theta.current,p,sigma,maxIter=maxIter,conv=conv,D=D,Fisher=T)
    
    theta.current<-estimates$theta
    theta.history[j,ni.given,]<-theta.current
    se.history[j,ni.given,]<-estimates$SE
    
    if (stoppingCriterion=="CONJUNCTIVE") { se.met<-all(estimates$SE<=minSE)
    } else if (stoppingCriterion=="COMPENSATORY") {
      V<-solve(-estimates$Hessian)
      if (selectionMethod=="C") cV<-sqrt(abs(matrix(c.weights,nrow=1)%*%V%*%matrix(c.weights,ncol=1)))
      else cV<-sqrt(abs(matrix(rep(1/p,p),nrow=1)%*%V%*%matrix(rep(1/p,p),ncol=1)))
      
      se.met<-cV<=minSE
    }
    
    if (ni.given>=max.to.administer || (se.met && ni.given>=minNI)) {
      crit.met<-TRUE
      theta.CAT[j,]<-round(estimates$theta,4)
      se.CAT[j,]<-round(estimates$SE,4)
      ni.administered[j]<-ni.given
    }
  }
  
  
  if (plot.audit.trail) {
    plot(1:maxNI,seq(minTheta,maxTheta,length=maxNI),main=paste("CAT Audit Trail - Examinee ",j,sep=""),xlab="Items Administered",ylab="Theta",type="n",las=1,bg="white")
    idx<-0
    for (h in p:1) {
      idx<-idx+1
      points(1:ni.given,theta.history[j,1:ni.given,h],type="b",pch=8+h,lty=h,col=h)
      abline(h=estimates.full$theta[j,h],lty=h,col=h)
      text(1,minTheta+0.3*(idx),paste(theta.labels[h]," : ",sprintf("%6.3f",theta.CAT[j,h])," SE: ",sprintf("%5.3f",se.CAT[j,h])),cex=0.8,adj=0);
    }
    if (stoppingCriterion=="COMPENSATORY") text(1,minTheta+0.3*(idx+1),paste("cSE=",round(cV,digits=3),sep=""),cex=0.8,adj=0)
    if (p>1) legend("bottomright",theta.labels,col=1:p,lty=1:p,pch=8+1:p,bg="white")
    else if (p==1) {
      for (i in 1:ni.given) {
        lines(rep(i,2),c(theta.history[j,i,1]-1.96*se.history[j,i,1],theta.history[j,i,1]+1.96*se.history[j,i,1]),col="blue")
      }
    }
    item.string<-paste(items.used[j,1:ni.given],collapse=",")
    text(1,maxTheta,paste("Items: ",item.string,sep=""),cex=0.7,adj=0)
    
    resp.string<-paste(selected.item.resp[j,1:ni.given],collapse=",")
    text(1,maxTheta-(maxTheta-minTheta)/20,paste("Responses: ",resp.string,sep=""),cex=0.7,adj=0)
    
    if (content.balancing) text(1,maxTheta-(maxTheta-minTheta)/10,paste("Content Distribution:",paste(round(current.content.dist,digits=2),collapse=",")),cex=0.7,adj=0)
    
  }
}

if (content.balancing) {
  
overall.content.dist<-overall.content.freq/sum(overall.content.freq)
content.dist<-rbind(target.content.dist,overall.content.dist)
par.mar<-par()$mar
par(xpd=T,mar=par()$mar+c(0,0,0,4))
barplot(content.dist,ylab="proportion",xlab="content category",las=1,names.arg=paste(1:ncc),col=c("black","grey"),beside=T)
legend(ncc*3+0.2,max(target.content.dist,current.content.dist)/2,c("Target","Current"),fill=c("black","grey"))
par(xpd=F,mar=par.mar)
}

dev.new(record=T,width=10,height=6.5)
par(mfrow=c(2,3))
for (h in 1:p) {
plot(estimates.full$theta[,h],theta.CAT[,h],xlim=c(minTheta,maxTheta),ylim=c(minTheta,maxTheta),xlab=paste("Full Bank",theta.labels[h]),ylab=paste("CAT Theta",theta.labels[h]),col="blue")
text(minTheta,maxTheta,paste("r =",sprintf("%5.3f",cor(estimates.full$theta[,h],theta.CAT[,h]))),adj=0)
abline(0,1)
}


out<-list(call=call,items.used=items.used,selected.item.resp=selected.item.resp,ni.administered=ni.administered,
                theta.CAT=theta.CAT,se.CAT=se.CAT,theta.history=theta.history,se.history=se.history,theta.Full=estimates.full$theta,se.Full=estimates.full$SE,ipar=ipar,p=p)

class(out)<-"MAT"

return(out)
}

