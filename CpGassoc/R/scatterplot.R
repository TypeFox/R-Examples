scatterplot <-
function(x,cpg.rank=NULL,cpg.name=NULL,save.plot=NULL,file.type="pdf",eps.size=c(5,5),popup.pdf=FALSE,beta.values=NULL,main.title=NULL,...)   {

  nam.ind<-x$info$Phenotype
  if(!((x$info$betainfo) %in% ls(.GlobalEnv)) & is.null(beta.values)){
     stop("\nbeta values is no longer in the environment. Either provide the new name of the beta values\n object (see help(scatterplot))",
          " or reload the object\n")
      }
  if(is.null(cpg.rank) & is.null(cpg.name)) {
    stop("No cpg sites given. Please specify what sites you want to plot\n")
      }
  if(is.null(beta.values)) {y=eval(parse(text=as.character(x$info$betainfo)))}
  else {y<-beta.values}
  if(!is.null(cpg.rank) & !is.null(cpg.name)) {
    cpg.name=NULL
    warning("Used cpg.rank\n")
    }
  if(!is.null(cpg.rank)) {
    x$results<-x$results[order(x$results$P.value),]
    place<-as.character(x$results[cpg.rank,1])
    if(!is.null(save.plot)) {
      save.plot<-paste(save.plot,cpg.rank,sep="_")
         }}
  if(!is.null(cpg.name)) {
    place<-cpg.name
    if(!is.null(save.plot)) {
      save.plot<-paste(save.plot,cpg.name,sep="_")
    }}


  holder<-which(x$results[,1] %in% place)
  if(!is.null(x$covariates)) {
    covar<-paste(names(x$covariates),collapse=", ")
    }
  if(!is.null(x$chip) & !x$info$random) {
    if(!is.null(x$covariates)) {
      covar<-paste(covar,x$info$chipinfo,sep=" and ")}
    else{
      covar<-x$info$chipinfo
    } }    
 if(x$info$logittran & !is.factor(x$indep)) {
  sampvalues<-seq(min(x$indep),max(x$indep),length.out=200)
  }

 betasval<-t(y[place,])
 h.r<-(length(place)==1 & is.matrix(y)) | is.null(dimnames(betasval))
 if(h.r) {
	betasval<-t(betasval)
  dimnames(betasval)<-list(NULL,place)
		}

   
 if(is.null(dimnames(betasval)) & length(place)!=1) {
    betasval<-t(betasval)
    dimnames(betasval)<-list(1:length(betasval),place)
    }

 
 for(i in 1:length(place)){
        if(!is.null(save.plot)){
          if(!(file.type %in% c("eps","pdf"))) {
              stop("Incorrect file type. Must be pdf or eps\n")
                } 
          if(file.type=="eps"){
              postscript(paste(save.plot[i],".eps",sep=""), horizontal = FALSE, 
                onefile = FALSE, paper = "special",width=eps.size[1],height=eps.size[2])
                }
          if(file.type=="pdf" & !popup.pdf) {
                pdf(file=paste(save.plot[i],".pdf",sep=""))
              }
                  }
        titleinfo<-paste(place[i]," modeled on ",nam.ind,sep="")
        if(!is.null(x$covariates) |(!x$info$random & !is.null(x$chip))) {
          if(!is.null(x$covariates) & x$info$Num.Cov >4) {
              titleinfo<-paste(titleinfo,"in the presense of\n",x$info$Num.Cov,"covariates")
              }
          else{
            titleinfo<-paste(titleinfo," in the presence of:\n",covar,sep="")
          }}
        if(!is.null(main.title)) {
                titleinfo<-paste(place[i],": ",main.title,sep="")
        
                }
       if(is.factor(x$indep)) {
         
         boxplot(betasval[,place[i]]~x$indep,xlab=nam.ind,ylab=place[i],main=titleinfo,
            sub=paste("P.value= ",format(x$results$P.value[holder[i]],digits=3)),...)
          }
        else{
          plot(betasval[,place[i]]~x$indep,xlab=nam.ind,ylab=place[i],main=titleinfo,
              sub=paste("P-value = ",format(x$results$P.value[holder[i]],digits=3)),...) 
          if(!x$info$logittran) {
          abline(a=x$coefficients[place[i],3],b=x$coefficients[place[i],4],col="red")
              }
          else {
             if(is.factor(place)) {place<-as.character(place)}
             explogitbeta<-x$coefficients[place[i],3]+x$coefficients[place[i],4]*sampvalues
             testbeta<-exp(explogitbeta)/(1+exp(explogitbeta))
          points(sampvalues,testbeta,type="l",col="red")
            }
         }
    if(!is.null(save.plot))  {
       if(file.type=="eps") {
           dev.off()
            }
      else{
        if(popup.pdf) {
          dev.copy2pdf(file=paste(save.plot[i],".pdf",sep=""))
            }
        else{
            dev.off()
            }
      
      }}
  cat("Press enter to continue\n")
  readline()
  if(i==length(place)) {
        cat("All ",length(place)," sites plotted\n",sep="")
          break}  
     }}
