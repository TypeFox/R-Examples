##  This file includes some ad hoc R functions to import Luminex data,
##  preprocess Luminex array data and some others.  In ADA, we make
##  most of these functions as internal functions and use them only in
##  preparing the data for analysis.
print.luminex <- function(x, ...){
  print(summary(as.data.frame(x[c("Treat0","Treat2","Treat1",
                                  "Treat3","Background")])),
        ...)
  invisible(x)
}

.Import.Luminex <- function(sid, path="../data/"){
  rtfile1 = paste(path,'rtcarda2.csv',sep='')
  rtfile2 = paste(path,'rtcardb2.csv',sep='')
  lufile = paste(path,'luminex.csv',sep='')
  rtdata1 = read.csv(rtfile1, header=TRUE)
  rtdata2 = read.csv(rtfile2, header=TRUE)
  rtdata3 = rbind(rtdata1[,c(1,2,6)], rtdata2[,c(2,3,6)])
  rtdata3 = data.frame(rtdata3, Card = c(rep("A",nrow(rtdata1)),rep("B",nrow(rtdata2))))
  ludata1 = read.csv(lufile, header=TRUE)
  ludata = data.frame(Gene=ludata1$Sample, Pool=ludata1$Pool,
	xc = ludata1[,sid*5-2], xbg = ludata1[,sid*5+2],
	xdox = ludata1[,sid*5-1], xcis = ludata1[,sid*5],
	xifo = ludata1[,sid*5+1])
  sele = as.numeric(substr(rtdata3$Sample,1,2)) == sid
  rtdata3 = rtdata3[sele,]
  rtdata = data.frame(Gene = substr(as.character(rtdata3$Detector),1,
                        nchar(as.character(rtdata3$Detector))-8),
	Treat = substr(as.character(rtdata3$Sample),4,
          nchar(as.character(rtdata3$Sample))),
	RQ = rtdata3$RQ, Card=rtdata3$Card)
  treats = as.character(rtdata$Treat)
  sele = treats == "Ifos";
  treats[sele] = "Ifo"
  sele = !is.na(rtdata$RQ)
  rtdata = data.frame(Gene=rtdata$Gene, Treat=treats,RQ=rtdata$RQ, Card=rtdata$Card)[sele,]
  list(lu=ludata,rt=rtdata) 
}

.Intra.Mean <- function(x){
  treats = levels(x$Treat);treats
  ref=treats[1]
  pools = as.numeric(levels(as.factor(x$Pool)));
  ntreat = length(treats);ntreat
  npool = length(pools);npool
  if(is.na(match(ref,treats)))
    stop("Reference treatment is not found.")
  X = x$Signal - x$Bg
  sele1 = x$Pool == 1;
  sele2 = x$Treat == ref;
  xref = X[sele1&sele2];xref
  ypool = NULL;ytreat=NULL;ycf = NULL;
  for(j in 1:ntreat){
    sele2 = x$Treat == treats[j];
    for(i in 1:npool){
      sele1 = x$Pool == pools[i];
      ratios = mean(xref/X[sele1&sele2]);
      ypool = c(ypool,pools[i]);
      ytreat = c(ytreat,treats[j])
      ycf = c(ycf,ratios)
    }
  }
  list(Thetas=data.frame(Pool=ypool,Treat=ytreat,CorFactor=1/ycf))
}

.Intra.Median <- function(x){
  treats = levels(x$Treat);treats
  ref=treats[1]
  pools = as.numeric(levels(as.factor(x$Pool)));
  ntreat = length(treats);ntreat
  npool = length(pools);npool
  if(is.na(match(ref,treats)))
    stop("Reference treatment is not found.")
  X = x$Signal - x$Bg
  sele1 = x$Pool == 1;
  sele2 = x$Treat == ref;
  xref = X[sele1&sele2];xref
  ypool = NULL;ytreat=NULL;ycf = NULL;
  for(j in 1:ntreat){
    sele2 = x$Treat == treats[j];
    for(i in 1:npool){
      sele1 = x$Pool == pools[i];
      ratios = median(xref/X[sele1&sele2]);
      ypool = c(ypool,pools[i]);
      ytreat = c(ytreat,treats[j])
      ycf = c(ycf,ratios)
    }
  }
  list(Thetas=data.frame(Pool=ypool,Treat=ytreat,CorFactor=1/ycf))
}

.Intra.Me <- function(x,alpha=0.05,cutoff=3,iter=100,tol=1.e-6){
  treats = levels(x$Treat);treats
  ref=treats[1]
  pools = as.numeric(levels(as.factor(x$Pool)));
  ntreat = length(treats);ntreat
  npool = length(pools);npool
  if(is.na(match(ref,treats)))
    stop("Reference treatment is not found.")
  l=sum(x$Signal<x$Bg); n =nrow(x);l;n;
  if(l<30){l=30}else{l=2*l}; l
  se = sd(x$Signal[x$Signal<quantile(x$Signal,l/n)]);se
  x$Signal = x$Signal-x$Bg;  ## subtract the backgrounds
  Thetas = matrix(0,nrow=ntreat,ncol=npool); #initialize 
  cpool = 1; ctreat = ref;ctreat# select reference pool
  selec = x$Treat==ctreat & x$Pool==cpool; sum(selec)
  nbeads = levels(as.factor(as.character(x$Gene))); nbeads
  ref = x$Signal[selec];ref;
  i=0;
  for(treat in treats){
    i=i+1;j=0;
    sele1 = x$Treat==treat
    for(pool in pools){
      j=j+1
      sele2 = x$Pool==pool;
      ref2 = (x$Signal)[sele1&sele2];ref2;
      selep = ref>0&ref2>0;
      if(length(ref)!=length(ref2))
        warning("Normalization beads not match!");
      nfactor = (prod(ref[selep]/ref2[selep]))^.25;nfactor
      x$Signal[sele1&sele2] = x$Signal[sele1&sele2]*nfactor;
      Thetas[i,j]=nfactor;
    }
  }
  ## Screen normalization beads by strength and dispersion
  sigs = NULL; nblist = NULL;  nb0=NULL; mu=NULL;
  for(nbead in nbeads){
    sele = x$Gene==nbead;
    xb = x[sele,]
    if(any(xb$Signal<cutoff*se)){
      warning(paste(nbead, ": weak signal(s)!"))
    }else{
      nblist = c(nblist,nbead);
      sigs = c(sigs,sd(log(xb$Signal)))
      mu = c(mu,mean(log(xb$Signal)))
    }
  }
  m= nrow(xb);m;sigs
  Sig = sqrt(mean(sigs^2));Sig
  sigs2 = NULL; nblist2 = NULL;mu2=NULL
  for(i in 1:length(sigs)){
    chisq = (m-1)*(sigs[i]/Sig)^2;
    alpha = min(alpha/2,.5-alpha/2);
    if(chisq<=qchisq(alpha,m-1)|chisq>=qchisq(1-alpha,m-1)){
      warning(paste(nblist[i], "large dispersion!"))
    }else{
      nblist2 = c(nblist2,nblist[i]);
      sigs2 = c(sigs2,sigs[i])
      mu2 = c(mu2,mu[i])
    }
  }
  if(length(nblist2)<1)stop("No stable normlization beads can be used!")
  Sig = sqrt(mean(sigs2^2));Sig

  Iter=1; dllk=1;rllk=1; # stopping criteria
  llk0 = 0; llk1 =0;
  while((dllk>tol|rllk>tol)&Iter<iter){
    Iter=Iter+1
    i=0;
    sigs=NULL;mu=NULL;
    for(k in 1:length(nblist2)){
      nbead = nblist2[k]
      sele3 = x$Gene==nbead;
      xb = x$Signal[sele3]
      sigs = c(sigs,sd(xb))
      mu = c(mu,mean(xb))
    }
    llk1=0.
    for(treat in treats){
      i=i+1;j=0;
      sele1 = x$Treat==treat
      for(pool in pools){
        j=j+1;
        sele2 = x$Pool==pool;
        nfactor = exp(mean(log(x$Signal[sele1&sele2]))
                -mean(log(x$Signal)))
        x$Signal[sele1&sele2] = x$Signal[sele1&sele2]*nfactor
        Thetas[i,j] = Thetas[i,j]*nfactor;
        for(k in 1:length(nblist2)){
          nbead = nblist2[k]
          sele3 = x$Gene==nbead;
          sele = sele1&sele2&sele3;
          llk1 = llk1 - sigs[k] - 0.5*(log(x$Signal[sele])
            -log(Thetas[i,j])-log(mu[k]))^2/(sigs[k])^2;
        }
      }
    }
    dllk=abs(llk1-llk0);
    rllk=abs(dllk/llk1);
    llk0 = llk1
  }
  
  ypool = NULL;ytreat=NULL;ycf = NULL;
  for(j in 1:ntreat){
    sele2 = x$Treat == treats[j];
    for(i in 1:npool){
      sele1 = x$Pool == pools[i];
      ypool = c(ypool,pools[i]);
      ytreat = c(ytreat,treats[j])
      ycf = c(ycf,Thetas[j,i])
    }
  }
  out = data.frame(Pool=ypool,Treat=ytreat,CorFactor=1/ycf)

  list(Thetas=out,Iter=Iter-1,Tol=max(dllk,rllk))
}

.Intra.Norm <- function(x,method='me',
                       alpha=0.05,cutoff=3,iter=100,tol=1.e-6){
  if(!is.data.frame(x)) stop("'x' must be a data frame.")
  if(ncol(x)!=5) stop("'x' shall have five columns.")
  if(any(!names(x) == c("Gene","Pool","Treat","Signal","Bg")))
    names(x) = c("Gene",  "Pool",  "Treat","Signal","Bg")
  method <- match.arg(tolower(method),
                      c("me","mean","median"))
  out <- switch(method,
                mean = .Intra.Mean(x),
                median = .Intra.Median(x),
                me = .Intra.Me(x,alpha,cutoff,iter,tol)
                )
  list(Thetas=out$Thetas,Iter=out$Iter);
}

.FCLoess <- function(y,x){
  M=log(y/x);A=(log(y)+log(x))/2.;
  loessm = loess(M~A);
  y/x*exp(loessm$res);
}
.FCLoessM <- function(y,x){
  M=log(y/x);A=(log(y)+log(x))/2.;
  loessm = loess(M~A);
  y/x*exp(loessm$res-median(M));
}
.FCQuantile <- function(y,x){
  ody = order(y);odx=order(x);
  avg = (sort(y)+sort(x))/2;
  avg[ody]/avg[odx];
}


.normalize.Luminex <- function(x,intra='me',inter='none',treat='Cis',cutoff = 2.0){
  ## get the profiles of the normalizers
  Treat0 = treat
  sele = substr(x$Gene,1,4) =="Norm"
  x2 = x[sele,];
  ndata = rbind(
    data.frame(Gene=x2$Gene,Pool=x2$Pool,Treat='Ctrl',Signal=x2$xc,Bg=x2$xbg),
    data.frame(Gene=x2$Gene,Pool=x2$Pool,Treat='Dox',Signal=x2$xdox,Bg=x2$xbg),
    data.frame(Gene=x2$Gene,Pool=x2$Pool,Treat='Cis',Signal=x2$xcis,Bg=x2$xbg),
    data.frame(Gene=x2$Gene,Pool=x2$Pool,Treat='Ifo',Signal=x2$xifo,Bg=x2$xbg)
    )
  method <- match.arg(tolower(intra),
                      c("me","mean","median"))
  y = .Intra.Norm(ndata,method=method)$Thetas 
  ##  filter the weak signals
  sbg = sd(x$xbg);sbg
  treats = names(x)
  treat = paste('x',tolower(Treat0),sep='')
  i = match(treat,treats);
  if(is.na(i)) stop("The treatment is not found.")
  if(i<5) stop("Invalid treatment.")
  sele1 = x$xc >= x$xbg+cutoff*sbg
  sele2 = x[,i] >= x$xbg+cutoff*sbg
  sele3 = substr(x$Gene,1,3) =="hsa"
  sele = sele1 & sele2 & sele3
  if(sum(sele)<1) stop("No record selected.")
  x = x[sele,]

  npools = as.numeric(levels(as.factor(x2$Pool)))
  x$xc = x$xc - x$xbg
  x[,i] = x[,i] - x$xbg
  
  for(j in npools){
    sele1x = x$Pool == j
    sele1y = y$Pool == j
    sele2y = tolower(as.character(y$Treat)) == tolower(Treat0)
    sele3y = y$Treat == "Ctrl"
    sele = sele1y&sele2y
    if(sum(sele)!=1) stop(paste("Correction factor not found:",j, Treat0))
    nfact = y$CorFactor[sele]
    x[sele1x,i] = x[sele1x,i]/nfact
    sele = sele1y&sele3y
    if(sum(sele)!=1) stop(paste("Correction factor not found:",j, "Ctrl"))
    nfact = y$CorFactor[sele]
    x$xc[sele1x] = x$xc[sele1x]/nfact
  }
  x3 = data.frame(Gene=x$Gene, Ctrl=x$xc, Treat=x[,i])

  y2 = x3$Treat; x2 = x3$Ctrl
  nmethod <- match.arg(tolower(inter),
                       c('none',"loess","loessm","quantile"))
  fc = switch(tolower(nmethod),
    loess = .FCLoess(y2,x2),
    loessm = .FCLoessM(y2,x2),
    quantile = .FCQuantile(y2,x2),
    none = y2/x2
    );
  data.frame(x3, FC=fc)
}
