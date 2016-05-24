####################################################################################
### Inequality index

calcSGini<-function(x, w=NULL, param=2)
{
  calcConc(x, y=NULL, w, param)
}


calcConc<-function(x, y, w, param)
{
  if (param<=0) return(NULL)
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  if (!is.null(y))
  {
    if(!is.numeric(y)) return(NULL)
    yNA<-sum(as.numeric(is.na(y)))
  }
  weighted<-FALSE
  wNA<-NULL
  if (is.null(w)) w<-rep(1, length(x))
  else 
  {
    if (!is.numeric(w)) return(NULL)
    weighted<-TRUE
    wNA<-sum(as.numeric(is.na(w)))
  }
  if (is.null(y)) df<-cbind("x"=x, "w"=w)
  else df<-cbind("x"=x, "y"=y, "w"=w)
  df<-df[complete.cases(df),, drop=FALSE]
  if (nrow(df)==0) return (NULL)
  if (any(df[,"x"]<0)) return(NULL)
  if (sum(df[,"x"])==0) return(NULL)
  index<-0
  names(param)<-"param"
  if (nrow(df)>1)
  {
    if (param != 1)
    {
      df[,"w"]<-df[,"w"]/sum(df[,"w"])
      xMean<-weighted.mean(df[,"x"],df[,"w"])
      if (is.null(y)) df<-df[order(df[,"x"]),]
      else df<-df[order(df[,"y"]),]
      sp<-cumsum(df[,"w"])
      sp[length(sp)]<-1
      sm<-c(0, sp[-length(sp)])
      wr<-(((1-sm)^param)-((1-sp)^param))/(param*df[,"w"])
      wvar<-sum(df[,"w"] * ((wr - weighted.mean(wr, df[,"w"]))^2))
      reg1<-coef(lm((-param*wvar*df[,"x"]/xMean)~wr, weights=df[,"w"]))
      names(reg1)<-NULL
      index<-reg1[2]
    }
  }
  if (is.null(y)) names(index)<-"SGini"
  else names(index)<-"SConc"
  SG<-list(ineq=   list(index=index,
                        parameter=param),
           nas=    NULL)
  if (is.null(y)) SG[["nas"]]<-list(xNA=xNA, wNA=wNA,
                                    totalNA=length(x)-nrow(df))
  else SG[["nas"]]<-list(xNA=xNA, yNA=yNA, wNA=wNA,
                          totalNA=length(x)-nrow(df))
  class(SG)<-"ICI"
  return(SG)
}


calcGEI<-function(x, w=NULL, alpha=1)
{
  xNA<-sum(as.numeric(is.na(x)))
  weighted<-FALSE
  wNA<-NULL
  if (is.null(w)) w<-rep(1, length(x))
  else 
  {
    if (!is.numeric(w)) return(NULL)
    weighted<-TRUE
    wNA<-sum(as.numeric(is.na(w)))
  }
  df<-cbind("x"=x,"w"=w)
  df<-df[complete.cases(df),, drop=FALSE]
  if(nrow(df)==0) return (NULL)
  if(any(df[,"x"]<0)) return(NULL)
  if(sum(df[,"x"])==0) return(NULL)
  if(any(df[,"x"]==0) && alpha==0) return(NULL)
  if(any(df[,"x"]==0) && alpha==1) return(NULL)
  index<-0
  names(index)<-"GEI"
  names(alpha)<-"alpha"
  GEI<-list(ineq=   list(index=index,
                        parameter=alpha),
            nas=    list(xNA=xNA, wNA=wNA, totalNA=length(x)-nrow(df)))
  class(GEI)<-"ICI"
  if(nrow(df)==1) return(GEI)  
  if (alpha==0) index<-calcTheil0(df)
  if (alpha==1) index<-calcTheil1(df)
  if (alpha!=1 && alpha!=0)
  {
    xMean<-weighted.mean(df[,"x"], df[,"w"])
    index<-weighted.mean(((df[,"x"]/xMean)^alpha)-1, df[,"w"])/(alpha*(alpha-1))
  }
  names(index)<-"GEI"
  GEI[["ineq"]][["index"]]<-index
  return(GEI)  
}


calcTheil1<-function(df)
{
  x<-df[,1]
  w<-df[,2]
  w<-w/sum(w)
  xMean<-weighted.mean(x,w)
  x<-x/xMean
  return(sum(w*x*log(x)))
}


calcTheil0<-function(df)
{
  x<-df[,1]
  w<-df[,2]
  w<-w/sum(w)
  xMean<-weighted.mean(x,w)
  x<-x/xMean
  return(-sum(w*log(x)))
}


calcAtkinson<-function(x, w=NULL, epsilon=1)
{
  if (epsilon<0) return(NULL)
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  weighted<-FALSE
  wNA<-NULL
  if (is.null(w)) w<-rep(1, length(x))
  else 
  {
    if (!is.numeric(w)) return(NULL)
    weighted<-TRUE
    wNA<-sum(as.numeric(is.na(w)))
  }
  df<-cbind("x"=x,"w"=w)
  df<-df[complete.cases(df),, drop=FALSE]
  if(nrow(df)==0) return (NULL)
  if(any(df[,"x"]<0) || sum(df[,"x"])==0) return(NULL)
  if (any(df[,"x"]==0) && epsilon==1) return(NULL)
  index<-0
  names(index)<-"Atk"
  names(epsilon)<-"epsilon"
  Atk<-list(ineq=   list(index=index,
                        parameter=epsilon),
            nas=    list(xNA=xNA, wNA=wNA, totalNA=length(x)-nrow(df)))
  class(Atk)<-"ICI"
  if(nrow(df)==1) return(Atk)
  if (epsilon==1)
  {
    w<-df[,"w"]/sum(df[,"w"])
    xMean<-weighted.mean(df[,"x"],w)
    index<-1-(prod(exp(w*log(df[,"x"])))/xMean)
  }
  else
  {
    xMean<-weighted.mean(df[,"x"], df[,"w"])
    x1<-df[,"x"]/xMean
    w<-df[,"w"]/sum(df[,"w"])
    param<-1-epsilon
    index<-1-(weighted.mean((x1^param), w)^(1/param))
  }
  names(index)<-"Atk"
  Atk[["ineq"]][["index"]]<-index
  return(Atk)
}



####################################################################################
### Concentration index


calcSConc<-function(x, y, w=NULL, param=2)
{
  calcConc(x, y, w, param)
}



####################################################################################
### Curves
#   ordinary Lorenz curve, generalized LC, concentration curve 

curveLorenz<-function(x, w=NULL, gener=FALSE,
                      xlab=NA, ylab=NA, add=FALSE, grid=0, ...)
{
  curveLoCo(x, y=NULL, w, gener, xlab, ylab, add, grid, ...)
}


curveConcent<-function(x, y, w=NULL,
                      xlab=NA, ylab=NA, add=FALSE, grid=0, ...)
{
  curveLoCo(x, y, w, gener=FALSE, xlab, ylab, add, grid, ...)
}


curveLoCo<-function(x, y, w, gener, xlab, ylab, add, grid, ...)
{
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  if (!is.null(y))
  {
    if (!is.numeric(y)) return(NULL)
    yNA<-sum(as.numeric(is.na(y)))
  }
  weighted<-FALSE
  wNA<-NULL
  if (is.null(w)) w<-rep(1, length(x))
  else 
  {
    if (!is.numeric(w)) return(NULL)
    weighted<-TRUE
    wNA<-sum(as.numeric(is.na(w)))
  }
  if (is.null(y)) df<-cbind("x"=x,"w"=w)
  else df<-cbind("x"=x, "y"=y, "w"=w)
  df<-df[complete.cases(df),]
  if(any(df[,"x"]<0)) return(NULL)
  if(sum(df[,"x"])==0) return(NULL)
  if(nrow(df)==1) return(NULL)
  df[,"w"]<-df[,"w"]/sum(df[,"w"])
  if (is.null(y)) df<-df[order(df[,"x"]),]
  else df<-df[order(df[,"y"]),]
  if (gener) ymaxi<-weighted.mean(df[,"x"], df[,"w"])
  else ymaxi<-1
  if (!add)
  {
    plot(x=c(0,0), y=c(1,1), type = "l", lty = 1, las=1,
         xlim=c(0,1), ylim=c(0,ymaxi),
         xaxs = "i", yaxs = "i", cex.axis=0.75, cex.lab=1, xlab=xlab, ylab=ylab)
    if (grid>0)
    {
      abline(v=c((0:grid)/grid)[-c(1, grid+1)], col = "gray80")
      abline(h=c((0:grid*ymaxi)/grid)[-c(1, grid+1)], col = "gray80")

    }
    abline(0,ymaxi)
  }
  lines(cumsum(df[,"w"]),
        cumsum(ymaxi*(df[,"w"]*df[,"x"])/sum(df[,"w"]*df[,"x"])),
        type = "l", ...)
}


 
####################################################################################
### Decomposition of index
#

decompSGini<-function(x, z, w=NULL, param=2, decomp="BM", ELMO=TRUE)
{
  if (param<=0) return(NULL)
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  if (!is.factor(z)) return(NULL)
  zNA<-sum(as.numeric(is.na(z)))
  weighted<-FALSE
  wNA<-NULL
  if (is.null(w)) w<-rep(1, length(x))
  else 
  {
    if (!is.numeric(w)) return(NULL)
    weighted<-TRUE
    wNA<-sum(as.numeric(is.na(w)))
  }
  df<-data.frame("x"=x, "z"=z, "w"=w)
  df<-df[complete.cases(df),, drop=FALSE]
  if (nrow(df)==0) return (NULL)
  if (any(df[,"x"]<0)) return(NULL)
  if (sum(df[,"x"])==0) return(NULL)
  if (nrow(df)==1) return(NULL)
  names(param)<-"param"
  lx<-length(x)
  df[, "z"]<-factor(df[,"z"], exclude=NULL)
  df[, "w"]<-df[, "w"]/sum(df[, "w"])
  df<-df[order(df[,"x"]),]
  dfSplit<-split(df[,c("x","w")], df[,"z"])
  xMean<-weighted.mean(df[,"x"],df[,"w"])
  xMeanIntra<-sapply(dfSplit,
                  function(df) weighted.mean(df[,"x"],df[,"w"]), simplify=TRUE)
  wIntra<-sapply(dfSplit, function(df) sum(df[,"w"]), simplify=TRUE)
  sIntra<-wIntra*xMeanIntra/xMean
  if (weighted)
  {
    SGiniIntra<-sapply(dfSplit,
          function(df) calcSGini(df[,"x"],df[,"w"], param)[["ineq"]][["index"]],
          simplify=FALSE)
    SGini<-calcSGini(df[,"x"],df[,"w"], param)[["ineq"]][["index"]]
  }
  else
  {
    SGiniIntra<-sapply(dfSplit,
            function(df) calcSGini(df[,"x"],NULL, param)[["ineq"]][["index"]],
            simplify=FALSE)
    SGini<-calcSGini(df[,"x"],NULL, param)[["ineq"]][["index"]]
  }
  SGiniIntra<-unlist(SGiniIntra)
  names(SGiniIntra)<-names(wIntra)
  names(SGini)<-"SGini"
  if (ELMO)
  {
    nc<-outer(cumsum(df[, "w"]), cumsum(wIntra[order(xMeanIntra)]),
              FUN=function(x,y) {x-y>.Machine$double.eps^0.5})
    nc<-factor(apply(nc, 1, sum))
    dfSplitE<-split(df[,c("x","w")], nc)
    xMeanIntraE<-sapply(dfSplitE, function(df) weighted.mean(df[,"x"],df[,"w"]),
                        simplify=TRUE)
    wIntraE<-sapply(dfSplitE, function(df) sum(df[,"w"]), simplify=TRUE)
    elmo<-calcSGini(xMeanIntraE, wIntraE, param)[["ineq"]][["index"]]
    names(elmo)<-NULL
  }
  else elmo<-NULL
  if (decomp=="BM")
  {
    SGiniIntraContrib<-wIntra*sIntra*SGiniIntra
    SGiniIntraTot<-sum(SGiniIntraContrib)
    SGiniInter<- calcSGini(xMeanIntra, wIntra, param)[["ineq"]][["index"]]
    SGiniResid<- SGini - SGiniIntraTot - SGiniInter
    names(SGiniInter)<-NULL
    names(SGiniResid)<-NULL
    SGD<-list(ineq=   list(index=SGini,
                          parameter=param),
              decomp= list(within=SGiniIntraTot,
                          between=SGiniInter,
                          overlap=SGiniResid,
                          betweenELMO=elmo),
              intra=  list(SGiniGroups=SGiniIntra,
                          contribSGiniGroups=SGiniIntraContrib),
              stratif=NULL,              
              ws=     list(wIntra=wIntra,
                          sIntra=sIntra),
              nas=    list(xNA=xNA, zNA=zNA, wNA=wNA, totalNA=lx-nrow(df)))
  }
  if (decomp=="YL")
  {
    SGiniIntraContrib<-sIntra*SGiniIntra
    SGiniIntraTot<-sum(SGiniIntraContrib)
    sp<-cumsum(df[,"w"])
    sp[length(sp)]<-1
    sm<-c(0, sp[-length(sp)])
    wr<-(((1-sm)^param)-((1-sp)^param))/(param*df[,"w"])
    wrSplit<-split(wr, df[,"z"])
    wrMeanIntra<-rep(0, length(wrSplit)) 
    SGiniStratif<-rep(0, length(wrSplit))
    for (i in 1: length(wrSplit))
    {
      dfS<-dfSplit[[i]]
      dfS[,"w"]<-dfS[,"w"]/sum(dfS[,"w"])
      wrMeanIntra[i]<-weighted.mean(wrSplit[[i]], dfS[,"w"])
      if (nrow(dfS) >1)
      {
        spS<-cumsum(dfS[,"w"])
        spS[length(spS)]<-1
        smS<-c(0, spS[-length(spS)])
        wrS<-(((1-smS)^param)-((1-spS)^param))/(param*dfS[,"w"])
        b1<-coef(lm(wrS~dfS[,"x"], weights=dfS[,2]))[2]
        b2<-coef(lm(wrSplit[[i]]~dfS[,"x"], weights=dfS[,"w"]))[2]
        SGiniStratif[i]<-(1-(b2/b1))/(1-wIntra[i])
      }
    }  
    wr<-wrMeanIntra
    wrmean<-weighted.mean(wr, wIntra)
    wrvar<-sum(wIntra * ((wr - wrmean)^2))
    SGiniInter<-coef(lm((-param*wrvar*xMeanIntra/xMean)~wr, weights=wIntra))[2]
    names(SGiniInter)<-NULL
    names(SGiniStratif)<-names(wIntra)
    SGiniStratifContrib<-sIntra*SGiniIntra*SGiniStratif*(wIntra-1)
    SGiniStratifTot<-sum(SGiniStratifContrib)
    names(param)<-"param"
    SGD<-list(ineq=   list(index=SGini,
                          parameter=param),
              decomp= list(within=SGiniIntraTot,
                          between=SGiniInter,
                          stratif=SGiniStratifTot,
                          betweenELMO=elmo),
              intra=  list(SGiniGroups=SGiniIntra,
                          contribSGiniGroups=SGiniIntraContrib),
              stratif=list(stratifGroups=SGiniStratif,
                            contribStratifGroups=SGiniStratifContrib),              
              ws=     list(wIntra=wIntra,
                          sIntra=sIntra),
              nas=    list(xNA=xNA, zNA=zNA, wNA=wNA, totalNA=lx-nrow(df)))
  }
class(SGD)<-"ICI"
return(SGD)  
}


decompGEI<-function(x, z, w=NULL, alpha=1, ELMO=TRUE)
{
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  if (!is.factor(z)) return(NULL)
  zNA<-sum(as.numeric(is.na(z)))
  weighted<-FALSE
  wNA<-NULL
  if (is.null(w)) w<-rep(1, length(x))
  else 
  {
    if (!is.numeric(w)) return(NULL)
    weighted<-TRUE
    wNA<-sum(as.numeric(is.na(w)))
  }
  df<-data.frame("x"=x, "z"=z, "w"=w)
  df<-df[complete.cases(df),, drop=FALSE]
  if (nrow(df)==0) return (NULL)
  if (any(df[,"x"]<0)) return(NULL)
  if (sum(df[,"x"])==0) return(NULL)
  if(any(df[,"x"]==0) && alpha==0) return(NULL)
  if(any(df[,"x"]==0) && alpha==1) return(NULL)
  if (nrow(df)==1) return(NULL)
  names(alpha)<-"alpha"
  lx<-length(x)
  df[, "z"]<-factor(df[,"z"], exclude=NULL)
  df[, "w"]<-df[, "w"]/sum(df[, "w"])
  dfSplit<-split(df[,c("x","w")], df[,"z"])
  xMean<-weighted.mean(df[,"x"],df[,"w"])
  xMeanIntra<-sapply(dfSplit, function(df) weighted.mean(df[,"x"],df[,"w"]),
                      simplify=TRUE)
  wIntra<-sapply(dfSplit, function(df) sum(df[,"w"]), simplify=TRUE)
  sIntra<-wIntra*xMeanIntra/xMean
  GEIw<-sapply(dfSplit,
          function(df) calcGEI(df[,"x"], df[,"w"], alpha)[["ineq"]][["index"]],
          simplify=FALSE)
  GEIw<-unlist(GEIw)
  names(GEIw)<-NULL
  GEIwContrib<-GEIw*(sIntra^alpha)*(wIntra^(1-alpha))
  GEIwTot<-sum(GEIwContrib)
  GEIb<-calcGEI(xMeanIntra, wIntra, alpha)[["ineq"]][["index"]]
  index<-GEIwTot+GEIb
  names(index)<-"GEI"
  names(GEIb)<-NULL
  names(GEIw)<-names(wIntra)
  if (ELMO)
  {
    df<-df[order(df[,"x"]),]
    nc<-outer(cumsum(df[, "w"]), cumsum(wIntra[order(xMeanIntra)]),
              FUN=function(x,y) {x-y>.Machine$double.eps^0.5})
    nc<-factor(apply(nc, 1, sum))
    dfSplit<-split(df[,c("x","w")], nc)
    xMeanIntraE<-sapply(dfSplit, function(df) weighted.mean(df[,"x"],df[,"w"]),
                        simplify=TRUE)
    wIntraE<-sapply(dfSplit, function(df) sum(df[,"w"]), simplify=TRUE)
    elmo<-calcGEI(xMeanIntraE, wIntraE, alpha)[["ineq"]][["index"]]  
    names(elmo)<-NULL
  }
  else elmo<-NULL
  GEID<-list(ineq=   list(index=index,
                        parameter=alpha),
            decomp= list(within=GEIwTot,
                        between=GEIb,
                        betweenELMO=elmo),
            intra=  list(GEIGroups=GEIw,
                        contribGEIGroups=GEIwContrib),
            ws=     list(wIntra=wIntra,
                        sIntra=sIntra),
            nas=    list(xNA=xNA, zNA=zNA, wNA=wNA, totalNA=lx-nrow(df)))
  class(GEID)<-"ICI"
  return(GEID)
}


decompAtkinson<-function(x, z, w=NULL, epsilon=1, decomp="BDA", ELMO=TRUE)
{
  if (epsilon<0) return(NULL)
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  if (!is.factor(z)) return(NULL)
  zNA<-sum(as.numeric(is.na(z)))
  weighted<-FALSE
  wNA<-NULL
  if (is.null(w)) w<-rep(1, length(x))
  else 
  {
    if (!is.numeric(w)) return(NULL)
    weighted<-TRUE
    wNA<-sum(as.numeric(is.na(w)))
  }
  df<-data.frame("x"=x, "z"=z, "w"=w)
  df<-df[complete.cases(df),, drop=FALSE]
  if (nrow(df)==0) return (NULL)
  if (any(df[,"x"]<0)) return(NULL)
  if (any(df[,"x"]==0) && epsilon==1) return(NULL)
  if (sum(df[,"x"])==0) return(NULL)
  if (nrow(df)==1) return(NULL)
  names(epsilon)<-"epsilon"
  lx<-length(x)
  df[, "z"]<-factor(df[,"z"], exclude=NULL)
  df[, "w"]<-df[, "w"]/sum(df[, "w"])
  dfSplit<-split(df[,c("x","w")], df[,"z"])
  xMean<-weighted.mean(df[,"x"],df[,"w"])
  xMeanIntra<-sapply(dfSplit, function(df) weighted.mean(df[,"x"],df[,"w"]), simplify=TRUE)
  wIntra<-sapply(dfSplit, function(df) sum(df[,"w"]), simplify=TRUE)
  sIntra<-wIntra*xMeanIntra/xMean
  index<-calcAtkinson(df[,"x"], df[,"w"], epsilon)[["ineq"]][["index"]]
  names(index)<-"Atk"
  Atkw<-sapply(dfSplit,
    function(df) calcAtkinson(df[,"x"], df[,"w"], epsilon)[["ineq"]][["index"]],
    simplify=FALSE)
  Atkw<-unlist(Atkw)
  names(Atkw)<-names(wIntra)
  Atkb<-calcAtkinson(xMeanIntra, wIntra, epsilon)[["ineq"]][["index"]]
  names(Atkb)<-NULL
  if (ELMO)
  {
    df<-df[order(df[,"x"]),]
    nc<-outer(cumsum(df[, "w"]), cumsum(wIntra[order(xMeanIntra)]),
              FUN=function(x,y) {x-y>.Machine$double.eps^0.5})
    nc<-factor(apply(nc, 1, sum))
    dfSplit<-split(df[,c("x","w")], nc)
    xMeanIntraE<-sapply(dfSplit, function(df) weighted.mean(df[,"x"],df[,"w"]),
                        simplify=TRUE)
    wIntraE<-sapply(dfSplit, function(df) sum(df[,"w"]), simplify=TRUE)
    elmo<-calcAtkinson(xMeanIntraE, wIntraE, epsilon)[["ineq"]][["index"]]  
    names(elmo)<-NULL
  }
  else elmo<-NULL
  if (decomp=="BDA")
  {
    if (epsilon==1) AtkwTot<-1-prod(exp(wIntra*log(1-Atkw)))
    else
    {
      wse<-(wIntra^epsilon)*(sIntra^(1-epsilon))
      AtkwTot<-1-((sum(wse*((1-Atkw)^(1-epsilon)))/sum(wse))^(1/(1-epsilon)))
      names(AtkwTot)<-NULL
    }
    Atkc<-AtkwTot*Atkb
    AtkD<-list(ineq=   list(index=index,
                          parameter=epsilon),
              decomp= list(within=AtkwTot,
                          between=Atkb,
                          cross=Atkc,
                          betweenELMO=elmo),
              intra=  list(AtkGroups=Atkw),
              ws=     list(wIntra=wIntra,
                          sIntra=sIntra),
              nas=    list(xNA=xNA, zNA=zNA, wNA=wNA, totalNA=lx-nrow(df)))
  }
  if (decomp=="DP")
  {
    if (epsilon==1)
    {
      AtkwTot<-prod(exp(wIntra[Atkw!=0]*log((xMeanIntra[Atkw!=0]/xMean)*Atkw[Atkw!=0])))
    }
    else
    {
      wse<-(wIntra^epsilon)*(sIntra^(1-epsilon))
      AtkwTot<-(sum(wse[Atkw!=0]*(Atkw[Atkw!=0]^(1-epsilon))))^(1/(1-epsilon))
    }
    names(AtkwTot)<-NULL
    Atkr<-index-AtkwTot-Atkb
    names(Atkr)<-NULL
    AtkD<-list(ineq=   list(index=index,
                          parameter=epsilon),
              decomp= list(within=AtkwTot,
                          between=Atkb,
                          residual=Atkr,
                          betweenELMO=elmo),
              intra=  list(AtkGroups=Atkw),
              ws=     list(wIntra=wIntra,
                          sIntra=sIntra),
              nas=    list(xNA=xNA, zNA=zNA, wNA=wNA, totalNA=lx-nrow(df)))
  }
  class(AtkD)<-"ICI"
  return(AtkD)
}


summary.ICI<-function(object, ...)
{
  x<-object
  if (is.null(x)) return(NULL)
  mc<-match.call()
  digits<-ifelse("digits" %in% names(mc),
                  mc[["digits"]], max(5, getOption("digits") - 5))
  names(x[["ineq"]])<-NULL
  xi<-unlist(x[["ineq"]])
  print.default(xi, print.gap=8, quote=FALSE, digits=digits)
  if ("decomp" %in% names(x))
  {
    typeDecomp<-""
    if (names(xi)[1]=="SGini" && length(x[["decomp"]])>2)
    {
      typeDecomp<- " (BM)"
      if (names(x[["decomp"]])[3]=="stratif") typeDecomp<- " (YL)"
    }
    if (names(xi)[1]=="Atk" && length(x[["decomp"]])>2)
    {
      typeDecomp<- " (DP)"
      if (names(x[["decomp"]])[3]=="cross") typeDecomp<- " (BDA)"
    }
    cat(paste("\nDecomposition", typeDecomp, ":\n", sep=""))
    decom<- unlist(x[["decomp"]])
    if (names(decom)[length(decom)]!="betweenELMO") print.default(decom,
                                        print.gap=2, quote=FALSE, digits=digits)
    else
    {
      print.default(decom[-length(decom)], print.gap=2,
                    quote=FALSE, digits=digits)
      print.default(decom[length(decom)], print.gap=2,
                    quote=FALSE, digits=digits)
    }
  }
}