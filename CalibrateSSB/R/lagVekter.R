

lagVekter = function(
  modellen,
  brutto,
  pop = NULL,
  min_w = -Inf,
  max_w = Inf,
  totalReturn=0,
  maxiter=10,
  samplingWeights=NULL,
  returnwGross=FALSE)
{
  z = brutto

  lm_modell = lm(as.formula(modellen),data=z)


  if(is.numeric(pop) & is.vector(pop))
  {
    if(length(pop)==1)
    {
       totPop=pop
    }
    else
    {
      totPop=setTotal(pop,lmObject=lm_modell)
    }
  }
  else
  {
    if(is.null(pop))
    {
      tp = getTotal(z,lm_modell,samplingWeights)
    }
    else
    {
      tp = getTotal(pop,lm_modell)
    }

    if(totalReturn>1) return(tp)
    totPop =  tp$N*tp$colSum/tp$colN
  }
  if(totalReturn>0) return(totPop)



  x_model = model.matrix(lm_modell)


  r_model = model.frame(lm_modell)[1][,]



  a = ImpVekt(x_model,r_model ,totPop=totPop,totPopReturn=TRUE,w=samplingWeights)

  w0=a$vekt



  w = w0[r_model==1]
  x_netto = x_model[r_model==1,]

  minw_eps = 1E-6 * max(min_w,1)
  maxwPlusEps = (1+1E-6) * max_w   # Avoid Inf-Inf

  fortsett = (sum(w<(min_w+minw_eps)|w>maxwPlusEps )>0)
  iter=0
  if(fortsett & !is.null(samplingWeights)) stop("Limits&SamplingWeights not implemented")
  while(fortsett)
  {
    iter=iter+1
    wFixed = NA+w
    wFixed[wFixed>min_w] = NA
    wFixed[w<min_w] = min_w
    wFixed[w>max_w] = max_w
    w = ImpVektFixed(x_netto,totPop=a$totPop,wFixed=wFixed)

    fortsett = (sum(w<(min_w-minw_eps)|w>maxwPlusEps )>0)

    if(fortsett)
      if((iter>=maxiter)|(sum(w[!is.na(wFixed)]<(min_w-minw_eps)|w[!is.na(wFixed)]>maxwPlusEps )>0))
      {
        w = w0[r_model==1]
        fortsett = FALSE
      }
  }

  w1 = w0
  w1[r_model==1] = w
  if(returnwGross){
    return(list(w=w1,wGross=a$wGross))

  }
  w1
}


setTotal = function(total,lmObject)
{
  x   = model.matrix(lmObject)[1,]
  x[1:length(x)] = NA
  varnames = names(total)[names(total) %in% names(x)]
  for(i in 1:length(varnames))
    x[names(x)==varnames[i]] = total[names(total)==varnames[i]]
  x
}


getTotal = function(data,lmObject,w=NULL)
{
  x=model.frame(lmObject)
  x1=dim(x)[1]
  x2=dim(x)[2]
  d1=dim(data)[1]
  while(dim(x)[1]<x1+d1)
  {
      diffn = x1+d1-dim(x)[1]
      if(diffn>dim(x)[1]) x=rbind(x,x)
      else x=rbind(x,x[1:diffn,])
  }
  x[(x1+1):(x1+d1),]=NA
  varnames =names(data)[names(data) %in% names(x)]
  if(length(varnames)>0) for(i in 1:length(varnames))
    x[(x1+1):(x1+d1),names(x)==varnames[i]] = data[,names(data)==varnames[i],drop=FALSE]
  m=model.matrix(lmObject,data=x,na.action=NULL)
  m=m[(x1+1):(x1+d1), ,drop=FALSE]
  x=NULL
  if(is.null(w))
  {
    x$colSum =  colSums(m,na.rm=TRUE)
    x$colN   =  colSums(!is.na(m))
    x$N = dim(m)[1]
  } else
  {
    x$colSum =  colSums(w*m,na.rm=TRUE)
    x$colN   =  colSums(w*(!is.na(m)))
    x$N = sum(w)
  }
  x
}


sumList = function(x,xNew)
{
  for(i in 1:length(x)) x[[i]] =  x[[i]]+ xNew[[i]]
  x
}






