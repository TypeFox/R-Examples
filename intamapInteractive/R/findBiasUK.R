findBiasUK = function (object)
  {

  locations = coordinates(object$observations)
  value = object$observations$value 
  borders = object$observations$regCode
  regCode = object$observations$regCode

  categ = findCateg(object)  
  conti = findConti(object)  
  fixed = findFixed(object)  
  
  if ((length(conti) > 0) & (length(fixed)==0))
    {
    mainDF=data.frame(
                      x=locations[,1],y=locations[,2],
                      area=value,
                      regCode=regCode,
                      conti=conti,
                      categ=categ,
                      borders=borders,
                      value=value,
                      value.weight=value,
                      res=value)
    }
  if ((length(conti) == 0) & (length(fixed)>0))
    {
    mainDF=data.frame(
                      x=locations[,1],y=locations[,2],
                      area=value,
                      regCode=regCode,
                      fixed=fixed,
                      categ=categ,
                      borders=borders,
                      value=value,
                      value.weight=value,
                      res=value)
    }

  if ((length(categ) > 0) & (length(conti) == 0) & (length(fixed)==0))
    {
    mainDF=data.frame(
                      x=locations[,1],y=locations[,2],
                      area=value,
                      regCode=regCode,
                      categ=categ,
                      borders=borders,
                      value=value,
                      value.weight=value,
                      res=value)
    }
    
  if ((length(categ) == 0) & (length(conti) == 0) & (length(fixed)==0))
    {
    mainDF=data.frame(
                      x=locations[,1],y=locations[,2],
                      area=value,
                      regCode=regCode,
                      borders=borders,
                      value=value,
                      value.weight=value,
                      res=value)
    }
  mainDF=na.omit(mainDF)
  mainDF=mainDF[duplicated(mainDF[,1:2])==F,]

  x_s=mainDF$x              
  y_s=mainDF$y
                             
  coordinates(mainDF)=~x+y
  nl=length(mainDF$value)
  V=diag(1,nl,nl)
  Vmes=diag(0,nl,nl)

# OLS
  if ((length(conti) > 0) & (length(fixed) == 0))
    {
    t1=BLUE(  data=mainDF, 
              fixed2=fixed, conti2="conti", categ2="categ", bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  if ((length(categ) > 0) & (length(conti) == 0) & (length(fixed) == 0))
    {
    t1=BLUE(  data=mainDF, 
              fixed2=fixed, conti2=conti, categ2="categ", bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  if ((length(categ) == 0) & (length(conti) == 0) & (length(fixed) == 0))
    {
    t1=BLUE(  data=mainDF, 
              fixed2=fixed, conti2=conti, categ2=categ, bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  if ((length(conti) == 0) & (length(fixed) > 0))
    {
    t1=BLUE(  data=mainDF, 
              fixed2="fixed", conti2=conti, categ2="categ", bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  mainDF$res=mainDF$value-c(t1$drifts%*%na.omit(t1$bias$coef))
  vResFit = autofitVariogram(res ~ 1, mainDF)
  Dx=Dy=D=matrix(nrow=nl,ncol=nl)
  x_s=coordinates(mainDF)[,1]
  y_s=coordinates(mainDF)[,2]
  sqDiff=function(x,y) (x-y)^2
  for (i in seq(1,nl))
    	{
     	Dx[i,]=sqDiff(x=x_s[i],y=x_s)
     	Dy[i,]=sqDiff(x=y_s[i],y=y_s)
     	}
  D=sqrt(Dx+Dy)
  for (i in 1:dim(D)[2]) 
    {
    V[,i]=variogramLine(vResFit$var_model,dist_vector=D[,i], debug.level = 0, covariance=T)[,2]
    }

# First GLS
  if ((length(conti) > 0) & (length(fixed) == 0))
    {
    t2=BLUE(  data=mainDF, 
              fixed2=fixed, conti2="conti", categ2="categ", bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  if ((length(categ) > 0) & (length(conti) == 0) & (length(fixed) == 0))
    {
    t2=BLUE(  data=mainDF, 
              fixed2=fixed, conti2=conti, categ2="categ", bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  if ((length(categ) == 0) & (length(conti) == 0) & (length(fixed) == 0))
    {
    t2=BLUE(  data=mainDF, 
              fixed2=fixed, conti2=conti, categ2=categ, bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  if ((length(conti) == 0) & (length(fixed) > 0))
    {
    t2=BLUE(  data=mainDF, 
              fixed2="fixed", conti2=conti, categ2="categ", bias.factor="regCode", 
              n=nl, V=V, Vmes=Vmes)
    }
  if ((length(categ) == 0) & (length(conti) == 0) & (length(fixed) == 0))
    {
    t2$bias$coef=t2$bias$coef-mean(t2$bias$coef)
    }
  if (length(categ)>0) return(list(regionalBias = data.frame(regCode=unique(regCode), cbias=t2$bias[,1])))
  else return(list(regionalBias = data.frame(regCode=unique(regCode), cbias=t2$bias[,1]-mean(t2$bias[,1]))))
  }


##################### Internal functions ########################33

findCateg = function (object)  
  {
  return(NULL)
  }

findConti = function (object)  
  {
  return(NULL)
  }

findFixed = function (object)  
  {
  return(NULL)
  }

factorDisj = function (obj, colnames)
# This function transforms a factorial drift column into a disjunctive matrix
# INPUT
#
#     obj     data frame containing the factors and drifts
#     colname name of the column to use
	{
  return(sapply(unique(obj[[colnames[1]]]),function(x){ ifelse(obj[[colnames[1]]]== x,1,0)}))
  }

BLUE = function(data, fixed2, conti2, categ2, condition.categ2,
                bias.factor, condition.bias, n, V, Vmes)
  {
# Choose all conditions between an intercept or not

	if (length(categ2) == 1 & length(bias.factor) == 0)
    {
    return(BLUE1(data, fixed2, conti2, categ2, condition.categ2,
                bias.factor, condition.bias, n, V, Vmes))
    } else {
      if (length(bias.factor) == 1 & length(categ2) == 0)
      {
      return(BLUE1(data, fixed2, conti2, categ2, condition.categ2,
                  bias.factor, condition.bias, n, V, Vmes))
      } 
      else {
            return(BLUE2(data, fixed2, conti2, categ2, condition.categ2,
                    bias.factor, condition.bias, n, V, Vmes))
           }
      }
  }


BLUE1 = function(data, fixed2, conti2, categ2, condition.categ2,
                bias.factor, condition.bias, n, V, Vmes)
	{
# Fixed applies only to continuous independent variables
  nc = length(categ2)
  nco = length(conti2)
  nb = length(bias.factor)
  
  if (length(fixed2) == 0)
    {value=data$value}
  if (length(fixed2) > 0)
    {
    value=data$value-data[[fixed2[1]]]
    if (length(fixed2) > 1)
      {
      for (i in 2:length(fixed2))
        {
        value=value-data[[fixed2[i]]]
        }
      }
    }                             
    
  if (nc>0) 
    {
    cf = factorDisj(data,categ2[1])
    }
  if (nb>0) 
    {
    cf = factorDisj(data,bias.factor[1])
    }
   
  e=NULL
  if (length(conti2) > 0)
    {
    e = data[[conti2[1]]]
    if (nco>1)
      {
      for (i in 2:nco)
        {
        e = cbind(e,data[[conti2[i]]])
        }
      }
    }
  
  f = cbind(cf,e)
  V=V+Vmes
  invV=solve(V)
	z=t(f)%*%invV%*%f
  v.p = diag(solve(z))
	p = solve(z)%*%t(f)%*%invV%*%value
  return(list(bias=data.frame(coef=p,var=v.p),drifts=f))
  }

  BLUE2 = function(data, fixed2, conti2, categ2, condition.categ2,
                  bias.factor, condition.bias, n, V, Vmes)
  	{
    nc = length(categ2)
    nco = length(conti2)
    if (length(fixed2) == 0)
      {value=data$value}
    if (length(fixed2) > 0)
      {value=data$value-data[[fixed2]]}
    
    cf = factorDisj(data,categ2[1])
    c = factorDisj(data,categ2[1])[,2:length(levels(data[[categ2[1]]]))]
    if (nc>1)
      {
      for (i in 2:nc)
        {
        cf = cbind(cf,factorDisj(data,categ2[i]))
        c = cbind(c,factorDisj(data,categ2[i])[,2:length(levels(data[[categ2[i]]]))]-cf[,1])
        }
      }
    c = c-cf[,1]
  
    nb = length(bias.factor)
  
    df = factorDisj(data,bias.factor[1])
    d = factorDisj(data,bias.factor[1])[,2:length(levels(data[[bias.factor[1]]]))]
    if (nb>1)
      {
      for (i in 2:nb)
        {
        df = factorDisj(data,bias.factor[i])
        d = cbind(d,factorDisj(data,bias.factor[i])[,2:length(levels(data[[bias.factor[1]]]))]-df[,1])
        }
      }
    d = d-df[,1]
    
    if (length(conti2) > 0)
      {
      e = data[[conti2[1]]]
      if (nco>1)
        {
        for (i in 2:nco)
          {
          e = cbind(e,data[[conti2[i]]])
          }
        }
      f = cbind(1,d,c,e)
      f2 = cbind(1,df,cf,e)
      }
    
    if (length(conti2) == 0)
      {
      f = cbind(1,d,c)
      f2 = cbind(1,df,cf)
      }

    V=V+Vmes
    invV=solve(V)
  	z=t(f)%*%invV%*%f
    v.p = diag(solve(z))
  	p = solve(z)%*%t(f)%*%invV%*%value
  
  # Conditionning the result like Sum of the biases equals 0
  # by default mean of one bias factor equals 0
  # in the case of m multiple categ2orical variables in the drift there should be m-1 conditions
  # indicating the simple arithmetic mean or a value to one of the categories
  # A condition is writen as c("Factor", "Mean", "Value") or c("Factor", "Level", "Value")
  
    ndis = nb+nc
    ndis2 = dim(c)[2] + dim(d)[2]
  
    r = c(rep(-1,dim(d)[2]),rep(0,dim(c)[2]))
    r = rbind(r,c(rep(0,dim(d)[2]),rep(-1,dim(c)[2])))
  
    a = r%*%p[2:(dim(c)[2] + dim(d)[2] + 1)]
  
    v.a=solve(z)
  
    v.d = sum(solve(z)[2:(dim(d)[2]+1),2:(dim(d)[2]+1)])
    v.c = sum(solve(z)[(dim(d)[2]+2):(dim(d)[2]+dim(c)[2]+1),(dim(d)[2]+2):(dim(d)[2]+dim(c)[2]+1)])
  
    return(list(bias=data.frame(coef=c(p[1],a[1],p[2:(dim(d)[2]+1)],a[2],p[(dim(d)[2]+2):(dim(d)[2]+dim(c)[2]+2)]),
                      var=c(v.p[1],v.d,v.p[2:(dim(d)[2]+1)],v.c,v.p[(dim(d)[2]+2):(dim(d)[2]+dim(c)[2]+2)])),
                      drifts=f2)) 
    }


##################### End of internal functions ##########################