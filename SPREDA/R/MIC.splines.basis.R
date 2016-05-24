MIC.splines.basis <-
function(x, df = NULL, knots = NULL,
                           boundary.knots=NULL,type="Ms",degree = 3,delta=0.01,eq.alloc=F)
{
  if(is.null(df)&is.null(knots))
  {
    stop("either df or knots needs to be supplied")
  }
  
  if(!is.null(df))
  {
    cc=df-degree+1
    knots=quantile(x,1:(cc-1)/cc)
  }
  
  if(is.null(boundary.knots))
  {
    boundary.knots=range(x)
    boundary.knots[2]=boundary.knots[2]+0.0001
  }
  
  if(eq.alloc==T)
  {
    if(is.null(df))
    {
      stop("df needs to be supplied for equal allocation")
    }
    cc=df-degree+1
    knots=boundary.knots[1]+(1:(cc-1))*(boundary.knots[2]-boundary.knots[1])/cc
  }
  
  tt=c(rep(boundary.knots[1],degree),knots,rep(boundary.knots[2],degree))
  nn=length(x)
  mm=degree+length(knots)
  mat=matrix(0, nrow=nn,ncol=mm)
  
  #M splines basis
  if(type=="Ms")
  {
    for(i in 1:nn)
    {
      for(j in 1:mm)
      {
        mat[i,j]=m.spline.x(x[i],tt,j,k=degree)
      }
    }
  }
  
  #I spline
  if(type=="Is")
  {
    for(j in 1:mm)
    {
      mat[,j]=i.spline.x(x,tt,j,k=degree,delta=delta,Cs=F)
    }
    mat=sweep(mat,2,colMeans(mat),"-")
  }
  
  #C spines
  if(type=="Cs")
  {
    for(j in 1:mm)
    {
      mat[,j]=i.spline.x(x,tt,j,k=degree,delta=delta,Cs=T)
    }
    x1=x-mean(x)
    x1=x1/sqrt(sum(x1*x1))
    mat=sweep(mat,2,colMeans(mat),"-")
    for(j in 1:mm)
    {
      mat[,j]=mat[,j]-x1*sum(mat[,j]*x1)
    }
    mat=cbind(x1,mat)
  }
  ###return results
  res=list(mat=mat,x=x, df=df, knots=knots,boundary.knots=boundary.knots,
           type=type,degree=degree,delta=delta)
  attr(res,"class")="MICsplines"
  return(res)
}
