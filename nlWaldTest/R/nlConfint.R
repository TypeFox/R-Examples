nlConfint <-
function(obj=NULL,texts,level = 0.95,coeff=NULL,Vcov=NULL,
                      df2=NULL,x=NULL)
{
  if (!is.null(obj))
  {
    co = try(coef(obj),silent =T)
    cond = attr(co,"condition")
    if (is.null(coeff)&&(is.null(cond))) 
      coeff = co
    vc = try(vcov(obj),silent=T)
    cond2 = attr(vc,"condition")
    if (is.null(Vcov)&&(is.null(cond2))) 
      Vcov = vc
  }
  if (is.null(coeff)) 
  {
    if (is.null(obj))
      mess = "Both  'obj' and 'coeff' are missing"
    else
    {
      clm = class(obj)
      part1 = "There are no coef() methods for model objects of class \""
      mess = paste0(part1,clm,"\".\nInput the 'coeff' parameter.")
    }
    stop(mess)
  }
  
  if (is.null(Vcov)) 
  {
    if (is.null(obj))
      mess = "Both  'obj' and 'Vcov' are missing"
    else
    {
      clm = class(obj)
      part1 = "There are no vcov() methods for model objects of class \""
      mess = paste0(part1,clm,"\".\nInput the 'Vcov' parameter.")
    }
    stop(mess)
  }
  if (length(texts)>1)
    kkk = texts[1]
  else
    kkk = strsplit(texts[1],";")[[1]]
  # kkk = strsplit(kkk[1],"=")[[1]]
  kkkfl = as.formula(paste("~",kkk[1]))
  vvss = setdiff(all.vars(kkkfl),"x")
  texts = .smartsub(vvss,"b",texts)
  if (length(texts) > 1)
    ltext0 = texts
  else
    ltext0=strsplit(texts,";")[[1]]
  texts1=gsub("[","",texts,fixed=T)
  texts1=gsub("]","",texts1,fixed = T)
  if (length(texts1)>1)
    ltext = texts1
  else
    ltext=strsplit(texts1,";")[[1]]
  r = length(ltext)
  n = length(coeff)
  namess = paste0("b",1:n)
  for (j in 1L:n)
    assign(namess[j],coeff[j])
  if (!is.null(x))
  {nx=length(x)
  namesx = paste0("x",1:nx)
  for (j in 1L:nx)
    assign(namesx[j],x[j])  
  }
  grad=c()
  hess=c()
  for (i in 1L:r)
  {
    fli <- as.formula(paste("~",ltext[i]))
    z=try(deriv(as.formula(fli),namess),silent = T)
    if (class(z)=="try-error")
    {
      tei = as.character(i)
      tri2 = ", numerical derivatives were used in delta-method"
      wate = paste0("Note: For function ",i,tri2)
      message(wate)
      ez=numericDeriv(quote(eval(parse(text=ltext[i]))),namess)
    }
    else
      ez=eval(z)
    hessj=attr(ez,"gradient")
    grad=rbind(grad,ez[1])
    hess=rbind(hess,hessj)
  }
  Rb=grad
  ddd=hess%*%Vcov%*%t(hess)
  matr = chol2inv(chol(ddd))
  ses=sqrt(diag(ddd))
  trydf = identical(df2,T)
  if (trydf)
  {
    isdf = try(df.residual(obj),silent =T)
    df2 = isdf
    if (is.null(df2))
    {
      wn="Note: Failed to extract df for denomenator; z-intervals applied"
      message(wn)
    }
  }
  .getint(as.numeric(Rb),ltext0,ses, level, df = df2) 
  
}
