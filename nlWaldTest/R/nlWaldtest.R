nlWaldtest <-
function(obj=NULL,texts,rhss=NULL,coeff=NULL,Vcov=NULL,
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
  kkk = strsplit(kkk[1],"=")[[1]]
  kkkfl = as.formula(paste("~",kkk[1]))
  vvss = setdiff(all.vars(kkkfl),"x")
  texts = .smartsub(vvss,"b",texts)
  texts1=gsub("[","",texts,fixed=T)
  texts1=gsub("]","",texts1,fixed = T)
  if (length(texts1)>1)
    ltext = texts1
  else
    ltext=strsplit(texts1,";")[[1]]
  r = length(ltext)
  for (i in 1L:r)
  {
    ltexti = ltext[i]
    if(grepl("=",ltexti))
      ltext[i] = paste0(gsub("=","-(",ltexti),")")
  }
  
  n = length(coeff)
  namess = paste0("b",1:n)
  if (is.null(rhss))
    rhss=rep(0,r)
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
      wate = paste0("Note: For restriction ",i,tri2)
      message(wate)
      ez=numericDeriv(quote(eval(parse(text=ltext[i]))),namess)
    }
    else
      ez=eval(z)
    hessj=attr(ez,"gradient")
    grad=rbind(grad,ez[1])
    hess=rbind(hess,hessj)
  }
  Rb=grad-rhss
  ddd=hess%*%Vcov%*%t(hess)
  matr = chol2inv(chol(ddd))
  stat0 = as.numeric(t(Rb)%*%matr%*%Rb)
  df1 = r
  trydf = identical(df2,T)
  if (trydf)
  {
    isdf = try(df.residual(obj),silent =T)
    df2 = isdf
    if (is.null(df2))
    {
      wn="Note: Failed to extract df for denomenator. Chi-square test is applied"
      message(wn)
    }
  }
  if (r==1)
    gramma = "a restriction"
  else
    gramma ="restrictions"
  if (is.null(df2)) 
  {
    method <- paste0("Wald Chi-square test of ",gramma," on model parameters")
    LM = stat0
    names(LM) <- "Chisq"
    df <- r
    names(df) <- "df"
    RVAL <- list(statistic = LM, parameter = df, method = method,
                 p.value = pchisq(LM, df, lower.tail = FALSE))
  }  
  else
  { stat1 = stat0/r
    p2 <- pf(stat1, df1, df2,lower.tail =F)
    method <- paste0("Wald F test of ",gramma," on model parameters")
    LM <- stat1
    names(LM) <- paste("F")
    df1 <- r
    names(df1) <- "df1"
    names(df2) <- "df2"
    RVAL <- list(statistic = LM, parameter = c(df1,df2), method = method,
                 p.value = p2)
  }
  if (!is.null(obj))
    RVAL$data.name = deparse(substitute(obj))
  else
    RVAL$data.name = "model coefs and vcov were supplied separately"
  class(RVAL) <- "htest"
  return(RVAL) 
}
