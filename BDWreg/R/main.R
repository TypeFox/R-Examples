# main function
bdw<-function(data               , formula = NA      ,
              reg.q =FALSE       , reg.b = FALSE     ,
              logit = TRUE       ,
              initial = c(.5,1)  , iteration = 25000 ,
              v.scale = 0.1      ,
              RJ     = FALSE     ,
              dist.q = imp.d     , dist.b = imp.d    ,
              q.par = c(0,0)     , b.par = c(0,0)    ,
              penalized = FALSE  ,
              dist.l = imp.d     ,  l.par  = c(0,0)  ,
              bi.period =.25     , cov  = 1          ,
              sampling = c('bin') ,
              est = Mode         ,
              fixed.l= -1        , ...
)
{
  chain  =  par.bayesian.tot.dw(data = data            ,formula = formula    ,
                                para.q = reg.q         , para.b = reg.b      ,
                                q.par = q.par          , b.par = b.par       ,
                                dist.q = dist.q        , dist.b = dist.b     ,
                                dist.l = dist.l        , l.par = l.par       ,
                                initials = initial     , burn.in=bi.period   ,
                                iterations = iteration , fixed.l = fixed.l   ,
                                v.scale=v.scale        , cov.m = cov         ,
                                sampling = sampling    , penalized = penalized,
                                logit = logit          , RJ = RJ
  )
  res  =  apply(chain$chain,2,est)
  out  = list(res=res , chain=chain)
  class(out) = 'bdw'
  return(out)
}

