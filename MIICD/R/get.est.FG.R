
get.est.FG <-
  function( x , data ,  status , trans , cens , formula ){
  df1 <- data.frame( data , 'ns' = x )
  FGR1<-paste("FGR( Hist( time = ns , event = ",status,", cens.code = '", cens ,"' )", deparse(formula) ,", data = df1 , cause = '",trans,"' )" ,sep='' )
  FGR2<-eval(parse(text=FGR1)) 
  vc1<-FGR2$crrFit$var
  #vc1<-vcov(FGR2$crrFit)
  coef1<-FGR2$crrFit$coef
  return( list( sigma = vc1 , beta = coef1 ) )
}


