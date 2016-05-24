library(survey)
set.seed(42)

df<-data.frame(x=exp(rnorm(1000)))
df$y<-round(df$x,1)
ddf<-svydesign(id=~1,data=df)
rdf<-as.svrepdesign(ddf)

SE(svyquantile(~x,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE))

SE(svyquantile(~x,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE))


svyquantile(~y,ddf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE,ties="rounded",interval.type="betaWald")

svyquantile(~y,rdf, c(0.01,0.1,0.5,0.9,0.99),ci=TRUE)



