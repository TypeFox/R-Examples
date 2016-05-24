setMethod("get.results","resultsClassTest",function(object)
{
 p.val=slot(object,"p.value")
 delta=slot(object,"delta")
 genelist=slot(object,"class.genes")
 return(list(p.value=p.val,delta=delta,class.genes=genelist))
})

setMethod("get.results","resultsIndTest",function(object)
{
 p.val=slot(object,"p.values")
 test.stat=slot(object,"d")
 gene.names=names(p.val)
 p=length(p.val)
 out=data.frame(d=test.stat,p.value=p.val)
 return(data.frame(out[order(p.val),]))
})

setMethod("get.results","resultsModTest",function(object)
{
 p.val=slot(object,"p.value")
 N=slot(object,"N")
 mod1=slot(object,"modules1")
 mod2=slot(object,"modules2")
 return(list(p.value=p.val,N=N,modules1=mod1,modules2=mod2))
})

