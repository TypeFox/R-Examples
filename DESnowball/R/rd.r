fs.rd <-
function(x,
	 method="mcd",
	 df=1)
{
    rd <- robust.distance.signed(x,method=method)
    pval <- pchisq(rd,df=df,lower.tail=F)
    full.list <- data.frame(rd=rd,
                            pval=pval,
                            positive=x$positives)
    row.names(full.list) <- row.names(x)
    full.list[order(pval,decreasing=F),]
}
