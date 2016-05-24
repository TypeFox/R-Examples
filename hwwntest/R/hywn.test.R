hywn.test <-
function(x, filter.number=10, family="DaubExPhase"){


p.val.a <- hwwn.test(x)$p.value
p.val.b <- genwwn.test(x, filter.number=filter.number, family=family)$p.value
p.val.c <- bartlettB.test(x)$p.value
p.val.d <- Box.test(x, lag=20)$p.value

p.val <- min(p.adjust(c(p.val.a, p.val.b, p.val.c, p.val.d), method="bonferroni"))

ll <- list(p.value=p.val, method="Hybrid Test")
class(ll) <- "htest"
return(ll)
}
