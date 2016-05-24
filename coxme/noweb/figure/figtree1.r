pdf('figtree1.pdf', width=7, height=6)
ptplot(y ~ x1 + (x3 + x4)* x2)
dev.off()
