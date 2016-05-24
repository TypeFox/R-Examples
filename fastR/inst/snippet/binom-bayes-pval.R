1- pbeta(0.5,38+1,62+1)          # 1-sided Bayesian p-value
binom.test(38,100,alt="less")$p.value      # for comparison
prop.test(38,100,alt="less")$p.value       # for comparison
