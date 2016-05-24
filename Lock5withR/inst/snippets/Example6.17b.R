favstats(~BodyTemp, data = BodyTemp50)
t.test(~BodyTemp, mu = 98.6, data = BodyTemp50)
pval(t.test(~BodyTemp, mu = 98.6, data = BodyTemp50)) # to find the p-value directly

