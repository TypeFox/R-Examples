# one-sided test manually and using binom.test()
1-pbinom(15,50,1/6);
binom.test(16,50,1/6,alternative="greater");
