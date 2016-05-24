row.perc( xtabs(~t2d + genotype, fusion1m) )
chisq.test( xtabs(~t2d + genotype, fusion1m) )
chisq.test( xtabs(~t2d + (Tdose >=1), fusion1m) )
chisq.test( xtabs(~t2d + (Tdose <=1), fusion1m) )
