# "exact" p-value
binom.test(60,100);                             
# approximate p-value
z <- ( 0.6 - 0.5 ) / sqrt(0.5 * 0.5/100); z;
2 * (1 - pnorm(z) )                         

# approximate p-value with continuity correction
z <- ( 0.595 - 0.5 ) / sqrt(0.5 * 0.5 / 100); z;    #0.595 = 59.5/100
2 * (1 - pnorm(z) )                         

# R can automate the approximate version too:
prop.test(60,100)            # uses continuity correction by default
prop.test(60,100,correct=FALSE)     # turn off continuity correction
