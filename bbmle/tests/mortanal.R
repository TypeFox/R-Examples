library(bbmle)

## goby data in dump format

x <- structure(list(indiv = structure(as.integer(c(20, 77, 79, 21, 
33, 40, 11, 28, 43, 85, 56, 49, 29, 37, 57, 36, 66, 65, 19, 69, 
47, 60, 23, 25, 39, 84, 12, 5, 76, 55, 32, 10, 75, 4, 78, 80, 
86, 48, 54, 22, 18, 61, 41, 74, 68, 14, 53, 45, 30, 17, 62, 3, 
7, 50, 34, 82, 8, 70, 38, 52, 2, 63, 81, 15, 44, 58, 13, 26, 
73, 83, 59, 42, 72, 67, 35, 16, 1, 46, 27, 64, 51, 24, 71, 6, 
9, 31)), .Label = c("f10al1", "f10al2", "f10al3", "f10r1", "f10r2", 
"f11al1", "f11al2", "f11al3", "f11al4", "f11r1", "f11r2", "f11r3", 
"f12al1", "f12al2", "f12al3", "f12al4", "f12al5", "f12r1", "f12r2", 
"f12r3", "f12r4", "f12r5", "f12r6", "f13al1", "f13r1", "f14al1", 
"f14al2", "f14r1", "f14r2", "f15al1", "f15al2", "f15r1", "f15r2", 
"f18al1", "f18al2", "f18r1", "f18r2", "f19al1", "f19r1", "f19r2", 
"f1al1", "f1al2", "f1r1", "f20al1", "f20al2", "f20al3", "f20r1", 
"f20r2", "f20r3", "f2al1", "f2al2", "f2al3", "f2al4", "f2r1", 
"f2r2", "f2r3", "f2r4", "f3al1", "f3al2", "f3r1", "f3r2", "f4al1", 
"f5al1", "f5al2", "f5r1", "f5r2", "f6al1", "f6al2", "f6r1", "f7al1", 
"f7al2", "f7al3", "f7al4", "f7al5", "f7r1", "f7r2", "f7r3", "f7r4", 
"f7r5", "f7r6", "f9al1", "f9al2", "f9al4", "f9r1", "f9r2", "f9r3"
), class = "factor"), group = structure(as.integer(c(5, 5, 5, 
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 
3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)), .Label = c("AL", 
"AL-Rat5th", "AL-RatOv", "R", "R-ALat5th"), class = "factor"), 
    lifespan = as.integer(c(391, 370, 346, 341, 334, 320, 319, 
    317, 314, 307, 295, 260, 30, 10, 397, 380, 364, 355, 352, 
    341, 340, 339, 336, 320, 314, 312, 308, 302, 296, 290, 284, 
    267, 263, 263, 255, 253, 242, 222, 220, 181, 64, 36, 192, 
    192, 189, 186, 183, 181, 180, 176, 173, 171, 170, 169, 166, 
    11, 247, 235, 234, 233, 232, 224, 221, 220, 215, 210, 210, 
    204, 202, 17, 13, 301, 300, 296, 281, 271, 253, 250, 241, 
    239, 232, 221, 220, 214, 33, 30))), .Names = c("indiv", "group", 
"lifespan"), class = "data.frame", row.names = c("1", "2", "3", 
"4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 
"16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", 
"27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
"38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", 
"49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", 
"60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", 
"71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", 
"82", "83", "84", "85", "86"))

mlife <- log(mean(x$lifespan))
Bm0w <- mle2(lifespan~dweibull(scale=exp(llambda),shape=alpha),
             start=list(llambda=mlife,alpha=1),
             data=x)
Bm1w <- mle2(lifespan~dweibull(scale=exp(llambda),shape=alpha),
           start=list(llambda=mlife,alpha=1),
           parameters=list(llambda~group),
           data=x)
Bm2w <- mle2(lifespan~dweibull(scale=exp(llambda),shape=alpha),
             start=list(llambda=mlife,alpha=1),
             parameters=list(llambda~group,alpha~group),
             data=x)           
Bm3w <- mle2(lifespan~dweibull(scale=exp(llambda),shape=alpha),
           start=list(llambda=mlife,alpha=3),
           parameters=list(alpha~group),
           data=x)
anova(Bm0w,Bm1w)
anova(Bm0w,Bm1w,Bm2w)
anova(Bm0w,Bm3w,Bm2w)
AICctab(Bm0w,Bm1w,Bm2w,Bm3w,sort=TRUE,nobs=nrow(x),delta=TRUE)

