 Sys.getenv("R_LIBS")
 library()
 require("GPArotation")
 search()
 Sys.info()

 fuzz <- 1e-6 
 all.ok <- TRUE  

 data(WansbeekMeijer, package="GPArotation")

 fa.none  <- factanal(factors=2,  covmat=NetherlandsTV, rotation="none")

 tst <- t(matrix(c(
       0.6972803, -0.3736554,
       0.7774628, -0.3184149,
       0.6832300, -0.3620428,
       0.6612198,  0.2361132,
       0.6972393,  0.3026050,
       0.7100285,  0.4059509,
       0.6353584,  0.3526947
      ), 2, 7))

 if( fuzz < max(abs(fa.none$loadings - tst))) {
    cat("Calculated value is not the same as test value in test WansbeekMeijer 1. Value:\n")
    print(fa.none$loadings, digits=18)
    cat("difference:\n")
    print(fa.none$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 



 fa.varimax <- GPFoblq(fa.none$loadings, method="varimax", normalize=TRUE)
 
# with eps=1e-8
# tst <- t(matrix(c(
#	0.229695829694226694, -0.757005882905721683,
#	0.325474298411086493, -0.774533969509160203,
#	0.227951538606475851, -0.738861531224136225,
#	0.634850649690308022, -0.299876110481063607,
#	0.707312661165822032, -0.278246783076943283,
#	0.789359884149245072, -0.214120439603779994,
#	0.698885205896135120, -0.199081171877497243
#	), 2, 7))

# with eps=1e-5
 tst <- t(matrix(c(
  0.229698038368303409, -0.757005212686898243,
  0.325476558225504142, -0.774533019824047542,
  0.227953694341768043, -0.738860866094951829,
  0.634851524619887475, -0.299874258087383661,
  0.707313472988376213, -0.278244719250824557,
  0.789360508873491518, -0.214118136377292989,
  0.698885786741510029, -0.199079132641678647
  ), 2, 7))

 if( fuzz < max(abs(fa.varimax$loadings - tst))) {
    cat("Calculated value is not the same as test value in test WansbeekMeijer 2. Value:\n")
    print(fa.varimax$loadings, digits=18)
    cat("difference:\n")
    print(fa.varimax$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 



 fa.oblimin <- GPFoblq(fa.none$loadings, method="oblimin", normalize=TRUE)

# with eps=1e-8
# tst <- t(matrix(c(
#      -0.0244898894997362740, -0.8055076884898763057,
#	0.0821776433220552660, -0.7883517482514345032,
#      -0.0194442483441249758, -0.7847120136813017233,
#	0.6350106056917923514, -0.1038114236654337219,
#	0.7293893902400611085, -0.0495156037400738894,
#	0.8517915457391848078,  0.0588983480418694277,
#	0.7504355940804637859,  0.0408946221245683056
#	), 2, 7))

# with eps=1e-5
 tst <- t(matrix(c(
  -0.0244886312423446446, -0.8055069385602275922,
   0.0821788889356081659, -0.7883509906546982693,
  -0.0194430219824419312, -0.7847112821295906260,
   0.6350108529538124325, -0.1038111848933331444,
   0.7293895650539216069, -0.0495153948664520185,
   0.8517915670863017708,  0.0588984825074335624,
   0.7504356301074717184,  0.0408947509009953206
   ), 2, 7))

 if( fuzz < max(abs(fa.oblimin$loadings - tst))) {
    cat("Calculated value is not the same as test value in test WansbeekMeijer 3. Value:\n")
    print(fa.oblimin$loadings, digits=18) 
    cat("difference:\n")
    print(fa.oblimin$loadings - tst, digits=18)
    all.ok <- FALSE  
    } 


cat("tests completed.\n")

if (! all.ok) stop("some tests FAILED")

