##TFM-Pvalue
g++ -I/Users/gtan/src/R-devel/include -DNDEBUG  -I/usr/local/include -I"/Users/gtan/src/R-devel/library/Rcpp/include"   -fPIC  -g -O2  -c ArgumentException.cpp -o ArgumentException.o
g++ -I/Users/gtan/src/R-devel/include -DNDEBUG  -I/usr/local/include -I"/Users/gtan/src/R-devel/library/Rcpp/include"   -fPIC  -g -O2  -c FileException.cpp -o FileException.o
g++ -I/Users/gtan/src/R-devel/include -DNDEBUG  -I/usr/local/include -I"/Users/gtan/src/R-devel/library/Rcpp/include"   -fPIC  -g -O2  -c ParseException.cpp -o ParseException.o
g++ -I/Users/gtan/src/R-devel/include -DNDEBUG  -I/usr/local/include -I"/Users/gtan/src/R-devel/library/Rcpp/include"   -fPIC  -g -O2  -c Matrix.cpp -o Matrix.o
g++ -I/Users/gtan/src/R-devel/include -DNDEBUG  -I/usr/local/include -I"/Users/gtan/src/R-devel/library/Rcpp/include"   -fPIC  -g -O2  -c TFMpvalue.cpp -o TFMpvalue.o
g++ -I/Users/gtan/src/R-devel/include -DNDEBUG  -I/usr/local/include -I"/Users/gtan/src/R-devel/library/Rcpp/include"   -fPIC  -g -O2  -c TFMMain.cpp -o TFMMain.o
g++ -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -o TFMPvalue.so TFMMain.o TFMpvalue.o ArgumentException.o FileException.o ParseException.o Matrix.o -L/Users/gtan/src/R-devel/lib -lR -dylib_file libRblas.dylib:/Users/gtan/src/R-devel/lib/libRblas.dylib -Wl,-framework -Wl,CoreFoundation

library(Rcpp)                                    
mat <- matrix(c(3, 5, 4, 2, 7, 0, 3, 4, 9, 1, 1, 3, 3, 6, 4, 1, 11, 0, 3, 0, 11, 0, 2, 1, 11, 0, 2, 1, 3, 3, 2, 6, 4, 1, 8, 1, 3, 4, 6, 1, 8, 5, 1, 0, 8, 1, 4, 1, 9,0,2,3,9,5,0,0,11,0,3,0,2,7,0,5), nrow = 4, 
                    dimnames = list(c("A","C","G","T")
                                    ))
bg=c(A=0.25, C=0.25, G=0.25, T=0.25)
score=8.77
type="PFM"
pvalue=1e-5
dyn.load("/Users/gtan/Repositories/Bitbucket/TFMPvalue/src/TFMPvalue.so")
.Call("sc2pv", mat, score, bg, type)
.Call("pv2sc", mat, pvalue, bg, type)
.Call("lazyScore", mat, pvalue, bg, type, 1e-5)

./TFMpvalue-fastpvalue -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m MA0045.pfm -s 8.77 -G 1e-5
./TFMpvalue-lazydistrib -a 0.25 -t 0.25 -c 0.25 -g 0.25 -m MA0045.pfm -p 1e-5 -G 1e-5

R CMD build TFMPvalue
R CMD check --as-cran TFMPvalue_0.0.5.tar.gz
R CMD install TFMPvalue_0.0.5.tar.gz

library(TFMPvalue)
library(RUnit)
TFMsc2pv(pwm, score, type="PWM")



