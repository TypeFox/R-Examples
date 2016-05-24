#-----------------------------------------------------------------------
#     Copyright (C) 2012-2014  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
# test wrap, as and NA values (mostly of myself, should be enhanced...)
#
# This file cannot be set in the tests directory:
#
# Error in compileCode(f, code, language = language, verbose = verbose) : 
#  Compilation ERROR, function(s)/method(s) not created! Error in file(filename, "r", encoding = encoding) : 
#  cannot open the connection
# Calls: local ... eval.parent -> eval -> eval -> eval -> eval -> source -> file
# In addition: Warning message:
# In file(filename, "r", encoding = encoding) :
#  cannot open file 'startup.Rs': No such file or directory
# Execution halted
# Calls: cxxfunction -> compileCode
# Execution halted
#
if (require("inline"))
{
#library(rtkpp)
body <- '
  // wrap SEXP
  NumericMatrix RData1(tab1);
  NumericMatrix RData2(tab2);
  NumericMatrix RData3(tab3);

  // wrap Rcpp matrix
  RMatrix<double> data1(RData1);
  // start tests
  for (int i=data1.beginRows(); i< data1.endRows(); ++i)
  {
    for (int j= data1.beginCols(); j < data1.endCols(); ++j)
    {
      if (Arithmetic<Real>::isNA(data1(i,j)))
      { data1(i,j)= 100*i+j;}
    }
  }
  data1(0,0) = Arithmetic<double>::NA();
  RDataHandler handler(RData2, "model1", "gaussian_sjk");
  handler.addData(RData3,  "model2", "gaussian_sj");

  RMatrix<double> data2;
  int nbVar;
  handler.getData("model1", data2, nbVar);

  data2(0,0) = 10;

  Rcpp::List l;
  l.push_back(tab1, "1");
  data1(1,1) = 11;
  l.push_back(tab1, "2");
  data1(2,2) = 22;
  l.push_back(tab2, "3");

  RMatrix<double> m = as< RMatrix<double> >(tab1);
  m(1,1) = Arithmetic<double>::NA();
  // return final output
  return l;
'

fx <- cxxfunction( signature(tab1 = "matrix", tab2 = "matrix", tab3 = "matrix" ) , body, plugin = "rtkpp", verbose = TRUE )


data(iris)
mat1 <- as.matrix(iris[1:4])[1:10,]
mat2 <- as.matrix(iris[1:4])[11:20,]
mat3 <- as.matrix(iris[1:4])[21:30,]
mat1[10,4] = NA;
mat2[10,4] = NA;
mat3[10,4] = NA;

xem <- fx( mat1, mat2, mat3)

mat1
mat2
mat3
xem
}
