expect.v <-
function(am, siga, dme.explist){
# expect.v() - expected value of V matrix given siga estimates
    v <- matrix(0,am$l * am$n, am$l * am$n)
    if (any(am$components == "VarE(I)")) {
     v <- v + kronecker(matrix(siga["VarE(I)", ], am$l, am$l), 
                        matrix(dme.explist$vmat[ ,"VarE(I)"],am$n,am$n),
                        make.dimnames=T)
    }
    if (any(am$components == "VarG(Ia)")) {
     v <- v + kronecker(matrix(siga["VarG(Ia)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarG(Ia)"],am$n,am$n),
                        make.dimnames=T)
    }
    if (any(am$components == "VarG(Id)")) {
     v <- v + kronecker(matrix(siga["VarG(Id)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarG(Id)"],am$n,am$n),
                        make.dimnames=T)
    }
    if (any(am$components == "VarG(Ia:a)")) {
     v <- v + kronecker(matrix(siga["VarG(Ia:a)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarG(Ia:a)"],am$n,am$n),
                        make.dimnames=T)
    }
    if (any(am$components == "VarG(Ia:d)")) {
     v <- v + kronecker(matrix(siga["VarG(Ia:d)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarG(Ia:d)"],am$n,am$n),
                        make.dimnames=T)
    }
    if (any(am$components == "VarG(Id:d)")) {
     v <- v + kronecker(matrix(siga["VarG(Id:d)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarG(Id:d)"],am$n,am$n),
                        make.dimnames=T)
    }
    if (any(am$components == "VarGs(Ia)")) {
     v <- v + kronecker(matrix(siga["VarGs(Ia)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarGs(Ia)"],am$n,am$n),
                        make.dimnames=T)
    }
    if (any(am$components == "VarE(M)")) {
     v <- v + kronecker(matrix(siga["VarE(M)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarE(M)"],am$n,am$n),
                        make.dimnemes=T) 
    }
    if (any(am$components == "VarE(C)")) {
     v <- v + kronecker(matrix(siga["VarE(C)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarE(C)"],am$n,am$n),
                        make.dimnemes=T) 
    }
    if (any(am$components == "VarE(M&!C)")) {
     v <- v + kronecker(matrix(siga["VarE(M&!C)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarE(M&!C)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "VarE(M&C")) {
     v <- v + kronecker(matrix(siga["VarE(M&C)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"VarE(M&C)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "VarG(Ma")) {
     v <- v + kronecker(matrix(siga["VarG(Ma)", ], am$l, am$l), 
                        matrix(dme.explist$vmat[ ,"VarG(Ma)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "VarG(Md")) {
     v <- v + kronecker(matrix(siga["VarG(Md)", ], am$l, am$l), 
                        matrix(dme.explist$vmat[ ,"VarG(Md)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "VarG(Ma:a")) {
     v <- v + kronecker(matrix(siga["VarG(Ma:a)", ], am$l, am$l), 
                        matrix(dme.explist$vmat[ ,"VarG(Ma:a)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "VarG(Ma:d")) {
     v <- v + kronecker(matrix(siga["VarG(Ma:d)", ], am$l, am$l), 
                        matrix(dme.explist$vmat[ ,"VarG(Ma:d)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "VarG(Md:d")) {
     v <- v + kronecker(matrix(siga["VarG(Md:d)", ], am$l, am$l), 
                        matrix(dme.explist$vmat[ ,"VarG(Md:d)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "VarGs(Ma")) {
     v <- v + kronecker(matrix(siga["VarGs(Ma)", ], am$l, am$l), 
                        matrix(dme.explist$vmat[ ,"VarGs(Ma)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovE(I,M)")) {
     v <- v + kronecker(matrix(siga["CovE(I,M)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovE(I,M)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovE(M,I)")) {
     v <- v + kronecker(matrix(siga["CovE(M,I)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovE(M,I)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovE(I,M&!C)")) {
     v <- v + kronecker(matrix(siga["CovE(I,M&!C)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovE(I,M&!C)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovE(M&!C,I)")) {
     v <- v + kronecker(matrix(siga["CovE(M&!C,I)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovE(M&!C,I)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Ia,Ma)")) {
     v <- v + kronecker(matrix(siga["CovG(Ia,Ma)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Ia,Ma)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Ma,Ia)")) {
     v <- v + kronecker(matrix(siga["CovG(Ma,Ia)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Ma,Ia)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Id,Md)")) {
     v <- v + kronecker(matrix(siga["CovG(Id,Md)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Id,Md)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Md,Id)")) {
     v <- v + kronecker(matrix(siga["CovG(Md,Id)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Md,Id)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Ia:a,Ma:a)")) {
     v <- v + kronecker(matrix(siga["CovG(Ia:a,Ma:a)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Ia:a,Ma:a)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Ma:a,Ia:a)")) {
     v <- v + kronecker(matrix(siga["CovG(Ma:a,Ia:a)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Ma:a,Ia:a)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Ia:d,Ma:d)")) {
     v <- v + kronecker(matrix(siga["CovG(Ia:d,Ma:d)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Ia:d,Ma:d)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Ma:d,Ia:d)")) {
     v <- v + kronecker(matrix(siga["CovG(Ma:d,Ia:d)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Ma:d,Ia:d)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Id:d,Md:d)")) {
     v <- v + kronecker(matrix(siga["CovG(Id:d,Md:d)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Id:d,Md:d)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovG(Md:d,Id:d)")) {
     v <- v + kronecker(matrix(siga["CovG(Md:d,Id:d)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovG(Md:d,Id:d)"],am$n,am$n),
                        make.dimnemes=T)
    }

    if(any(am$components ==  "CovGs(Ia,Ma)")) {
     v <- v + kronecker(matrix(siga["CovGs(Ia,Ma)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovGs(Ia,Ma)"],am$n,am$n),
                        make.dimnemes=T)
    }
    if(any(am$components ==  "CovGs(Ma,Ia)")) {
     v <- v + kronecker(matrix(siga["CovGs(Ma,Ia)", ], am$l, am$l),
                        matrix(dme.explist$vmat[ ,"CovGs(Ma,Ia)"],am$n,am$n),
                        make.dimnemes=T)
    }

    return(v)
}
