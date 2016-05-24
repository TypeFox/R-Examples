# 9 values of M, 2000 MCMC cycles:  1 minute  20100518
# Default for breast cancer data: matrices are 2001 x 21
bf1 <- function(data,seed=1,ncycles=2000,
                d=c(.1,.1,0,1000),K=10,burnin=1000)
  {
   if (!is.numeric(data) || !is.matrix(data))
     stop("data must be numeric matrix")
   psi.hat <- data[,1]; se.psi.hat <- data[,2]
   if (!all(se.psi.hat >= 0)) {
      stop("negative standard error in column 2 of data.\n")
   }
   if (!is.numeric(d) || length(d) != 4){
     stop("hyperparameter must be numeric of length 4")
   }
   else if (d[1] <= 0 || d[2] <= 0 || d[4] <= 0)
   {
      stop("Gamma shape, Gamma scale, and normal variance hyperparameters
           must all be greater than zero\n")
   }   
   set.seed(seed)
   M00.25 <- dirichlet.o(data, ncycles=ncycles, M=.25,d=d,K=K)$chain
#   system("echo 'done with run' | mail -s 'bsp' XXX")
   set.seed(seed);
   M00.5 <- dirichlet.o(data, ncycles=ncycles, M=.5,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   set.seed(seed);
   M001  <- dirichlet.o(data, ncycles=ncycles, M=1,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   set.seed(seed);
   M002  <- dirichlet.o(data, ncycles=ncycles, M=2,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   set.seed(seed);
   M004  <- dirichlet.o(data, ncycles=ncycles, M=4,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   set.seed(seed);
   M008  <- dirichlet.o(data, ncycles=ncycles, M=8,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   set.seed(seed);
   M016  <- dirichlet.o(data, ncycles=ncycles, M=16,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   set.seed(seed);
   M032  <- dirichlet.o(data, ncycles=ncycles, M=32,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   set.seed(seed);
   M064  <- dirichlet.o(data, ncycles=ncycles, M=64,d=d,K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
#  set.seed(seed);
#  M128  <- dirichlet.o(data, M=128, ncycles=ncycles,  K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
#  set.seed(seed);
#  M256  <- dirichlet.o(data, M=256, ncycles=ncycles,  K=K)$chain
#  system("echo 'done with run' | mail -s 'bsp' dburrdoss@gmail.com")
   rowl <- burnin+1; rowu <- ncycles+1
   xxx <- list(M00.25=M00.25[-(1:rowl),], M00.5=M00.5[-(1:rowl),],
     M001=M001[-(1:rowl),], M002=M002[-(1:rowl),],
     M004=M004[-(1:rowl),], M008=M008[-(1:rowl),],
     M016=M016[-(1:rowl),], M032=M032[-(1:rowl),],
     M064=M064[-(1:rowl),])
   invisible(xxx)
 }
#put next line in the function, maybe
#save(M00.25, M00.5, M001, M002,
#     M004, M008, M016, M032, M064,
#        file="breast-rdat-seed0", compress=TRUE)

bf2 <- function(chain.list)
{
#    aaa <- names(chain.list)
    nM <- length(chain.list)
    if (!is.list(chain.list) || nM != 9) {stop("'chain.list' must be a
     list of 9 matrices of MCMC output")}
    dOut <- dim(chain.list[[1]])
    ncycles <- dOut[1] - 1
    nparam <- dOut[2]
    for (ii in 1:nM) 
      {
        dOutii <- dim(chain.list[[ii]])
        if (dOutii[1] -1 != ncycles || dOutii[2] != nparam)
          {stop("all matrices in chain.list must have same size\n")}
      }   
    mat00.25 <- chain.list[1]
    mat00.5 <- chain.list[2]
    mat001 <- chain.list[3]
    mat002 <- chain.list[4]
    mat004 <- chain.list[5]
    mat008 <- chain.list[6]
    mat016 <- chain.list[7]
    mat032 <- chain.list[8]
    mat064 <- chain.list[9]
    cc <- numeric(9)
    # (lr of ) 4/2*2/1*1/.5*.5/.25 = lr of 4/.25
    cc[1] <- 
      mean(lr(Mnew=016/4, Mold=008/4, cc=1, mat.list=mat002)) *
      mean(lr(Mnew=008/4, Mold=004/4, cc=1, mat.list=mat001)) *
      mean(lr(Mnew=004/4, Mold=002/4, cc=1, mat.list=mat00.5)) *
      mean(lr(Mnew=002/4, Mold=001/4, cc=1, mat.list=mat00.25))
    cc[2] <-
      mean(lr(Mnew=016/4, Mold=008/4, cc=1, mat.list=mat002)) *
      mean(lr(Mnew=008/4, Mold=004/4, cc=1, mat.list=mat001)) *
      mean(lr(Mnew=004/4, Mold=002/4, cc=1, mat.list=mat00.5))
    # lr of 4/1
    cc[3] <-
      mean(lr(Mnew=016/4, Mold=008/4, cc=1, mat.list=mat002)) *
      mean(lr(Mnew=008/4, Mold=004/4, cc=1, mat.list=mat001))
    # lr of 4/2
    cc[4] <-
      mean(lr(Mnew=016/4, Mold=008/4, cc=1, mat.list=mat002))
    # lr of 4/4 = 1
    cc[5] <- 1
    # lr of 4/8
    cc[6] <-
      mean(lr(Mnew=016/4, Mold=032/4, cc=1, mat.list=mat008))
    # lr of 4/16
    cc[7] <-
      mean(lr(Mnew=016/4, Mold=032/4, cc=1, mat.list=mat008)) * 
      mean(lr(Mnew=032/4, Mold=064/4, cc=1, mat.list=mat016))
    # lr of 4/32
    cc[8] <-
      mean(lr(Mnew=016/4, Mold=032/4, cc=1, mat.list=mat008)) * 
      mean(lr(Mnew=032/4, Mold=064/4, cc=1, mat.list=mat016)) * 
      mean(lr(Mnew=064/4, Mold=128/4, cc=1, mat.list=mat032))
    # lr of 4/64
    cc[9] <-
      mean(lr(Mnew=016/4, Mold=032/4, cc=1, mat.list=mat008)) * 
      mean(lr(Mnew=032/4, Mold=064/4, cc=1, mat.list=mat016)) * 
      mean(lr(Mnew=064/4, Mold=128/4, cc=1, mat.list=mat032)) * 
      mean(lr(Mnew=128/4, Mold=256/4, cc=1, mat.list=mat064)) 
    cc
  }

