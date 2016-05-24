#
# 13_07_03 13_07_04
#
# defining objects for the /rbmn/ package
#
set.seed(1234);
####################################################
#
# <1> DEFINING SOME INITIAL OBJECTS
#
####################################################
rbmn0chain.01 <-
     list(names = c("A", "B", "C"),
          mu    = c(  0,   0,   0),
          sigma = c(  1,   2,   3),
          corre = c(0.5, 0.5), 
          roots = c("A"), colliders = NULL
         );
####################################################
rbmn0chain.02 <-
     list(names = letters[1:5],
          mu    = 1:5,
          sigma = 5:1,
          corre = rep(1/3,4), 
          roots = c("a","d"), colliders = c("c")
         );
####################################################
rbmn0chain.03 <- generate8chain(10);
####################################################
rbmn0adja.04 <- matrix(c(rep(0,4),1,1,rep(0,4),1,
                      rep(0,5),1,1,0,1,rep(0,5)),
                    5,dimnames=
                    list(c("1.1","1.2","2.1","2.2","C"),
                         c("1.1","1.2","2.1","2.2","C")));
####################################################
#
# <2> TRANSFORMING THEM INTO OTHER TYPES
#
####################################################
rbmn0nbn.01 <- chain2nbn(rbmn0chain.01);
rbmn0nbn.02 <- chain2nbn(rbmn0chain.02);
rbmn0nbn.03 <- chain2nbn(rbmn0chain.03);
rbmn0nbn.04 <-  adja2nbn( rbmn0adja.04);
####################################################
rbmn0adja.01 <- adja4nbn(rbmn0nbn.01);
rbmn0adja.02 <- adja4nbn(rbmn0nbn.02);
rbmn0adja.03 <- adja4nbn(rbmn0nbn.03);
####################################################
rbmn0mn.01 <- nbn2mn(rbmn0nbn.01);
rbmn0mn.02 <- nbn2mn(rbmn0nbn.02);
rbmn0mn.03 <- nbn2mn(rbmn0nbn.03);
rbmn0mn.04 <- nbn2mn(rbmn0nbn.04);
####################################################
rbmn0gema.01 <- nbn2gema(rbmn0nbn.01);
rbmn0gema.02 <- nbn2gema(rbmn0nbn.02);
rbmn0gema.03 <- nbn2gema(rbmn0nbn.03);
rbmn0gema.04 <- nbn2gema(rbmn0nbn.04);
####################################################
#
# <3> PREPARING A CROSSED STRUCTURE
#
####################################################
uu <- vv <- rbmn0adja.01;
ss <- c("T","L","A"); cc <- c("L","F","B");
dimnames(uu) <- list(ss,ss);
dimnames(vv) <- list(cc,cc);
uun <- adja2nbn(uu);
vvn <- adja2nbn(vv);
ww <- as.vector(outer(ss,cc,paste,sep=""));
rbmn0crarc.05 <- arcs4nbn1nbn(uun,vvn,nona=ww);
rbmn0nbn.05 <- crossed4nbn1nbn(uun,vvn,nona=ww);
rm(list=c("uu","vv","cc","ss","uun","vvn","ww"));
#
rbmn0adja.05 <- adja4nbn(rbmn0nbn.05);
rbmn0gema.05 <- nbn2gema(rbmn0nbn.05);
rbmn0mn.05   <- nbn2mn (rbmn0nbn.05);
