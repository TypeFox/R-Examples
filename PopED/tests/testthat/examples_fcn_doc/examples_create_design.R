library(PopED)

xt1 <- list(c(1,2,3),c(1,2,3,4))
xt4 <- list(c(1,2,3,4,5),c(1,2,3,4))
xt2 <- rbind(c(1,2,3,4),c(1,2,3,4))
xt3 <- c(1,2,3,4)

design_1 <- create_design(xt=xt1,groupsize=20)
design_2 <- create_design(xt=xt4,groupsize=20)
design_3 <- create_design(xt=xt2,groupsize=20)
design_4 <- create_design(xt=xt3,groupsize=20)

design_5 <- create_design(xt=xt3,groupsize=20,m=3)

design_6 <- create_design(xt=xt1,groupsize=20,model_switch=ones(2,4))

design_7 <-create_design(xt=xt1,groupsize=20,a=c(2,3,4))
design_8 <-create_design(xt=xt1,groupsize=20,a=rbind(c(2,3,4),c(4,5,6)))
design_9 <-create_design(xt=xt1,groupsize=20,a=list(c(2,3,4,6),c(4,5,6)))
design_10 <-create_design(xt=xt1,groupsize=20,a=list(c(2,3,4),c(4,5,6)))

design_11 <-create_design(xt=c(0,1,2,4,6,8,24),
                         groupsize=50,
                         a=c(WT=70,DOSE=1000))

design_12 <-create_design(xt=c(0,1,2,4,6,8,24),
                         groupsize=50,
                         a=c(WT=70,DOSE=1000),m=2)

design_13 <-create_design(xt=c(0,1,2,4,6,8,24),
                         groupsize=50,
                         a=list(c(WT=70,DOSE=1000),c(DOSE=90,WT=200,AGE=45)),m=2)

design_14 <-create_design(xt=c(0,1,2,4,6,8,24),
                         groupsize=50,
                         a=list(list(WT=70,DOSE=1000),list(DOSE=90,WT=200,AGE=45)),m=2)

design_15 <-create_design(xt=xt4,
                          groupsize=c(50,20),
                          a=rbind(c("DOSE"=2,"WT"=3,"AGE"=4),
                                  c(4,5,6)))

