context("create design, design space and PopED databases")

test_that("designs are created correctly using create_design()", {
  
  source("examples_fcn_doc/examples_create_design.R")
  
  expect_true(all(size(design_1$xt)==c(2,4)))  
  expect_true(is.na(design_1$xt[1,4]))
  
  expect_true(all(size(design_2$xt)==c(2,5)))
  expect_true(is.na(design_2$xt[2,5]))
  
  expect_true(all(size(design_5$xt)==c(3,4)))
  expect_true(!is.null(colnames(design_5$xt)))
  expect_true(!is.null(rownames(design_5$xt)))
  
  expect_true(!is.null(rownames(design_5$xt)))
  
  expect_true(!any(is.na(design_6$model_switch)))
  
  expect_null(colnames(design_7$a))
  expect_true(all(size(design_7$a)==c(2,3)))
  
  expect_true(all(size(design_9$a)==c(2,4)))
  expect_null(colnames(design_9$a))
  expect_true(is.na(design_9$a[2,4]))
  
  expect_equivalent(design_10,design_8)
  
  expect_true(!is.null(colnames(design_11$a)))
  
  expect_true(!is.null(colnames(design_12$a)))
  expect_true(all(size(design_12$a)==c(2,2)))
  
  expect_true(all(size(design_13$a)==c(2,3)))
  expect_true(is.na(design_13$a[1,3]))
    
  expect_equivalent(design_14,design_13)
  
  expect_true(!is.null(colnames(design_15$a)))
  
  expect_error(create_design(xt=xt1,groupsize=20,m=3))
  
  
})

test_that("design spaces are created correctly using create_design_space()", {
  
  source("examples_fcn_doc/examples_create_design_space.R")
  
  expect_equivalent(ds_1$design,design_1)
  expect_equivalent(ds_1$design_space$mina,design_1$a)
  expect_equivalent(ds_1$design_space$maxa,design_1$a)
  expect_equivalent(ds_1$design_space$minxt,design_1$xt)
  expect_equivalent(ds_1$design_space$maxxt,design_1$xt)
  expect_equivalent(ds_1$design_space$mingroupsize,design_1$groupsize)
  expect_equivalent(ds_1$design_space$maxgroupsize,design_1$groupsize)
  expect_true(length(unique(c(ds_1$design_space$grouped_a)))==length(ds_1$design_space$grouped_a))
  expect_true(length(unique(c(ds_1$design_space$grouped_xt)))==length(ds_1$design_space$grouped_xt))
    
  expect_true(all(ds_2$design_space$maxxt==10))
  expect_true(all(ds_2$design_space$minxt==0))
  expect_true(all(ds_2$design_space$maxni==10))
  
  expect_error(create_design_space(design_1,maxni=10,minni=11))
  expect_error(create_design_space(design_1,minni=15)) 
  expect_error(create_design_space(design_1,maxni=10,mingroupsize=30))
  expect_error(create_design_space(design_1,maxni=10,mingroupsize=20))
  expect_error(create_design_space(design_1,maxni=10,mingroupsize=20,maxxt=10))
  expect_error(create_design_space(design_1,maxni=10,mingroupsize=20,minxt=0))
    
  expect_equivalent(ds_3$design_space$mintotgroupsize,40)
    
  expect_equivalent(ds_4$design_space$maxa,rbind(c(100,2000),c(100,2000)))
    
  expect_equivalent(ds_5$design_space$mina,rbind(c(10,20),c(10,20)))
    
  expect_true(length(ds_6$design_space$x_space)==4)  
  expect_true(length(ds_6$design_space$x_space[[1,1]])==1)  
    
  expect_true(length(ds_7$design_space$x_space)==4)  
  expect_true(length(ds_7$design_space$x_space[[1,2]])==length(seq(100,400,by=20)))  
  
  expect_error(create_design_space(design_2,x_space=list(SEX=c(0,2),DOSE_discrete=seq(100,400,by=20)))) 
    
  expect_equivalent(ds_8$design_space$grouped_xt,rbind(c(1,2,3,4,5),c(1,2,3,4,5))) 
  expect_error(create_design_space(design_2,x_space=list(SEX=c(1,2),DOSE_discrete=seq(100,400,by=20)),grouped_xt=c(1,2,3,4,6)))
  
  expect_true(all(ds_9$design_space$grouped_xt[1,]==ds_9$design_space$grouped_xt[2,],na.rm=TRUE))
  
  expect_true(all(ds_10$design_space$grouped_a[1,]==ds_10$design_space$grouped_a[2,],na.rm=TRUE))
    
  expect_equivalent(ds_11$design_space$grouped_a, rbind(c(1,2),c(3,2)))
  
  expect_true(all(ds_12$design_space$grouped_x[1,]==ds_12$design_space$grouped_x[2,],na.rm=TRUE))
    
  expect_equivalent(ds_13$design_space$grouped_x, rbind(c(1,2),c(3,2)))
  
})

test_that("create.poped.database works for different inputs", {
  
  source("examples_fcn_doc/examples_create.poped.database.R")
  
  ## evaluate initial design
  output <- calc_ofv_and_fim(poped.db,ofv_calc_type=1)
  crit <- ofv_criterion(output$ofv,size(output$fim,1),poped.db,ofv_calc_type=1)
  
  expect_equal(crit,1794.658,tolerance=1e-3)
  expect_true(!is.null(dimnames(poped.db$design$xt)))
  
  output_2 <- calc_ofv_and_fim(poped.db,ofv_calc_type=4)
  crit_2 <- ofv_criterion(output_2$ofv,size(output_2$fim,1),poped.db,ofv_calc_type=4)
  expect_equal(crit_2,1.76066,tolerance=1e-3)
  
  poped.db_1 <- create.poped.database(ff_file="ff.PK.1.comp.oral.sd.CL",
                                      fg_file="sfg",
                                      fError_file="feps.prop",
                                      bpop=c(CL=0.15, V=8, KA=1.0, Favail=1), 
                                      notfixed_bpop=c(1,1,1,0),
                                      d=c(CL=0.07, V=0.02, KA=0.6), 
                                      sigma=0.01,
                                      groupsize=16,
                                      xt=list(c(0.5,1,2,6,24,36,72,120),
                                              c(0.5,1,2,6,24,36,72,120)),
                                      minxt=0,
                                      maxxt=120,
                                      a=c(DOSE=70))
  
  output_1 <- calc_ofv_and_fim(poped.db_1)
  
  expect_equivalent(output_1,output_2)
  expect_true(!is.null(dimnames(poped.db$design$a)))
  
})
