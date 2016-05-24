context("FIM calculation")

test_that("RSE from evaluate.fim", {
  
  source("examples_fcn_doc/examples_evaluate.fim.R")
  
  expected.reduced <- c(4.7,2.8,13.9,25.6,30.3,25.8,11.2)
  expected.full <- c(3.6, 2.6,4.8,26.3, 30.9, 26.5, 12.4)
  comp.red.1 <- get_rse(FIM.1,poped.db)
  comp.red.4 <- get_rse(FIM.4,poped.db,fim.calc.type=4)
  comp.full.0 <- get_rse(FIM.0,poped.db)
  
  
  for(i in 1:length(expected.reduced)){
    expect_that(round(comp.red.1[[i]],digits=1), equals(expected.reduced[i], 
                                        tolerance = 0.01, scale = expected.reduced[i]))
  }
  for(i in 1:length(expected.full)){
    expect_that(round(comp.full.0[[i]],digits=1), equals(expected.full[i], 
                                         tolerance = 0.01, scale = expected.full[i]))
  }
  for(i in 1:(length(expected.reduced)-1)){
    expect_that(round(comp.red.4[[i]],digits=1), equals(expected.reduced[i], 
                                        tolerance = 0.01, scale = expected.reduced[i]))
  }
  
})

test_that("det(FIM) using evaluate.fim, approx and derivative types", {
  
  source("examples_fcn_doc/warfarin_basic.R")
  
  ## check all FIM calculations
  df <- c()
  for(i in c(0,1,4,5,6,7)){
    for(j in c(0,1)){
      FIM <- evaluate.fim(poped.db,fim.calc.type=i,deriv.type=j) 
      tmp <- data.frame("fim.calc.type"= i, "deriv.type"=j, "det.FIM"=det(FIM))
      df <- rbind(df,tmp)
    }
  }
  #print(df,digits=3,row.names=F,print.gap=3)

  for(i in c(0,1,4,5,6,7)){
    expect_that(df[df$fim.calc.type==i & df$deriv.type==0,]$det.FIM, 
                equals(df[df$fim.calc.type==i & df$deriv.type==1,]$det.FIM,
                       tolerance=1e-5))
  }
  
  expect_that(df[df$fim.calc.type==0 & df$deriv.type==0,]$det.FIM, 
              equals(df[df$fim.calc.type==5 & df$deriv.type==0,]$det.FIM,
                     tolerance=1e-5))
  expect_that(df[df$fim.calc.type==0 & df$deriv.type==0,]$det.FIM, 
              equals(df[df$fim.calc.type==6 & df$deriv.type==0,]$det.FIM,
                     tolerance=1e-5))
  expect_that(df[df$fim.calc.type==1 & df$deriv.type==0,]$det.FIM, 
              equals(df[df$fim.calc.type==7 & df$deriv.type==0,]$det.FIM,
                     tolerance=1e-5))
  
})

test_that("ofv calculation", {
  
  source("examples_fcn_doc/warfarin_optimize.R")
  
  FIM <- evaluate.fim(poped.db) # new name for function needed
  
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=1), equals(det(FIM)))
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=2), equals(1/trace_matrix(inv(FIM))))
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=4), equals(log(det(FIM)))) 
  expect_that(ofv_fim(FIM,poped.db,ofv_calc_type=7), equals(1/sum(get_rse(FIM,poped.db,use_percent=FALSE))))
})

  
test_that("internal FIM calculations", {
  
  source("examples_fcn_doc/warfarin_optimize.R")
  
  source("examples_fcn_doc/examples_mf.R")    
  expect_that(det(output$ret*32), is_identical_to(det(evaluate.fim(poped.db,fim.calc.type=0))))
  source("examples_fcn_doc/examples_mf3.R")
  expect_that(det(output$ret*32), is_identical_to(det(evaluate.fim(poped.db,fim.calc.type=1))))
  source("examples_fcn_doc/examples_mf5.R")
  expect_that(det(output$ret*32), is_identical_to(det(evaluate.fim(poped.db,fim.calc.type=4))))
  source("examples_fcn_doc/examples_mf6.R")
  expect_that(det(output$ret*32), is_identical_to(det(evaluate.fim(poped.db,fim.calc.type=5))))
  source("examples_fcn_doc/examples_mf7.R")
  expect_that(det(output$ret*32), is_identical_to(det(evaluate.fim(poped.db,fim.calc.type=6))))
  source("examples_fcn_doc/examples_mf8.R")
  expect_that(det(output$ret*32), is_identical_to(det(evaluate.fim(poped.db,fim.calc.type=7))))
  
})

