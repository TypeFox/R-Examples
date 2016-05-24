require(SensMixed)
checkMAM <- TRUE

if(checkMAM){
  
  #convert some variables to factors in TVbo
  TVbo <- convertToFactors(TVbo, c("Assessor", "Repeat", "Picture"))
  ##
  
  result <- sensmixed(names(TVbo)[5:(ncol(TVbo))],
                      Prod_effects=c("TVset", "Picture"),
                      replication="Repeat", individual="Assessor", data=TVbo, 
                      calc_post_hoc = TRUE)
  
  result
  result$fixed
  
  result_MAM <- sensmixed(names(TVbo)[5:(ncol(TVbo))],
                          Prod_effects=c("TVset", "Picture"),
                          replication="Repeat", individual="Assessor", 
                          data=TVbo,
                          MAM = TRUE)
  
  result_MAM
  
  result_MAM_mult <- sensmixed(names(TVbo)[5:(ncol(TVbo))],
                               Prod_effects=c("TVset", "Picture"),
                               replication="Repeat", individual="Assessor", 
                               data=TVbo,
                               MAM = TRUE, mult.scaling = TRUE)
  
  result_MAM_mult
  
  atts <- names(TVbo)[5:10]
  res_MAM <- sensmixed(atts, Prod_effects = c("TVset"), 
                       individual="Assessor", data=TVbo, 
                       MAM_PER=TRUE)
  
  ## selecting only 3 attributes
  atts <- names(TVbo)[5:10]
  res <- sensmixed(atts, Prod_effects=c("TVset"), individual="Assessor", 
                   data=TVbo, MAM=TRUE, reduce.random=FALSE, 
                   calc_post_hoc = TRUE)
  
  res_MAM <- sensmixed(atts, Prod_effects=c("TVset"), individual="Assessor", 
                       data=TVbo, MAM_PER=TRUE)
  
  #####################################################
  ## bug for Noise attribute - different with MAManalysis
  TVbo$x <- rep(NA, nrow(TVbo))
  x.prd <- scale(predict(lm(Noise ~ TVset, data=TVbo)), scale=FALSE)
  notNA <- rownames(x.prd)
  TVbo[notNA, "x"] <- x.prd
  
  m.noisemam <- lmer(Noise ~ TVset + x:Assessor +
                            (1|Assessor) + (1|TVset:Assessor), data=TVbo )
  all.equal(anova(m.noisemam, ddf="lme4")[,"F value"], 
            anova(m.noisemam, type = 1)[, "F.value"])
  #######################################################
  
  TOL <- 1e-2
  for(i in 1:length(atts)){
    ## ERROR: different results for the Noise attribute
    ## because of Assessor:TVset std dev = 0
    if(i==3)
      tools::assertError(stopifnot(all.equal(res_MAM[[3]][, , i][2:3,"F"], 
                                             c(res$fixed$Fval[, i],
                                               res$scaling$FScaling[, i]), 
                                             tol=TOL, check.attributes = FALSE)))
    else
      stopifnot(all.equal(res_MAM[[3]][, , i][2:3,"F"], c(res$fixed$Fval[, i],
                                                          res$scaling$FScaling[, i]), 
                          tol=TOL, check.attributes = FALSE))
  }
  
  
  TOL <- 1e-2
  ## diff lsmeans agree
  for(i in 1:length(atts)){
    stopifnot(all.equal(res_MAM[[4]][, , i][1, 1], res$post_hoc[[i]][1, 1], tol=TOL, 
                        check.attributes = FALSE))
    stopifnot(all.equal(res_MAM[[4]][, , i][1, 2], res$post_hoc[[i]][2, 1], tol=TOL, 
                        check.attributes = FALSE))
    stopifnot(all.equal(res_MAM[[4]][, , i][3, 2], res$post_hoc[[i]][3, 1], tol=TOL, 
                        check.attributes = FALSE))
  }
  
  ## selecting another attributes to compare
  atts <- names(TVbo)[10:12]
  res <- sensmixed(atts, Prod_effects=c("TVset"), individual="Assessor", 
                   data=TVbo, MAM=TRUE, reduce.random=FALSE, 
                   calc_post_hoc = TRUE)
  
  res_MAM <- sensmixed(atts, Prod_effects=c("TVset"), individual="Assessor", 
                       data=TVbo, MAM_PER=TRUE)
  TOL <- 1e-1
  for(i in 1:length(atts)){
    print(i)
    stopifnot(all.equal(res_MAM[[3]][, , i][2:3,"F"], c(res$fixed$Fval[, i],
                                                        res$scaling$FScaling[, i]), 
                        tol=TOL, check.attributes = FALSE))
  }
  
  TOL <- 1e-1
  for(i in 1:length(atts)){
    stopifnot(all.equal(res_MAM[[4]][, , i][1, 1], res$post_hoc[[i]][1, 1], 
                        tol=TOL, 
                        check.attributes = FALSE))
    stopifnot(all.equal(res_MAM[[4]][, , i][1, 2], res$post_hoc[[i]][2, 1], 
                        tol=TOL, 
                        check.attributes = FALSE))
    stopifnot(all.equal(res_MAM[[4]][, , i][3, 2], res$post_hoc[[i]][3, 1], 
                        tol=TOL, 
                        check.attributes = FALSE))
  }
  
  ## check post-hoc for the multiway case
#   atts <- names(TVbo)[5:19]
#   res <- sensmixed(atts, Prod_effects=c("TVset", "Picture"), individual="Assessor",                  
#                    data=TVbo, MAM=TRUE, calc_post_hoc = TRUE)
  
}
