context("Examples")

test_that("The Examples run", {
  
  if(skip_examples) skip("Examples with long run times")
  
  run_dir <- file.path("..","..","inst","examples")
  source(file.path(run_dir,"ex.1.a.PK.1.comp.oral.md.intro.R"))  
  source(file.path(run_dir,"ex.1.b.PK.1.comp.oral.md.re-parameterize.R"))  
  source(file.path(run_dir,"ex.1.c.PK.1.comp.oral.md.ODE.R"))  
  source(file.path(run_dir,"ex.2.a.warfarin.evaluate.R"))  
  source(file.path(run_dir,"ex.2.b.warfarin.optimize.R"))  
  source(file.path(run_dir,"ex.2.c.warfarin.discrete.optimization.R"))  
  
  source(file.path(run_dir,"ex.2.d.warfarin.ODE.compiled.R"),chdir=T)  
  expect_equal(det(FIM),det(FIM.compiled))
  expect_equal(det(FIM)^(1/size(FIM,1)), det(FIM.1)^(1/size(FIM.1,1)),tol=1e-2)
  
  source(file.path(run_dir,"ex.2.e.warfarin.ED.R"))  
  source(file.path(run_dir,"ex.2.f.warfarin.Ds.R"))  
  source(file.path(run_dir,"ex.3.a.PKPD.1.comp.oral.md.imax.D-opt.R"))  
  source(file.path(run_dir,"ex.3.b.PKPD.1.comp.oral.md.imax.ED-opt.R"))  
  source(file.path(run_dir,"ex.4.PKPD.1.comp.emax.R"))  
  source(file.path(run_dir,"ex.5.PD.emax.hill.R"))  
  source(file.path(run_dir,"ex.6.PK.1.comp.oral.sd.R"))  
  source(file.path(run_dir,"ex.7.PK.1.comp.maturation.R"))  
  source(file.path(run_dir,"ex.8.a.tmdd_qss_one_target.R"))  
  source(file.path(run_dir,"ex.8.b.tmdd_qss_one_target_compiled.R"),chdir=T)  
  source(file.path(run_dir,"ex.9.PK.2.comp.oral.md.ode.compiled.R"),chdir=T)  
  
  source(file.path(run_dir,"ex.10.PKPD.HCV.compiled.R"),chdir=T)
  expect_equal(crit,crit_reference_reduced,tolerance=0.01,scale=crit_reference_reduced)
  for(i in 1:length(rse)) expect_equal(rse[i], rse_reference_reduced[i],tolerance=0.05,scale=rse_reference_reduced[i],check.names=F)
  expect_equal(crit_full,crit_reference_full,tolerance=0.01,scale=crit_reference_full)
  for(i in 1:length(rse_full)) expect_equal(rse_full[i],rse_reference_full[i],tolerance=0.05,scale=rse_reference_full[i],check.names=F)
  
  
})