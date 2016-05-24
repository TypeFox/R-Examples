## ---- eval=FALSE, install------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("soodoku/captr")

## ----eval=FALSE, load_lib------------------------------------------------
#  library(captr)

## ---- eval=FALSE, set_token----------------------------------------------
#  # Not a real token
#  set_token("6dbee39a047c4de2b576b966")

## ---- eval=FALSE, create_batch-------------------------------------------
#  batch <- create_batch("wisc_ads")
#  batch$id

## ---- eval=FALSE, upload_image-------------------------------------------
#  path <- system.file("extdata/wisc_ads", package = "captr")
#  files <- dir(path, full.names = TRUE)
#  upimage <- lapply(files, upload_image, batch_id = batch$id)
#  
#  names(upimage[[5]])

## ---- eval=FALSE, test_readiness-----------------------------------------
#  tester <- test_readiness(batch_id=batch$id)
#  tester$errors

## ---- eval=FALSE, batch_price--------------------------------------------
#  price <- batch_price(batch_id=batch$id)
#  price$total_user_cost_in_cents

## ---- eval=FALSE, submit_batch-------------------------------------------
#  submit <- submit_batch(batch_id=batch$id)
#  submit$related_job_id

## ---- eval=FALSE, track_progress-----------------------------------------
#  progress <- track_progress(submit$related_job_id)
#  progress$percent_completed

## ---- eval=FALSE, list_instance_sets-------------------------------------
#  list_instances <- list_instance_sets(job_id=submit$related_job_id)
#  list_instances$id

## ---- eval=FALSE, get_instance_set---------------------------------------
#  res1 <- get_instance_set(instance_set_id=list_instances$id[1])
#  res1$best_estimate

## ---- eval=FALSE, get_all------------------------------------------------
#  get_all(job_id=submit$related_job_id)

