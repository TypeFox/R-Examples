## library(growfunctions, quietly = TRUE)

context("test gpdpgrow() returns correct objects")

##
## Load cps dataset
##
data(cps)

## take a portion of the data matrix for compute speed
y_short             <- cps$y[,(cps$yr_label %in% c(2012:2013))]
T                   <- ncol(y_short)
N                   <- nrow(y_short)

##
## estimation model
##
mod <- function(y, gp_cov, sn_order=NULL, niter, nburn, nthin){
          gpdpgrow(y = y, jitter = 0.001, gp_cov = gp_cov, 
                   sn_order = sn_order, n.iter = niter, n.burn = nburn, 
                   n.thin = nthin, progress = FALSE, n.tune = 0)   
	 }

test_that("gpdpgrow() with no missing data and one covariance term returns correct objects", {
	P              <- 3 ## number of covariance parameters associated with a rational quadratic
	niter		<- 5
	nburn		<- 1
	nthin		<- 1
	GPDP		     <- mod(y = y_short, gp_cov = "rq",niter = niter,nburn = nburn, nthin =nthin)
     ## perform plots of posterior mean function values vs. data and functions grouped to cluster
	plots_gp       <- cluster_plot( object = GPDP,  units_name = "state", units_label = cps$st, 
	                                    single_unit = FALSE, credible = TRUE )

	## evaluating class
	expect_that(GPDP,is_a("gpdpgrow"))
	## evaluating MCMC sample results
	expect_that(ncol(GPDP$bb), equals(N*T))
	## expect_that(GPDP$f, is_a("list"))
	expect_that(nrow(GPDP$f), equals(1))
	expect_that(ncol(GPDP$Theta), equals(N*P))
     ## check cluster_plot()
	expect_that(plots_gp$p.cluster, is_a("ggplot"))
	expect_that(plots_gp$p.fit, is_a("ggplot"))
	expect_that(nrow(plots_gp$map), equals(N))
})


test_that("gpdpgrow() with some missing data and one covariance term returns correct objects", {
     ## insert missing values in observed data matrix, y
     # randomly assign missing positions in y.
     # assume every unit has equal number of missing positions
     # randomly select number of missing observations for each unit
     m_factor  = .1
     M         = floor(m_factor*N*T)
     m_vec     = rep(floor(M/N),N)
     if( sum(m_vec) < M )
     {
          m_left              <- M - sum(m_vec)
          pos_i               <- sample(1:N, m_left, replace = FALSE)
          m_vec[pos_i]        <- m_vec[pos_i] + 1
     } # end conditional statement on whether all missing cells allocated
     ## randomly select missing positions for each unit
     pos       <- matrix(0,N,T)
     for( i in 1:N )
     {
          sel_ij              <- sample(3:(T-3), m_vec[i], replace = FALSE) ## avoid edge effects
          pos[i,sel_ij]       <- 1
     }
     
     ## blank cells in response corresponding to missing positions
     y_obs               <- y_short
     y_obs[pos == 1]     <- NA       
     
     ## set number of iterations
     P              <- 3 ## number of covariance parameters associated with a rational quadratic
     niter		<- 5
     nburn		<- 1
     nthin		<- 1
     GPDP_m		<- mod(y = y_obs,gp_cov = "rq",niter = niter,nburn = nburn,nthin = nthin)
     
     ## evaluating class
     expect_that(GPDP_m,is_a("gpdpgrow"))
     ## evaluating MCMC sample results
     expect_that(ncol(GPDP_m$bb), equals(N*T))
     expect_that(ncol(GPDP_m$Theta), equals(N*P))
})


test_that("gpdpgrow() with no missing data and two covariance term returns correct objects", {
     P              <- 6 ## number of covariance parameters associated with a rational quadratic
     niter		<- 5
     nburn		<- 1
     nthin		<- 1
     gp_cov_2       <- c("rq","sn")
     GPDP_2		<- mod(y_short,gp_cov_2,3,niter,nburn,nthin)
     
     ## evaluating class
     expect_that(GPDP_2,is_a("gpdpgrow"))
     ## evaluating MCMC sample results
     expect_that(ncol(GPDP_2$bb), equals(N*T))
     expect_that(nrow(GPDP_2$f), equals(length(gp_cov_2)))
     expect_that(ncol(GPDP_2$Theta), equals(N*P))
})



