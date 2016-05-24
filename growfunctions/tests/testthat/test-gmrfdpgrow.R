## library(growfunctions, quietly = TRUE)

context("test gmrfdpgrow() returns correct objects")

##
## Load cps dataset
##
data(cps)

## take a portion of the data matrix for compute speed
y_short             <- cps$y[,(cps$yr_label %in% c(2011:2013))]
T                   <- ncol(y_short)
N                   <- nrow(y_short)

##
## estimation model
##
mod <- function(y, q_type, q_order, niter, nburn, nthin){
               gmrfdpgrow(y = y, q_type = q_type, q_order = q_order,
                         n.iter = niter, n.burn = nburn, n.thin = nthin)   
}

test_that("gmrfdpgrow() with no missing data and one covariance term returns correct objects", {
     q_type         <- "tr"
     q_order        <- 2
     K              <- length(q_type)
     niter		<- 6
     nburn		<- 1
     nthin		<- 1
     GMRFDP		<- mod(y_short,q_type,q_order,niter,nburn,nthin)
     
     ## evaluating class
     expect_that(GMRFDP,is_a("gmrfdpgrow"))
     ## evaluating MCMC sample results
     expect_that(ncol(GMRFDP$bb), equals(N*T))
     expect_that(nrow(GMRFDP$f), equals(length(q_type)))
     expect_that(ncol(GMRFDP$Kappa), equals(N*K))
})


test_that("gmrfdpgrow() with some missing data and one covariance term returns correct objects", {
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
     
     q_type         <- "tr"
     q_order        <- 2
     K              <- length(q_type)
     niter     	<- 6
     nburn		<- 1
     nthin		<- 1
     GMRFDP_m		<- mod(y_obs,q_type,q_order,niter,nburn,nthin)
     ## perform plots of posterior mean function values vs. data and functions grouped to cluster
     plots_gmrf     <- cluster_plot( object = GMRFDP_m,  units_name = "state", units_label = cps$st, 
                                     single_unit = FALSE, credible = TRUE )
     
     ## evaluating class
     expect_that(GMRFDP_m,is_a("gmrfdpgrow"))
     ## evaluating MCMC sample results
     expect_that(ncol(GMRFDP_m$bb), equals(N*T))
     expect_that(nrow(GMRFDP_m$f), equals(length(q_type)))
     expect_that(ncol(GMRFDP_m$Kappa), equals(N*K))
     ## check cluster_plot()
     expect_that(plots_gmrf$p.cluster, is_a("ggplot"))
     expect_that(plots_gmrf$p.fit, is_a("ggplot"))
     expect_that(nrow(plots_gmrf$map), equals(N))
})


test_that("fit_compare() that plots fit comparison of two objects returns desired plot.", {
     ## first model employs a single term
     q_type_1         <- "tr"
     q_order_1        <- 2
     ## second model employs 2 terms
     q_type_2         <- c("tr","sn")
     q_order_2        <- c(2,3)
     K                <- length(q_type_2)
     ## set run length
     niter          <- 6
     nburn		<- 1
     nthin		<- 1
     ## estimate models
     GMRFDP_1     	<- mod(y_short,q_type_1,q_order_1,niter,nburn,nthin)
     GMRFDP_2		<- mod(y_short,q_type_2,q_order_2,niter,nburn,nthin)
     ## generate clustering for comparison plot
     fit_2          <- cluster_plot( object = GMRFDP_1,  units_name = "state", units_label = cps$st, 
                                       single_unit = FALSE, credible = TRUE )
     ## generate fit comparison plot
     objects                       <- vector("list",2)
     objects[[1]]                  <- GMRFDP_1
     objects[[2]]                  <- GMRFDP_2
     label.object                  <- c("gmrf_rw2","gmrf_trsn")
     ## the existence of H is an indirect test of cluster_plot()
     H                             <- fit_2$map[order(fit_2$map$units_numeric),]$cluster
     plot_compare                  <- fit_compare( objects = objects, H = H, 
                                                   label.object = label.object,
                                                   y.axis.label = "normalized y values",
                                                   units_name = "state",
                                                   units_label = cps$st) 
     ## evaluating class
     expect_that(GMRFDP_2,is_a("gmrfdpgrow"))
     ## evaluating MCMC sample results
     expect_that(ncol(GMRFDP_2$bb), equals(N*T))
     expect_that(nrow(GMRFDP_2$f), equals(length(q_type_2)))
     ## expect_that(GMRFDP_2$f, is_a("list"))
     expect_that(ncol(GMRFDP_2$Kappa), equals(N*K))
     ## check fit_compare()
     expect_that(plot_compare$p.t, is_a("ggplot"))
     expect_that(nrow(plot_compare$map), equals(N))
})
