library(testthat)
library('gmum.r')
source('energy.R')

test_that("EllipseGauss: energy is correct", {
    data(cec_ellipse_gauss)
    dataset_points <- cec_ellipse_gauss.input
    
    t <- 5
    accepted <- 0
    nclusters <- 4
    npoints = dim(dataset_points)[1]
    for(i in 1:t)
    {
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='standard')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = standard_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)        
        
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='spherical')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = sphere_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
        
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='diagonal')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = diagonal_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
    }
    print("EllipseGauss: energy is correct")
})

test_that("mouse_1: energy is correct", {
    data(cec_mouse_1)
    dataset_points <- cec_mouse_1.input
    t <- 5
    accepted <- 0
    nclusters <- 3
    npoints <- dim(dataset_points)[1]
    for(i in 1:t)
    {
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='standard')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = standard_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
        
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='spherical')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = sphere_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
        
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='diagonal')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = diagonal_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
    }
    print("mouse_1: energy is correct")
})

test_that("mouse_1_spherical: energy is correct", {
    data(cec_mouse_1_spherical)
    dataset_points <- cec_mouse_1_spherical.input 
  
    t <- 5
    accepted <- 0
    nclusters <- 3
    npoints <- dim(dataset_points)[1]
    for(i in 1:t)
    {
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='standard')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = standard_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
        
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='spherical')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = sphere_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
        
        c <- CEC(k=nclusters, x=dataset_points, method.init='random', method.type='diagonal')
        energy_func_value <- cec_energy(dataset = dataset_points, clustering = c$clustering, entropy_func = diagonal_entropy)
        expect_equal(expected=c$energy,  object=energy_func_value, tolerance=.000001)
    }
    print("mouse_1_spherical: energy is correct")
})
