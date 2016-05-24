context("GLCM textures")

suppressMessages(library(raster))

# First run tests tests without block-by-block processing
rasterOptions(todisk=FALSE)

# Make a function to get 2d matrix from 3d matrix returned by glcm
get_pkg_glcm_texture <- function(statistic, window, shift) {
    if (length(statistic) != 1) {
        stop('length of statistic must be equal to 1')
    }
    # Note the na_val=0 is needed to match ENVI output
    texture <- glcm(test_raster, 32, window, shift, statistic, na_val=0)
    return(getValues(texture))
}

# Test all statistics that are available in EXELIS ENVI match the textures 
# output by pkg
test_that("glcm on 3x3 window with 1x1 shift works", {
    expect_equal(get_pkg_glcm_texture('mean_ENVI', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$mean_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('variance_ENVI', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$variance_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('homogeneity', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$homogeneity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('contrast', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$contrast),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('dissimilarity', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$dissimilarity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('entropy', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$entropy),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('second_moment', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$second_moment),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('correlation', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$correlation),
                 tolerance=.000001)
})

test_that("glcm on 5x7 window with 2x3 shift works", {
    expect_equal(get_pkg_glcm_texture('mean_ENVI', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$mean_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('variance_ENVI', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$variance_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('homogeneity', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$homogeneity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('contrast', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$contrast),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('dissimilarity', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$dissimilarity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('entropy', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$entropy),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('second_moment', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$second_moment),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('correlation', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$correlation),
                 tolerance=.000001)
})

## ENVI currently has a bug and does not properly handle negative shifts, so 
## don't run this test. The R glcm package properly handles these shifts.
# test_that("glcm on 5x3 window with -1,-2 shift works", {
#     expect_equal(get_pkg_glcm_texture('mean_ENVI', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$mean_ENVI),
#                  tolerance=.000001)
#     expect_equal(get_pkg_glcm_texture('variance_ENVI', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$variance_ENVI),
#                  tolerance=.000001)
#     expect_equal(get_pkg_glcm_texture('homogeneity', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$homogeneity),
#                  tolerance=.000001)
#     expect_equal(get_pkg_glcm_texture('contrast', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$contrast),
#                  tolerance=.000001)
#     expect_equal(get_pkg_glcm_texture('dissimilarity', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$dissimilarity),
#                  tolerance=.000001)
#     expect_equal(get_pkg_glcm_texture('entropy', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$entropy),
#                  tolerance=.000001)
#     expect_equal(get_pkg_glcm_texture('second_moment', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$second_moment),
#                  tolerance=.000001)
#     expect_equal(get_pkg_glcm_texture('correlation', c(5, 3), c(-1, -2)),
#                  expected=getValues(expected_textures_5x3_n1xn2$correlation),
#                  tolerance=.000001)
# })

# Test that glcm run on a raster matches the output from running glcm directly 
# on a matrix
glcm_corr_matrix <- glcm(raster::as.matrix(test_raster), 32, c(3, 3), c(1, 1), 'correlation', na_val=0)
glcm_corr_matrix <- matrix(glcm_corr_matrix, nrow=nrow(glcm_corr_matrix))
test_that("GLCM run on a matrix works correctly", {
    expect_equal(glcm_corr_matrix,
                 expected=raster::as.matrix(expected_textures_3x3_1x1$correlation),
                 tolerance=.000001)
})

glcm_corr_int <- round(glcm(test_raster, 32, c(3, 3), c(1, 1), 'correlation', na_val=0) * 1000)
test_that("GLCM scaling works correctly when run with scaling and rounding", {
    expect_equal(glcm(test_raster, 32, c(3, 3), c(1, 1), 'correlation', 
                      asinteger=TRUE, scale_factor=1000, na_val=0),
                 expected=glcm_corr_int,
                 tolerance=.000001)
})

# Re-run glcm tests with block-by-block processing
rasterOptions(todisk=TRUE)

# Test all statistics that are available in EXELIS ENVI match the textures 
# output by pkg
test_that("glcm on 3x3 window with 1x1 shift works", {
    expect_equal(get_pkg_glcm_texture('mean_ENVI', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$mean_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('variance_ENVI', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$variance_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('homogeneity', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$homogeneity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('contrast', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$contrast),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('dissimilarity', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$dissimilarity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('entropy', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$entropy),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('second_moment', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$second_moment),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('correlation', c(3, 3), c(1, 1)),
                 expected=getValues(expected_textures_3x3_1x1$correlation),
                 tolerance=.000001)
})

test_that("glcm on 5x7 window with 2x3 shift works", {
    expect_equal(get_pkg_glcm_texture('mean_ENVI', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$mean_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('variance_ENVI', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$variance_ENVI),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('homogeneity', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$homogeneity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('contrast', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$contrast),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('dissimilarity', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$dissimilarity),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('entropy', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$entropy),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('second_moment', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$second_moment),
                 tolerance=.000001)
    expect_equal(get_pkg_glcm_texture('correlation', c(5, 7), c(2, 3)),
                 expected=getValues(expected_textures_5x7_2x3$correlation),
                 tolerance=.000001)
})
