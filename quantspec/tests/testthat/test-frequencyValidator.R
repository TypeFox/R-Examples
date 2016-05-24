context("frequencyValidator")

test_that("frequencyValidator works as expected",{
			
			freq <- 2*pi*c(3,2,11,8,5)/10
			
			res <- frequenciesValidator(freq, N=10, steps=1:3)*10/(2*pi)
			expect_that(res, equals(c(3,2,1,8,5)))
			
			res <- frequenciesValidator(freq, N=10, steps=1:4)*10/(2*pi)
			expect_that(res, equals(c(3,2,1,2,5)))
			
			res <- frequenciesValidator(freq, N=10, steps=1:5)*10/(2*pi)
			expect_that(res, equals(c(3,2,1,5)))
			
			res <- frequenciesValidator(freq, N=10, steps=1:6)*10/(2*pi)
			expect_that(res, equals(c(1,2,3,5)))
			
			freq <- freq + 0.1 # No Fourier freq. any more!
			res <- (frequenciesValidator(freq, N=10, steps=6)-0.1)*10/(2*pi)
			expect_that(res, equals(c(2,3,5,8,11)))
			
			expect_that(res <- frequenciesValidator(freq, N=10, steps=3)*10/(2*pi),
					gives_warning())
			expect_that(res, equals(c(3,2,11,8,5)))
		}
)

