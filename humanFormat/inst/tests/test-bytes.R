library(testthat)
library(humanFormat)

test_that("A few small byte sizes render", {
	expect_equal("0", formatBytes(0))
	expect_equal("1", formatBytes(1))
	expect_equal("11", formatBytes(11))
})

test_that("SI sizes in a vector render", {
	expsSI <- c("2", "4", "8", "16", "32", "64", "128", "256", "512", "1.02KB", "2.05KB", 
		"4.10KB", "8.19KB", "16.38KB", "32.77KB", "65.54KB", "131.07KB", "262.14KB", "524.29KB", 
		"1.05MB", "2.10MB", "4.19MB", "8.39MB", "16.78MB", "33.55MB", "67.11MB", "134.22MB", 
		"268.44MB", "536.87MB", "1.07GB", "2.15GB", "4.29GB", "8.59GB", "17.18GB", "34.36GB", 
		"68.72GB", "137.44GB", "274.88GB", "549.76GB", "1.10TB", "2.20TB", "4.40TB", "8.80TB", 
		"17.59TB", "35.18TB", "70.37TB", "140.74TB", "281.47TB", "562.95TB", "1.13PB", "2.25PB", 
		"4.50PB", "9.01PB", "18.01PB", "36.03PB", "72.06PB", "144.12PB", "288.23PB", "576.46PB", 
		"1.15EB", "2.31EB", "4.61EB", "9.22EB", "18.45EB", "36.89EB", "73.79EB", "147.57EB", 
		"295.15EB", "590.30EB", "1.18ZB", "2.36ZB", "4.72ZB", "9.44ZB", "18.89ZB", "37.78ZB", 
		"75.56ZB", "151.12ZB", "302.23ZB", "604.46ZB", "1.21YB", "2.42YB", "4.84YB", "9.67YB", 
		"19.34YB", "38.69YB", "77.37YB", "154.74YB", "309.49YB")

	expect_equal(formatBytes(2^(1:88)), expsSI)
})

test_that("IEC sizes in a vector render", {

	expsIEC <- c("2", "4", "8", "16", "32", "64", "128", "256", "512", "1.00KiB", "2.00KiB", 
		"4.00KiB", "8.00KiB", "16.00KiB", "32.00KiB", "64.00KiB", "128.00KiB", "256.00KiB", 
		"512.00KiB", "1.00MiB", "2.00MiB", "4.00MiB", "8.00MiB", "16.00MiB", "32.00MiB", 
		"64.00MiB", "128.00MiB", "256.00MiB", "512.00MiB", "1.00GiB", "2.00GiB", "4.00GiB", 
		"8.00GiB", "16.00GiB", "32.00GiB", "64.00GiB", "128.00GiB", "256.00GiB", "512.00GiB", 
		"1.00TiB", "2.00TiB", "4.00TiB", "8.00TiB", "16.00TiB", "32.00TiB", "64.00TiB", "128.00TiB", 
		"256.00TiB", "512.00TiB", "1.00PiB", "2.00PiB", "4.00PiB", "8.00PiB", "16.00PiB", 
		"32.00PiB", "64.00PiB", "128.00PiB", "256.00PiB", "512.00PiB", "1.00EiB", "2.00EiB", 
		"4.00EiB", "8.00EiB", "16.00EiB", "32.00EiB", "64.00EiB", "128.00EiB", "256.00EiB", 
		"512.00EiB", "1.00ZiB", "2.00ZiB", "4.00ZiB", "8.00ZiB", "16.00ZiB", "32.00ZiB", 
		"64.00ZiB", "128.00ZiB", "256.00ZiB", "512.00ZiB", "1.00YiB", "2.00YiB", "4.00YiB", 
		"8.00YiB", "16.00YiB", "32.00YiB", "64.00YiB", "128.00YiB", "256.00YiB")

	expect_equal(formatIECBytes(2^(1:88)), expsIEC)
})