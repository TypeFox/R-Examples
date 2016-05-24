library(hashids)
context("Decoding Unit Tests")

test_that("Test decoding empty string", {
	h = hashid_settings(salt = '')
	expect_error(decode('', h),
	'decode: invalid hashid')
})

test_that("Decode with empty salt, single numbers", {
	h = hashid_settings(salt = '')
	expect_equal(decode('j0gW', h), 12345)
	expect_equal(decode('jR', h), 1)
	expect_equal(decode('Lw', h), 22)
	expect_equal(decode('Z0E', h), 333)
	expect_equal(decode('w0rR', h), 9999)
})

test_that("Decode with empty salt, vectors of numbers", {
	h = hashid_settings(salt = '')
	expect_equal(decode('vJvi7On9cXGtD', h), c(683, 94108, 123, 5))
	expect_equal(decode('o2fXhV', h), c(1, 2, 3))
	expect_equal(decode('xGhmsW', h), c(2, 4, 6))
	expect_equal(decode('3lKfD', h), c(99, 25))	
})

test_that("Decode with arbitrary string", {
	h = hashid_settings(salt = 'Arbitrary string')
	expect_equal(decode('QWyf8yboH7KT2', h), c(683, 94108, 123, 5))
	expect_equal(decode('neHrCa', h), c(1, 2, 3))
	expect_equal(decode('LRCgf2', h), c(2, 4, 6))
	expect_equal(decode('JOMh1', h), c(99, 25))	
	expect_equal(decode('l16D', h), 12345)
})

test_that("Decode with alphabet", {
	h = hashid_settings(
		salt = '',
		alphabet='!"#%&\',-/0123456789:;<=>ABCDEFGHIJKLMNOPQRSTUVWXYZ_`abcdefghijklmnopqrstuvwxyz~'
	)
	expect_equal(decode('_nJUNTVU3', h), c(2839, 12, 32, 5))
	expect_equal(decode('7xfYh2', h), c(1, 2, 3))
	expect_equal(decode('Z6R>', h), 23832)
	expect_equal(decode('AYyIB', h), c(99, 25))
})

test_that("Decode with different min_length", {
	h = hashid_settings(salt = '', min_length = 25)
	expect_equal(decode('pO3K69b86jzc6krI416enr2B5', h), c(7452, 2967, 21401))
	expect_equal(decode('gyOwl4B97bo2fXhVaDR0Znjrq', h), c(1, 2, 3))
	expect_equal(decode('Nz7x3VXyMYerRmWeOBQn6LlRG', h), 6097)
	expect_equal(decode('k91nqP3RBe3lKfDaLJrvy8XjV', h), c(99, 25))
})

test_that("decode with salt, custom alphabet", {
	h = hashid_settings(
		alphabet='!"#%&\',-/0123456789:;<=>ABCDEFGHIJKLMNOPQRSTUVWXYZ_`abcdefghijklmnopqrstuvwxyz~',
		salt = 'Arbitrary string'
	)
	expect_equal(decode("6'nHG/qf>oC8", h), c(683, 94108, 123, 5))
	expect_equal(decode('GEHMCe', h), c(1, 2, 3))
	expect_equal(decode('YlCrfQ', h), c(2, 4, 6))
	expect_equal(decode('jwe5', h), 23832)
	expect_equal(decode('Jz3fR', h), c(99, 25))
})

test_that("decode sequential and identical vectors", {
	h = hashid_settings(salt = 'this is my salt')
	expect_equal(decode('agHLu9hm', h), c(1, 2, 3, 4))
	expect_equal(decode('aQtEhBub', h), c(4, 3, 2, 1))
	expect_equal(decode('2bHEH5HY', h), c(1, 1, 1, 1))
})

test_that("Decode with all non-default parameters", {
	h = hashid_settings(salt = 'arbitrary salt', 
		min_length = 16,
		alphabet = 'abcdefghijklmnopqrstuvwxyz'
	)
	expect_equal(decode('wygqxeunkatjgkrw', h), c(7452, 2967, 21401))
	expect_equal(decode('pnovxlaxuriowydb', h), c(1, 2, 3))
	expect_equal(decode('jkbgxljrjxmlaonp', h), 60125)
	expect_equal(decode('erdjpwrgouoxlvbx', h), c(99, 25))
})

test_that("Test decoding invalid hashid", {
	h = hashid_settings(salt = '',
		alphabet = 'abcdefghijklmnop'
	)
	expect_error(decode('qrstuvwxyz', h),
	'decode: invalid hashid, cannot decode')
})

test_that("Test decode without standard separators", {
	h = hashid_settings(salt = '',
		alphabet='abdegjklmnopqrvwxyzABDEGJKLMNOPQRVWXYZ1234567890'
	)
	expect_equal(decode('X50Yg6VPoAO4', h), c(7452, 2967, 21401))
	expect_equal(decode('GAbDdR', h), c(1, 2, 3))
	expect_equal(decode('5NMPD', h), 60125)
	expect_equal(decode('yGya5', h), c(99, 25))
})

test_that("Test decode with two standard separators", {
	h = hashid_settings(salt = '',
		alphabet='abdegjklmnopqrvwxyzABDEGJKLMNOPQRVWXYZ1234567890uC'
	)
	expect_equal(decode('GJNNmKYzbPBw', h), c(7452, 2967, 21401))
	expect_equal(decode('DQCXa4', h), c(1, 2, 3))
	expect_equal(decode('38V1D', h), 60125)
	expect_equal(decode('373az', h), c(99, 25))	
})

test_that("Invalid hash id", {
	h = hashid_settings(salt = '', min_length = 6)
	hashed = encode(1, h)
	invalid_hash = substr(hashed, 2, nchar(hashed))
	expect_error(
		decode(invalid_hash, h),
		"decode: invalid hashid, cannot decode")
})

test_that("Test decode_hex", {
	h = hashid_settings(salt = '')
	h1 = hashid_settings(salt = '', min_length = 1000)
	expect_equal(decode_hex('y42LW46J9luq3Xq9XMly', h), '507f1f77bcf86cd799439011')
	expect_equal(
		decode_hex(
			encode_hex('507f1f77bcf86cd799439011', h), 
			h
		), 
		'507f1f77bcf86cd799439011'
	)
	expect_equal(
		decode_hex('WxMLpERDrmh25Lp4L3xEfM6WovWYO3IjkRMKR2ogCMVzn4zQlqt1WK8jKq7OsEpy2qyw1Vi2p', h),
		'f000000000000000000000000000000000000000000000000000000000000000000000000000000000000f'
	)
})

test_that("test invalid decode_hex", {
	h = hashid_settings(salt = '')
	expect_error(decode_hex('', h), 'decode: invalid hashid')
	expect_error(
		decode_hex('WxMLpERDrmh25Lp4L3xEfM6WovWYO3IjkRMKR2ogCMVlqt1WK8jKq7OsEp1Vi2p', h),
		'decode: invalid hashid, cannot decode'
	)
})