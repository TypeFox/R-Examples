library(hashids)
context("Encoding Unit Tests")

test_that("Test creating settings list", {
	expect_error(hashid_settings(salt = '', alphabet='abcabc'),
		'hashid_settings: Alphabet must be at least 16 unique characters.')
})

test_that("Test empty encode", {
	h = hashid_settings(salt = '')
	expect_error(encode('', h), 
		'encode: numbers must be non-negative')
})

test_that("No salt encode, single digit encode", {
	h = hashid_settings(salt = '')
	expect_equal(encode(12345, h), 'j0gW')
	expect_equal(encode(1, h), 'jR')
	expect_equal(encode(22, h), 'Lw')
	expect_equal(encode(333, h), 'Z0E')
	expect_equal(encode(9999, h), 'w0rR')
})

test_that("No salt encode, encode vectors", {
	h = hashid_settings(salt = '')
	expect_equal(encode(c(1,2,3), h), 'o2fXhV')
	expect_equal(encode(c(2,4,6), h), 'xGhmsW')
	expect_equal(encode(c(99, 25), h), '3lKfD')
	expect_equal(encode(c(683, 94108, 123, 5), h), 'vJvi7On9cXGtD')
})

test_that("encode with salt", {
	h = hashid_settings(salt = 'Arbitrary string')
	expect_equal(encode(c(683, 94108, 123, 5), h), 'QWyf8yboH7KT2')
	expect_equal(encode(c(1, 2, 3), h), 'neHrCa')
	expect_equal(encode(c(2, 4, 6), h), 'LRCgf2')
	expect_equal(encode(c(99, 25), h), 'JOMh1')
	expect_equal(encode(1337, h), 'VOd')
})

test_that("encode sequential and identical vectors", {
	h = hashid_settings(salt = 'this is my salt')
	expect_equal(encode(c(1, 2, 3, 4), h), 'agHLu9hm')
	expect_equal(encode(c(4, 3, 2, 1), h), 'aQtEhBub')
	expect_equal(encode(c(1, 1, 1, 1), h), '2bHEH5HY')
})

test_that("encode with no salt, custom alphabet", {
	h = hashid_settings(
		salt = '',
		alphabet='!"#%&\',-/0123456789:;<=>ABCDEFGHIJKLMNOPQRSTUVWXYZ_`abcdefghijklmnopqrstuvwxyz~'
	)
	expect_equal(encode(c(2839, 12, 32, 5), h), '_nJUNTVU3')
	expect_equal(encode(c(1, 2, 3), h), '7xfYh2')
	expect_equal(encode(23832, h), 'Z6R>')
	expect_equal(encode(c(99, 25), h), 'AYyIB')
})

test_that("encode with salt, custom alphabet", {
	h = hashid_settings(
		alphabet='!"#%&\',-/0123456789:;<=>ABCDEFGHIJKLMNOPQRSTUVWXYZ_`abcdefghijklmnopqrstuvwxyz~',
		salt = 'Arbitrary string'
	)
	expect_equal(encode(c(683, 94108, 123, 5), h), "6'nHG/qf>oC8")
	expect_equal(encode(c(1, 2, 3), h), 'GEHMCe')
	expect_equal(encode(c(2, 4, 6), h), 'YlCrfQ')
	expect_equal(encode(23832, h), 'jwe5')
	expect_equal(encode(c(99, 25), h), 'Jz3fR')
})

test_that("encode with set min length, no salt", {
	h = hashid_settings(min_length = 25, salt = '')
	expect_equal(encode(c(7452, 2967, 21401), h), 'pO3K69b86jzc6krI416enr2B5')
	expect_equal(encode(c(1, 2, 3), h), 'gyOwl4B97bo2fXhVaDR0Znjrq')
	expect_equal(encode(6097, h), 'Nz7x3VXyMYerRmWeOBQn6LlRG')
	expect_equal(encode(c(99, 25), h), 'k91nqP3RBe3lKfDaLJrvy8XjV')
})

test_that("encode with all parameters non-default", {
	h = hashid_settings(
		salt = 'arbitrary salt',
		min_length = 16,
		alphabet = 'abcdefghijklmnopqrstuvwxyz'
	)
	expect_equal(encode(c(7452, 2967, 21401), h), 'wygqxeunkatjgkrw')
	expect_equal(encode(c(1, 2, 3), h), 'pnovxlaxuriowydb')
	expect_equal(encode(60125, h), 'jkbgxljrjxmlaonp')
	expect_equal(encode(c(99, 25), h), 'erdjpwrgouoxlvbx')
})

test_that("encode without standard separators", {
	h = hashid_settings(
		salt = '',
		alphabet='abdegjklmnopqrvwxyzABDEGJKLMNOPQRVWXYZ1234567890'
	)
	expect_equal(encode(c(7452, 2967, 21401), h), 'X50Yg6VPoAO4')
	expect_equal(encode(c(1, 2, 3), h), 'GAbDdR')
	expect_equal(encode(60125, h), '5NMPD')
	expect_equal(encode(c(99, 25), h), 'yGya5')
})

test_that("encode with two standard separators", {
	h = hashid_settings(
		salt = '',
		alphabet='abdegjklmnopqrvwxyzABDEGJKLMNOPQRVWXYZ1234567890uC'
	)
	expect_equal(encode(c(7452, 2967, 21401), h), 'GJNNmKYzbPBw')
	expect_equal(encode(c(1, 2, 3), h), 'DQCXa4')
	expect_equal(encode(60125, h), '38V1D')
	expect_equal(encode(c(99, 25), h), '373az')
})

test_that("encode negative errors and float errors", {
	h = hashid_settings(salt = '')
	expect_error(encode(c(1, -2, 3), h),
	'encode: numbers must be non-negative')
	expect_error(encode(c(1, 2.5, 3), h),
	'encode: numbers must be integers')
})

test_that("test encode_hex", {
	h = hashid_settings(salt = '')
	h2 = hashid_settings(salt = '', min_length = 1000)
	expect_equal(encode_hex('507f1f77bcf86cd799439011', h), 'y42LW46J9luq3Xq9XMly')
	expect_equal(encode_hex('507f1f77bcf86cd799439011', h2), 
		'q0BxAAGOo0P13jn9Nw4WVkx5LlrLk5r398KOPj7JWvyDEgQnBGJ6QEg0kqnKYv71MpWZPw9YQ9KgREqwV51XzLpoA3Byj89PR2xo0yELBnD1vl7zqj5kN3lvYpmzOBVqE9QJLWZ1yvXnRJK20VZ5pAQ1N4YwLPkgNj7ApZ9v82QBRDo3wz4K0O8v0gwoPVj3NW9GpkDKBqRljxYAO6pgRk7EMzGvrW40Vkwl9xn70X4vBGZ1KWA3oEroMrAKX1gZLjp62EBwqGNyvYQLOnv6q9gzBwjNGmP712p5M7ABVE14W2p63kxvjlG9Z8KwpGP5yqjL9m0J7MEQg1zk8KAqnNGV604zEQw9lY5rxZ5nlRO6yQpvJYDBoW1LA8VnERgLq1V8yN0rXAmGW5vZDAy4582E3QDljzpVkY1Kg0nBPYjwmNVp4RQq3xWAEMkzDP5ZD9jnQkr7201yEw3mWvVyqgMDwokNLP0x6213EGOjZdy42LW46J9luq3Xq9XMlybJ4mz5nA98XlvrKQRpYWVB7xpgKAB8LMz4R6OJlNqYoXGKG7v8X1O06ZJgol2rL59ynPmO9rZBXNLJWoRwqG7Mx6vOPz3w4Y29JKjp7xB6MkQloE73Xq0kjPxM4rgzwK2NmG9gR2OpLjoZJyWvMm317DBXPDlXWoNkZ6r2RB3nvOV4AxYQP8RzOLrYywoKqmD0JnNXgyAXxEW5KMJ803o4kDrVZlR75x0WzR34YJVQkl9m8PODnVLMN6qDYOmjRP5Qygz8J2pPK3L958oyQDmnBJXNwZ21qMA7l6Q254nJxrEXLz1ymYZq5XGJnWVmYPkOrlyLE61MxGBEg6rmD79zojyOWq8xl3MXGP65AwxRn4rg78M0D2oKjkWgwAKQZOMrGVNY4JX3m6p72MJ6l0Drv4nNZOPxWmG8klBDoXRrm8A5xzON32yVjL4M40zGZXRw2qpNVlmx1Y6oA678YpRgmQyMEDZqvzJBKX2jw354'
	)
	expect_equal(encode_hex('f000000000000000000000000000000000000000000000000000000000000000000000000000000000000f', h), 
    	'WxMLpERDrmh25Lp4L3xEfM6WovWYO3IjkRMKR2ogCMVzn4zQlqt1WK8jKq7OsEpy2qyw1Vi2p'
	)
})

test_that("test invalid encode_hex", {
	h = hashid_settings(salt = '')
	expect_error(encode_hex('', h), 
	'encode_hex: invalid hex')
	expect_error(encode_hex('1234SGT8', h),
	'Invalid hex character')
})
