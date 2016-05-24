context("Test `qq` and `qqcat`")

test_that("Simple test", {
	a = 1
	expect_that(qq("this is @{a}"), 
				equals("this is 1"))
	expect_that(qqcat("this is @{a}"), 
				prints_text("this is 1"))
})			

test_that("pass a list as an environment", {       
	l = list(a = "a")
	expect_that(qq("this is @{a} in `l`", env = l),
				equals("this is a in `l`"))
	expect_that(qqcat("this is @{a} in `l`", env = l),
				prints_text("this is a in `l`"))
})

test_that("variables are multiple element vector", {
	a = 1:6
	expect_that(qq("@{a} is an @{ifelse(a %% 2, 'odd', 'even')} number\n"),
				equals("1 is an odd number\n2 is an even number\n3 is an odd number\n4 is an even number\n5 is an odd number\n6 is an even number\n"))
	expect_that(qq("@{a} is an @{ifelse(a %% 2, 'odd', 'even')} number\n", collapse = FALSE),
				equals(c("1 is an odd number\n", "2 is an even number\n", "3 is an odd number\n", "4 is an even number\n", "5 is an odd number\n", "6 is an even number\n")))
})

test_that("different code patterns", {
	expect_that(find_code("@\\{CODE\\}", "@{a}, @[b], @<c>, @(d), ${e}, `f`"),
				equals(list(template = "@{a}",
							code     =  "a")))
	expect_that(find_code("@\\[CODE\\]", "@{a}, @[b], @<c>, @(d), ${e}, `f`"),
				equals(list(template = "@[b]",
							code     =  "b")))
	expect_that(find_code("@<CODE>", "@{a}, @[b], @<c>, @(d), ${e}, `f`"),
				equals(list(template = "@<c>",
							code     =  "c")))       
	expect_that(find_code("@\\(CODE\\)", "@{a}, @[b], @<c>, @(d), ${e}, `f`"),
				equals(list(template = "@(d)",
							code     =  "d")))
	expect_that(find_code("\\$\\{CODE\\}", "@{a}, @[b], @<c>, @(d), ${e}, `f`"),
				equals(list(template = "${e}",
							code     =  "e")))
	expect_that(find_code("#\\{CODE\\}", "@{a}, @[b], @<c>, @(d), #{e}, `f`"),
				equals(list(template = "#{e}",
							code     =  "e"))) 
	expect_that(find_code("`CODE`", "@{a}, @[b], @<c>, @(d), #{e}, `f`"),
				equals(list(template = "`f`",
							code     =  "f")))                         
	expect_that(find_code("@\\[\\[CODE\\]\\]", "@{a}, @[b], @<c>, @(d), #{e}, `f`, @[[g]]"),
				equals(list(template = "@[[g]]",
							code     =  "g")))                         
})

test_that("simple template", {
	a = letters[1:3]
	b = 1:3
	expect_that(qq("`
	text = character(length(a))
	for(i in seq_along(a)) {
		text[i] = qq('<tr><td>@{a[i]}</td><td>@{b[i]}</td></tr>\n')
	}
	text
	`", code.pattern = "`CODE`"),
				equals("<tr><td>a</td><td>1</td></tr>\n<tr><td>b</td><td>2</td></tr>\n<tr><td>c</td><td>3</td></tr>\n"))
})

test_that("test `cat_prefix`", {
	qq.options("cat_prefix" = "INFO:")
	expect_that(qqcat("a"),
				prints_text("INFO:a"))
				
	qq.options("cat_prefix" = NULL)
	expect_that(qqcat("a"),
				prints_text("a"))
				
	qq.options("cat_prefix" = function() "DEBUG:")
	expect_that(qqcat("a"),
				prints_text("DEBUG:a"))
				
	qq.options("cat_prefix" = NULL)
	expect_that(qqcat("a"),
				prints_text("a"))

	qq.options("cat_prefix" = "INFO:", "cat_verbose" = FALSE)
	expect_silent(qqcat("a"))
				
	qq.options("cat_prefix" = "INFO:", "cat_verbose" = TRUE)
	expect_that(qqcat("a"),
				prints_text("INFO:a"))
	qq.options(RESET = TRUE)
	
	expect_that(qqcat("a", cat_prefix = "DEBUG:a"),
				prints_text("DEBUG:a"))
				
	qq.options("cat_prefix" = "INFO:")
	expect_that(qqcat("a", cat_prefix = "DEBUG:a"),
				prints_text("DEBUG:a"))
	expect_that(qqcat("a", cat_prefix = function() "DEBUG:a"),
				prints_text("DEBUG:a"))
	expect_that(qqcat("a",),
				prints_text("INFO:a"))
	qq.options(RESET = TRUE)
	
	op = qq.options(READ.ONLY = FALSE)
	qq.options(op)
	expect_that(qq.options("cat_prefix"), equals(""))
	qq.options(cat_prefix = function() "INFO:")
	op = qq.options(READ.ONLY = FALSE)
	qq.options(op)
})
