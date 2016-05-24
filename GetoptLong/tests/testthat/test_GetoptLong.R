
context("Test `GetoptLong`")

test_that("test simple", {
	GetoptLong("tag=i", "tag", argv_str = "--tag 1");         expect_that(tag, equals(1)); rm(tag)
})

test_that("test `tag=i`", {
	spec = c(
		"tag=i", "this is a description of tag which is long long and very long and extremly long..."
	)
	GetoptLong(spec, argv_str = "--tag 1");         expect_that(tag, equals(1)); rm(tag)
	tag = 2; GetoptLong(spec, argv_str = "");         expect_that(tag, equals(2)); rm(tag)
	tag = 2; GetoptLong(spec, argv_str = "--tag 1");         expect_that(tag, equals(1)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 --tag 2"); expect_that(tag, equals(2)); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag 0.1"), prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag a"),   prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag"),     prints_text("requires"))
	expect_that(GetoptLong(spec, argv_str = ""),          prints_text("mandatory"))
})

test_that("test `len|size=i`", {
	spec = c(
		"len|size=i", "length"
	)
	GetoptLong(spec, argv_str = "--len 1");      expect_that(len, equals(1)); rm(len)
	GetoptLong(spec, argv_str = "--size 1");     expect_that(len, equals(1)); rm(len)
	GetoptLong(spec, argv_str = "-l 1");         expect_that(len, equals(1)); rm(len)
	GetoptLong(spec, argv_str = "-s 1");         expect_that(len, equals(1)); rm(len)
	expect_that(GetoptLong(spec, argv_str = ""), prints_text("mandatory"))
})

# option name is a built-in function
test_that("test `length=i`", {
	spec = c(
		"length|size=i", "length"
	)
	length = function() {}
	GetoptLong(spec, argv_str = "--length 1");   expect_that(length, equals(1)); rm(length)
	expect_that(GetoptLong(spec, argv_str = ""), prints_text("mandatory"))
})

test_that("test `tag=s`", {
	spec = c(
		"tag=s", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1");   expect_that(tag, is_identical_to("1"));   rm(tag)
	GetoptLong(spec, argv_str = "--tag 0.1"); expect_that(tag, is_identical_to("0.1")); rm(tag)
	GetoptLong(spec, argv_str = "--tag a");   expect_that(tag, is_identical_to("a"));   rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag"), prints_text("requires"))
})

test_that("test `tag=f`", {
	spec = c(
		"tag=f", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1");   expect_that(tag, equals(1));   rm(tag)
	GetoptLong(spec, argv_str = "--tag 0.1"); expect_that(tag, equals(0.1)); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag a"),  prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag"),    prints_text("requires"))
})

test_that("test `tag=o`", {
	spec = c(
		"tag=o", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1");    expect_that(tag, equals(1));  rm(tag)
	GetoptLong(spec, argv_str = "--tag 0b11"); expect_that(tag, equals(3));  rm(tag)
	GetoptLong(spec, argv_str = "--tag 011");  expect_that(tag, equals(9));  rm(tag)
	GetoptLong(spec, argv_str = "--tag 0x1f"); expect_that(tag, equals(31)); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag 0.1"), prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag a"),   prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag"),     prints_text("requires"))
})

test_that("test `tag`", {
	spec = c(
		"tag", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1");   expect_that(tag, is_identical_to(TRUE)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 0");   expect_that(tag, is_identical_to(TRUE)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 0.1"); expect_that(tag, is_identical_to(TRUE)); rm(tag)
	GetoptLong(spec, argv_str = "--tag a");   expect_that(tag, is_identical_to(TRUE)); rm(tag)
	GetoptLong(spec, argv_str = "--tag");     expect_that(tag, is_identical_to(TRUE)); rm(tag)
	GetoptLong(spec, argv_str = "");          expect_that(tag, is_identical_to(FALSE)); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--no-tag"), prints_text("Unknown"))
})

test_that("test `tag=i@`", {
	spec = c(
		"tag=i@", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1");         expect_that(tag, is_identical_to(1));   rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 --tag 2"); expect_that(tag, equals(1:2)); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag 1 --tag a"), prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag 0.1"),       prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag a"),         prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag"),           prints_text("requires"))
	expect_that(GetoptLong(spec, argv_str = ""),                prints_text("mandatory"))
})

test_that("test `tag=i%`", {
	spec = c(
		"tag=i%", "desc"
	)
	GetoptLong(spec, argv_str = "--tag name=1"); expect_that(tag, is_identical_to(list(name = 1))); rm(tag)
	GetoptLong(spec, argv_str = "--tag name=1 value=2"); expect_that(tag, is_identical_to(list(name = 1))); rm(tag)
	GetoptLong(spec, argv_str = "--tag name=1 --tag value=2"); tag = tag[sort(names(tag))]; expect_that(tag, is_identical_to(list(name = 1, value = 2))); rm(tag)
	GetoptLong(spec, argv_str = "--tag name=1 --tag name=2"); expect_that(tag, is_identical_to(list(name = 2))); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag 1"),   prints_text("requires"))
	expect_that(GetoptLong(spec, argv_str = "--tag 0.1"), prints_text("requires"))
	expect_that(GetoptLong(spec, argv_str = "--tag a"),   prints_text("requires"))
	expect_that(GetoptLong(spec, argv_str = "--tag"),     prints_text("requires"))
	expect_that(GetoptLong(spec, argv_str = ""),          prints_text("mandatory"))
	
	
	expect_that({tag = 1;GetoptLong(spec, argv_str = "--tag name=1")},   prints_text("is a list"));rm(tag)
	expect_that({tag = list(name = list(1));GetoptLong(spec, argv_str = "--tag name=1")},   prints_text("containing simple vectors")); rm(tag)
	expect_that({tag = list(1);GetoptLong(spec, argv_str = "--tag name=1")},   prints_text("with names")); rm(tag)
	tag = list(name = 2)
	GetoptLong(spec, argv_str = ""); expect_that(tag, is_identical_to(list(name = 2))); rm(tag)
	tag = list(name = 2)
	GetoptLong(spec, argv_str = "--tag name=1"); expect_that(tag, is_identical_to(list(name = 1))); rm(tag)

	tag = list(name = 2, eman = 3)
	GetoptLong(spec, argv_str = "--tag name=1"); expect_that(tag, is_identical_to(list(name = 1, eman = 3))); rm(tag)
	tag = list(name = 2, eman = 3)
	GetoptLong(spec, argv_str = "--tag name=1 --tag eman=2 --tag sth=4"); expect_that(tag, is_identical_to(list(name = 1, eman = 2, sth = 4))); rm(tag)

})

test_that("test `tag=i{2}`", {
	spec = c(
		"tag=i{2}", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1 2"); expect_that(tag, equals(1:2)); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag 1"),   prints_text("Insufficient"))
	expect_that(GetoptLong(spec, argv_str = "--tag 1 --tag 2"),   prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag 0.1"), prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag a"),   prints_text("invalid"))
	expect_that(GetoptLong(spec, argv_str = "--tag"),     prints_text("requires"))
	expect_that(GetoptLong(spec, argv_str = ""),          prints_text("mandatory"))
})

test_that("test `tag=i{2,3}`", {
	spec = c(
		"tag=i{2,3}", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1 2");   expect_that(tag, equals(1:2)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2 3"); expect_that(tag, equals(1:3)); rm(tag)
	expect_that(GetoptLong(spec, argv_str = "--tag 1"), prints_text("Insufficient"))
	GetoptLong(spec, argv_str = "--tag 1 2 3"); expect_that(tag, equals(1:3)); rm(tag)
})

test_that("test `tag=i{2,}`", {
	spec = c(
		"tag=i{2,}", "desc"
	)
	expect_that(GetoptLong(spec, argv_str = "--tag 1"), prints_text("Insufficient"))
	GetoptLong(spec, argv_str = "--tag 1 2");     expect_that(tag, equals(1:2)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2 3");   expect_that(tag, equals(1:3)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2 3 4"); expect_that(tag, equals(1:4)); rm(tag)
})

test_that("test `tag=i{,3}`", {
	spec = c(
		"tag=i{,3}", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1");       expect_that(tag, equals(1));   rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2");     expect_that(tag, equals(1:2)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2 3");   expect_that(tag, equals(1:3)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2 3 4"); expect_that(tag, equals(1:3)); rm(tag)
})

test_that("test `tag=i{,}`", {
	spec = c(
		"tag=i{,}", "desc"
	)
	GetoptLong(spec, argv_str = "--tag 1");       expect_that(tag, equals(1));   rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2");     expect_that(tag, equals(1:2)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2 3");   expect_that(tag, equals(1:3)); rm(tag)
	GetoptLong(spec, argv_str = "--tag 1 2 3 4"); expect_that(tag, equals(1:4)); rm(tag)
})

test_that("test `verbose!`", {
	spec = c(
		"verbose!", "print messages"
	)
	GetoptLong(spec, argv_str = "--verbose");    expect_that(verbose, is_identical_to(TRUE)); rm(verbose)
	GetoptLong(spec, argv_str = "--no-verbose"); expect_that(verbose, is_identical_to(FALSE)); rm(verbose)
	GetoptLong(spec, argv_str = "--noverbose"); expect_that(verbose, is_identical_to(FALSE)); rm(verbose)
})

test_that("test other configurations", {
	GetoptLong.options("config" = "bundling")
	spec = c(
		"red|r", "using red",
		"blue|b", "using blue",
		"yellow|y", "using yellow"
	)
	GetoptLong(spec, argv_str = "--red --blue --yellow")
	expect_that(red,    is_identical_to(TRUE)); rm(red)
	expect_that(blue,   is_identical_to(TRUE)); rm(blue)
	expect_that(yellow, is_identical_to(TRUE)); rm(yellow)
	GetoptLong(spec, argv_str = "-rby")
	expect_that(red,    is_identical_to(TRUE)); rm(red)
	expect_that(blue,   is_identical_to(TRUE)); rm(blue)
	expect_that(yellow, is_identical_to(TRUE)); rm(yellow)
})

# with default value
test_that("test default values", {
	spec = c(
		"tag=i", "this is a description of tag which is long long and very long and extremly long..."
	)
	tag = 4
	GetoptLong(spec, argv_str = "--tag 1");  expect_that(tag, equals(1)); rm(tag)
	tag = 4
	GetoptLong(spec, argv_str = "");         expect_that(tag, equals(4)); rm(tag)
})

# message

GetoptLong.options(RESET = TRUE)

test_that("test `version` and `help` options", {

	spec = c(
		"verbose!", "print messages"
	)
	VERSION = "0.0.1"
	expect_that(GetoptLong(spec, argv_str = "--version"), prints_text("0.0.1"))
	rm(VERSION)
	expect_that(GetoptLong(spec, argv_str = "--version"), prints_text("No version information is found in source code."))
	expect_that(GetoptLong(spec, argv_str = "--version", version = FALSE), prints_text("Unknown"))
	
	expect_that(GetoptLong(spec, argv_str = "--help"), prints_text("Usage"))

	GetoptLong.options("startingMsg" = "
Description of this script
")

	GetoptLong.options("endingMsg" = "
Report bugs to xxx@xx.xx
")
	expect_that(GetoptLong(spec, argv_str = "--help"), prints_text("Report bugs"))
	expect_that(GetoptLong(spec, argv_str = "--help", help = FALSE), prints_text("Unknown"))
})

perl_bin = Sys.which("perl")
test_that(qq("--tag 1 -- @{perl_bin}"), {
	spec = c(
		"tag=i", "this is a description of tag which is long long and very long and extremly long..."
	)
	GetoptLong(spec, argv_str = qq("--tag 1 -- @{perl_bin}")); expect_that(tag, equals(1)); rm(tag)
})

test_that("no default value under interactive session", {
	spec = c(
		"tag=i", "tag"
	)
	expect_error(GetoptLong(spec))
})
