library(testthat)

test_that('Testing "returns" case using with_mock on a function in external package', {
  stub_builder <- stub(sub)
  stub_builder$returns('hang on')
  sub_stub <- stub_builder$f

  with_mock(sub = sub_stub,
            expect_equal(tools::file_path_sans_ext('test5'), 'hang on'),
            .env = 'tools')
})

test_that('Testing "throws" case using with_mock on a function in external package', {
  stub_builder <- stub(sub)
  stub_builder$throws('kkkkkkkkkk')
  sub_stub <- stub_builder$f

  with_mock(sub = sub_stub,
            expect_error(tools::file_path_sans_ext('dsfsd'), 'kkkkkkkkkk'),
            expect_error(tools::file_path_sans_ext('gfgxgx'), 'kkkkkkkkkk'),
            .env = 'tools')
})

test_that('Testing "expects" case using with_mock on a function in external package', {
  stub_builder <- stub(sub)
  stub_builder$strictlyExpects(pattern = "([^.]+)\\.[[:alnum:]]+$", replacement = '\\1', x = 'goo.goo',
                       ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)
  sub_stub <- stub_builder$f

  with_mock(sub = sub_stub,
            expect_error(tools::file_path_sans_ext('dsfsd')),
            expect_silent(tools::file_path_sans_ext('goo.goo')),
            .env = 'tools')
})

test_that('Testing non-simple cases using with_mock on a function in external package', {
  stub_builder <- stub(sub)

  stub_builder$onCall(1)$returns('yay!')

  stub_builder$onCall(2)$strictlyExpects(
    pattern = "([^.]+)\\.[[:alnum:]]+$", replacement = '\\1', x = 'test2',
    ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)$returns(10)

  stub_builder$onCall(3)$strictlyExpects(
    pattern = "([^.]+)\\.[[:alnum:]]+$", replacement = '\\1', x = 'test3',
    ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)$throws('err msg')

  stub_builder$onCall(4)$expects(x = 'test4')$returns('test4-res')

  stub_builder$withExactArgs(
    pattern = "[.](gz|bz2|xz)$", replacement = '', x = 'test5',
    ignore.case = FALSE, perl = FALSE, fixed = FALSE, useBytes = FALSE)$throws('test5-res')

  stub_builder$withArgs(x = 'test67')$returns('test67-res')

  stub_builder$withArgs(x = 'test67-res')$returns('test67-res-res')

  sub_stub <- stub_builder$f

  with_mock(sub = sub_stub,
            expect_equal(tools::file_path_sans_ext('dsfsdfs.gfg'), 'yay!'),
            .env = 'tools')

  with_mock(sub = sub_stub,
            expect_equal(tools::file_path_sans_ext('test2'), 10),
            .env = 'tools')

  with_mock(sub = sub_stub,
            expect_error(tools::file_path_sans_ext('test3'), 'err msg'),
            .env = 'tools')

  with_mock(sub = sub_stub,
            expect_equal(tools::file_path_sans_ext('test4'), 'test4-res'),
            .env = 'tools')

  with_mock(sub = sub_stub,
            expect_error(tools::file_path_sans_ext('test5', compression = TRUE), 'test5-res'),
            .env = 'tools')

  with_mock(sub = sub_stub,
            expect_equal(tools::file_path_sans_ext('test67'), 'test67-res'),
            .env = 'tools')

  with_mock(sub = sub_stub,
            expect_equal(tools::file_path_sans_ext('test67', compression = TRUE), 'test67-res-res'),
            .env = 'tools')

})