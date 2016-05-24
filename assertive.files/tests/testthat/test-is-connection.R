test_that("test.is_bzfile_connection.some_connections.returns_true_for_bzfile_connections",
  {
    bzcon <- bzfile("foo")
    on.exit(close(bzcon))
    expect_true(is_bzfile_connection(bzcon))
    fcon <- file("foo")
    on.exit(close(fcon), add = TRUE)
    expect_false(is_bzfile_connection(fcon))
  })

test_that("test.is_connection.a_connection.returns_true", {
  fcon <- file()
  on.exit(close(fcon))
  expect_true(is_connection(fcon))
})

test_that("test.is_connection.not_a_connection.returns_false", {
  expect_false(is_connection("not a connection"))
})

test_that("test.is_connection.std_connections.returns_true", {
  for (con in c(stdin, stdout, stderr)) {
    expect_true(is_connection(con()))
  }
})

test_that("test.is_file_connection.some_connections.returns_true_for_file_connections",
{
  fcon <- file("foo")
  on.exit(close(fcon))
  expect_true(is_file_connection(fcon))
  gcon <- gzfile("foo")
  on.exit(close(gcon), add = TRUE)
  expect_false(is_file_connection(gcon))
})

test_that("test.is_gzfile_connection.some_connections.returns_true_for_gzfile_connections",
{
  gcon <- gzfile("foo")
  on.exit(close(gcon))
  expect_true(is_gzfile_connection(gcon))
  fcon <- file("foo")
  on.exit(close(fcon), add = TRUE)
  expect_false(is_gzfile_connection(fcon))
})

test_that("test.is_incomplete_connection.a_closed_connection.returns_false", 
  {
    tcon <- textConnection("txt", "w", local = TRUE)
    close(tcon)
    expect_false(is_incomplete_connection(tcon))
  })

test_that("test.is_incomplete_connection.a_complete_connection.returns_false", 
  {
    tcon <- textConnection("txt", "w", local = TRUE)
    on.exit(close(tcon))
    cat("this has a final newline character\n", file = tcon)
    expect_false(is_incomplete_connection(tcon))
  })

test_that("test.is_incomplete_connection.an_incomplete_connection.returns_true", 
  {
    tcon <- textConnection("txt", "w", local = TRUE)
    on.exit(close(tcon))
    cat("this has no final newline character", file = tcon)
    expect_true(is_incomplete_connection(tcon))
  })

test_that("test.is_incomplete_connection.not_a_connection.returns_false", {
  expect_false(is_incomplete_connection("not a connection"))
})

test_that("test.is_open_connection.a_closed_connection.returns_false", {
  fcon <- file()
  close(fcon)
  expect_false(is_open_connection(fcon))
})

test_that("test.is_open_connection.an_open_readable_connection.returns_true", 
  {
    readable <- "r"
    file.create(tmp <- tempfile())
    fcon <- file(tmp, open = readable)
    on.exit({
      close(fcon)
      unlink(tmp)
    })
    expect_true(is_open_connection(fcon, readable))
  })

test_that("test.is_open_connection.an_open_writable_connection.returns_true", 
  {
    writable <- "w"
    file.create(tmp <- tempfile())
    fcon <- file(tmp, open = writable)
    on.exit({
      close(fcon)
      unlink(tmp)
    })
    expect_true(is_open_connection(fcon, writable))
  })

test_that("test.is_open_connection.not_a_connection.returns_false", {
  expect_false(is_open_connection("not a connection"))
}) 

test_that("test.is_pipe_connection.some_connections.returns_true_for_pipe_connections",
{
  #Insert Magritte jokes here.
  pcon <- pipe("foo")
  on.exit(close(pcon))
  expect_true(is_pipe_connection(pcon))
  fcon <- file("foo")
  on.exit(close(fcon), add = TRUE)
  expect_false(is_pipe_connection(fcon))
})

test_that("test.is_readable_connection.some_connections.returns_true_for_readable_connections",
{
  expect_true(is_readable_connection(stdin()))
  expect_false(is_readable_connection(stdout()))
  file.create(tmp <- tempfile())
  fcon <- file(tmp, "r")
  on.exit({
    close(fcon)
    unlink(tmp)
  })
  expect_true(is_readable_connection(fcon))
})

test_that("test.is_pipe_connection.some_connections.returns_true_for_pipe_connections",
{
  #Insert Magritte jokes here.
  pcon <- pipe("foo")
  on.exit(close(pcon))
  expect_true(is_pipe_connection(pcon))
  fcon <- file("foo")
  on.exit(close(fcon), add = TRUE)
  expect_false(is_pipe_connection(fcon))
})

test_that("test.is_stderr.some_connections.returns_true_for_stderr",
{
  expect_true(is_stderr(stderr()))
  expect_false(is_stderr(stdin()))
  expect_false(is_stderr(stdout()))
})

test_that("test.is_stdin.some_connections.returns_true_for_stdin",
{
  expect_false(is_stdin(stderr()))
  expect_true(is_stdin(stdin()))
  expect_false(is_stdin(stdout()))
})

test_that("test.is_stdout.some_connections.returns_true_for_stdout",
{
  expect_false(is_stdout(stderr()))
  expect_false(is_stdout(stdin()))
  expect_true(is_stdout(stdout()))
})

test_that("test.is_terminal_connection.some_connections.returns_true_for_terminal_connections",
{
  expect_true(is_terminal_connection(stderr()))
  expect_true(is_terminal_connection(stdin()))
  expect_true(is_terminal_connection(stdout()))
  fcon <- file("foo")
  on.exit(close(fcon))
  expect_false(is_terminal_connection(fcon))
})

test_that("test.is_text_connection.some_connections.returns_true_for_text_connections",
{
  tcon <- textConnection("foo")
  on.exit(close(tcon))
  expect_true(is_text_connection(tcon))
  fcon <- file("foo")
  on.exit(close(fcon), add = TRUE)
  expect_false(is_text_connection(fcon))
})

test_that("test.is_unz_connection.some_connections.returns_true_for_unz_connections",
{
  ucon <- unz("foo", "bar")
  on.exit(close(ucon))
  expect_true(is_unz_connection(ucon))
  fcon <- file("foo")
  on.exit(close(fcon), add = TRUE)
  expect_false(is_unz_connection(fcon))
})

test_that("test.is_url_connection.some_connections.returns_true_for_url_connections",
{
  ucon <- url("http://shop.oreilly.com/product/0636920028352.do")
  on.exit(close(ucon))
  expect_true(is_url_connection(ucon))
  fcon <- file("foo")
  on.exit(close(fcon), add = TRUE)
  expect_false(is_url_connection(fcon))
})

test_that("test.is_writable_connection.some_connections.returns_true_for_writable_connections",
{
  expect_false(is_writable_connection(stdin()))
  expect_true(is_writable_connection(stdout()))
  file.create(tmp <- tempfile())
  fcon <- file(tmp, "w")
  on.exit({
    close(fcon)
    unlink(tmp)
  })
  expect_true(is_writable_connection(fcon))
})

test_that("test.is_xzfile_connection.some_connections.returns_true_for_xzfile_connections",
{
  xcon <- xzfile("foo")
  on.exit(close(xcon))
  expect_true(is_xzfile_connection(xcon))
  fcon <- file("foo")
  on.exit(close(fcon), add = TRUE)
  expect_false(is_xzfile_connection(fcon))
})
