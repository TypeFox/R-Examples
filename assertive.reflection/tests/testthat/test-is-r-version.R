test_that("test.is_r_alpha.any_r.returns_true_if_is_r_alpha", 
{
  expected <- version$status == "alpha"
  actual <- is_r_alpha()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(
      cause(actual), 
      noquote(
        sprintf(
          "R's build type is %s, not alpha.",
          clean_status_string()
        )
      )
    )
  }
})

test_that("test.is_r_beta.any_r.returns_true_if_is_r_beta", 
{
  expected <- version$status == "beta"
  actual <- is_r_beta()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(
      cause(actual), 
      noquote(
        sprintf(
          "R's build type is %s, not beta.",
          clean_status_string()
        )
      )
    )
  }
})

test_that("test.is_r_devel.any_r.returns_true_if_is_r_devel", {
  expected <- version$status == "Under development (unstable)"
  actual <- is_r_devel()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(
      cause(actual), 
      noquote(
        sprintf(
          "R's build type is %s, not development.",
          clean_status_string()
        )
      )
    )
  }
})

test_that("test.is_r_patched.any_r.returns_true_if_is_r_patched", {
  expected <- version$status == "Patched"
  actual <- is_r_patched()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(
      cause(actual), 
      noquote(
        sprintf(
          "R's build type is %s, not patched.",
          clean_status_string()
        )
      )
    )
  }
})

test_that("test.is_r_release.any_r.returns_true_if_is_r_release", 
{
  expected <- version$status == ""
  actual <- is_r_release()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(
      cause(actual), 
      noquote(
        sprintf(
          "R's build type is %s, not release.",
          clean_status_string()
        )
      )
    )
  }
})

test_that("test.is_r_release_candidate.any_r.returns_true_if_is_r_release_candidate", 
{
  expected <- version$status == "RC"
  actual <- is_r_release_candidate()
  expect_equal(strip_attributes(actual), expected)
  if (!actual) {
    expect_equal(
      cause(actual), 
      noquote(
        sprintf(
          "R's build type is %s, not release candidate.",
          clean_status_string()
        )
      )
    )
  }
})
