context("sanitise")

test_that("sanitise correctly sanitises string", {
    expect_true(sanitise("$><|{}%&_#^~\\") == "\\$$>$$<$$|$\\{\\}\\%\\&\\_\\#\\verb|^|\\~{}$\\backslash$")  #tests existing stuff
    expect_true(sanitise("[]") == "{}[]")  #tests basic functionality we added
    expect_true(sanitise("{}[]") == "\\{\\}{}[]")  #tests the basic functionality added doesn't impact existing functionality
    expect_true(sanitise("\u00a3") == "\\pounds ")  #tests pounds
}) 
