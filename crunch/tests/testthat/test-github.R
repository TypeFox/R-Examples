context("GitHub version check")

test_that("version number comparison is correct", {
    expect_true(versionIsGreaterThan("1.2.2", "1.2.1"))
    expect_true(versionIsGreaterThan("1.2.2", "1.1.3"))
    expect_true(versionIsGreaterThan("1.2.2", "0.3.3"))
    expect_true(versionIsGreaterThan("2.0.0", "1.2.2"))
    expect_false(versionIsGreaterThan("1.2.2", "1.2.2"))
    expect_false(versionIsGreaterThan("1.2.1", "1.2.2"))
    expect_false(versionIsGreaterThan("1.1.3", "1.2.2"))
    expect_false(versionIsGreaterThan("1.2.2", "2.0.0"))
})

with_mock_HTTP({
    test_that("checkForNewVersion parses github json", {
        expect_identical(checkForNewVersion("github-versions.json", "1.5.3"),
            NULL)
        expect_identical(checkForNewVersion("github-versions.json", "1.5.1"),
            "1.5.3")
        expect_identical(checkForNewVersion("github-versions.json", "1.6.3"),
            NULL)
        ## Now that version is greater than 1.5.3:
        expect_identical(checkForNewVersion("github-versions.json"), NULL)
    })

    test_that("notifyIfNewVersion messages correctly", {
        expect_message(notifyIfNewVersion("github-versions.json", "1.5.1"),
            "There's a new version")
        expect_silent(notifyIfNewVersion("github-versions.json"))
    })
})

without_internet({
    test_that("notifyIfNewVersion doesn't hang if GitHub doesn't respond", {
        expect_silent(notifyIfNewVersion("github-versions.json", "1.5.1"))
    })
})

if (run.integration.tests) {

}
