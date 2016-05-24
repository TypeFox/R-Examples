library(Rexperigen)
context("Utility functions")

test_that("URL preparation works", {
    expect_equal(prepare.server.URL("http://alma.com/"), "http://alma.com/")
    expect_equal(prepare.server.URL("https://alma.com/"), "https://alma.com/")
    expect_equal(prepare.server.URL("http://alma.com"), "http://alma.com/")
    expect_equal(prepare.server.URL("alma.com"), "http://alma.com/")
    expect_equal(prepare.server.URL("alma.com/"), "http://alma.com/")
    expect_equal(create.API.request.URL("alma.com", "print", list(a=1, b="korte")),
                 "http://alma.com/print?a=1&b=korte")
})

test_that("API requests", {
    with_mock(
        "RCurl::getURL" = function(url, customrequest){
            expect_equal(url, "http://db.phonologist.org/version")
            expect_equal(customrequest, "GET")
        },
        API.request()
    )
    with_mock(
        "RCurl::getURL" = function(url, customrequest){
            expect_equal(url, "http://exp.erig.en/getdata?a=155&b=telephone")
            expect_equal(customrequest, "PUT")
        },
        API.request(server = "exp.erig.en",
                    request = "getdata",
                    params = list(a="155", b="telephone"),
                    method = "PUT"
                    )
    )
})

test_that("Authed API request", {
    options(Rexperigen.experimenter = "alma")
    options(Rexperigen.password = "korte")
    with_mock(
        "RCurl::getURL" = function(url, username, password, httpauth, customrequest){
            expect_equal(url, "http://exp.erig.en/login?a=me")
            expect_equal(customrequest, "GET")
            expect_equal(httpauth, RCurl::AUTH_DIGEST)
            expect_equal(username, "alma")
            expect_equal(password, "korte")
        },
        API.request(server = "exp.erig.en",
                    request = "login",
                    params = list(a="me"),
                    method = "GET",
                    auth = TRUE
                    )
    )    
    options(Rexperigen.experimenter = "")
    options(Rexperigen.password = "")
})

test_that("serverVersion", {
    with_mock(
        "RCurl::getURL" = function(...){
            "2.1.3"
        },
        "RCurl::url.exists" = function(url){
            expect_equal(url, "http://db.phonologist.org/version")
            TRUE
        },
        expect_equal(server.version(), "2.1.3")
    )

    with_mock(
        "RCurl::url.exists" = function(url){
            if(url == "http://a/version") FALSE
            else TRUE
        },
        expect_equal(server.version("a"), "1.0.0")
    )

    with_mock(
        "RCurl::url.exists" = function(url){
            FALSE
        },
        expect_equal(server.version(), NO_SERVER_ERROR)
    )
})

test_that("check authentication", {
    ## We start with 1.0.0 server
    r1 <- checkAuthentication("test", FALSE, 1)
    expect_equal(r1$request, "test.cgi")
    expect_equal(r1$auth, FALSE)
    expect_warning({
        r2 <- checkAuthentication("tost", TRUE, 1)
    }, "not supported")
    expect_equal(r2$request, "tost.cgi")
    expect_equal(r2$auth, FALSE)
    options(Rexperigen.server.version="4.6.1")
    r1 <- checkAuthentication("tist", FALSE, 1)
    expect_equal(r1$request, "tist")
    expect_equal(r1$auth, FALSE)
    r1 <- checkAuthentication("tast", TRUE)
    expect_equal(r1$request, "digest/tast", 1)
    expect_equal(r1$auth, TRUE)
    options(Rexperigen.server.version="1.0.0")
})

test_that("cleanurl", {
    expect_error(cleanURL("al/ba"), "not supported")
    options(Rexperigen.server.version="4.6.1")
    with_mock(
        "RCurl::getURL" = function(req, ...){
            expect_equal(req, "http://db.phonologist.org/cleanURL?sourceurl=al%2Fba")
            "al.ba"
        },
        expect_equal(cleanURL("al/ba"), "al.ba")
    )
})


options(Rexperigen.server.version="1.0.0")
