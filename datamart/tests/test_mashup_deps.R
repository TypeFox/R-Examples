#library(datamart)

test_mashup_deps <- function() {
    afun <- function(b, ...) b
    bfun <- function(...) "b"
    
    m <- datamart(
        resfunc(resource="a", depends="b", fun=afun),
        resfunc(resource="b", fun=bfun)
    )
    query(m, "a") == "b"

    # all(
        # class(query(m, "iris"))=="data.frame",
        # nrow(query(m, "iris"))==150,
        # all.equal(
            # query(m, "means"),
            # setNames(c(5.843333, 3.057333, 1.199333, 3.758000), c("Sepal.Length", "Sepal.Width", "Petal.Width", "Petal.Length")),
            # tolerance=10^(-4)
        # )
    # )
}


