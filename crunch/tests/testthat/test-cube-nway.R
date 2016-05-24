context("Cubes with >2 dimensions")

if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(newDatasetFromFixture("apidocs")), {
            names(categories(ds$q1)) <- LETTERS[1:5] ## To distinguish from other vars
            test_that("Cat x datetime x subvar", {
                kube <- crtabs(~ q1 + wave + petloc$Home, data=ds)
                expect_equivalent(as.array(kube)[,,"Cat"],
                    array(c(1, 1, 0,
                            2, 0, 0),
                        dim=c(3L, 2L),
                        dimnames=list(q1=c("A", "B", "C"),
                            wave=c("2014-12-01", "2015-01-01"))))
                expect_equivalent(as.array(kube)[,,"Dog"],
                    array(c(1, 0, 0,
                            0, 1, 1),
                        dim=c(3L, 2L),
                        dimnames=list(q1=c("A", "B", "C"),
                            wave=c("2014-12-01", "2015-01-01"))))
                expect_equivalent(as.array(kube)[,,"Bird"],
                    array(c(0, 0, 0,
                            1, 0, 0),
                        dim=c(3L, 2L),
                        dimnames=list(q1=c("A", "B", "C"),
                            wave=c("2014-12-01", "2015-01-01"))))
            })

            test_that("Categorical Array cube", {
                kube <- crtabs(~ petloc, data=ds)
                expect_equivalent(as.array(kube),
                    array(c(5, 6,
                            3, 4,
                            3, 6),
                        dim=c(2L, 3L),
                        dimnames=list(petloc=c("Home", "Work"),
                            petloc=c("Cat", "Dog", "Bird"))))
                kube@useNA <- "always"
                expect_equivalent(as.array(kube),
                    array(c(5, 6,
                            3, 4,
                            3, 6,
                            4, 3,
                            5, 1),
                        dim=c(2L, 5L),
                        dimnames=list(petloc=c("Home", "Work"),
                            petloc=c("Cat", "Dog", "Bird", "Skipped", "Not Asked"))))
            })
            test_that("CA x categorical", {
                kube <- crtabs(~ petloc + country, data=ds)
                expect_equivalent(as.array(kube)[,,"Belgium"],
                    array(c(2, 2,
                            0, 1,
                            0, 0),
                        dim=c(2L, 3L),
                        dimnames=list(petloc=c("Home", "Work"),
                            petloc=c("Cat", "Dog", "Bird"))))
            })
        })
    })
}
