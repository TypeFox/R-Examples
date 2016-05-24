context("Variable grouping and order setting")

test_that("VariableGroup and Order objects can be made", {
    expect_true(inherits(VariableGroup(group="group1", entities=""),
        "VariableGroup"))
    expect_true(inherits(VariableGroup(name="group1", entities=""),
        "VariableGroup"))
    vg1 <- VariableGroup(name="group1", entities="")
    expect_true(inherits(VariableOrder(vg1), "VariableOrder"))
    expect_true(inherits(VariableOrder(list(name="group1",
        entities=""), vg1), "VariableOrder"))
})

with_mock_HTTP({
    test.ds <- loadDataset("test ds")
    varcat <- allVariables(test.ds)

    test_that("ordering methods on variables catalog", {
        expect_true(inherits(ordering(variables(test.ds)), "VariableOrder"))
        expect_true(inherits(ordering(test.ds), "VariableOrder"))
        expect_identical(ordering(variables(test.ds)), ordering(test.ds))
    })

    test_that("relative URLs in hierarchical order", {
        vc <- varcat
        vc@orders$hier <- sub("hierarchical",
            "relative-hierarchical", vc@orders$hier)
        expect_identical(entities(vc@order),
            entities(VariableOrder(crGET(vc@orders$hier))))
    })

    test.ord <- ordering(test.ds)
    ent.urls <- urls(test.ord)
    nested.ord <- try(VariableOrder(
        VariableGroup(name="Group 1",
            entities=list(ent.urls[1],
                        VariableGroup(name="Nested", entities=ent.urls[2:3]),
                        ent.urls[4])),
        VariableGroup(name="Group 2", entities=ent.urls[5:6])))

    test_that("Can extract group(s) by name", {
        expect_identical(nested.ord[["Group 2"]],
            VariableGroup(name="Group 2", entities=ent.urls[5:6]))
        expect_identical(nested.ord$`Group 2`,
            VariableGroup(name="Group 2", entities=ent.urls[5:6]))
        expect_identical(nested.ord["Group 2"],
            VariableOrder(VariableGroup(name="Group 2",
            entities=ent.urls[5:6])))
    })
    test_that("Can create nested groups", {
        expect_true(inherits(nested.ord, "VariableOrder"))
        expect_identical(urls(nested.ord), ent.urls)
    })
    test_that("Can read nested groups from the API", {
        vc <- varcat
        vc@orders$hier <- sub("hierarchical",
            "nested-hierarchical", vc@orders$hier)
        expect_identical(entities(nested.ord),
            entities(VariableOrder(crGET(vc@orders$hier))))
    })
    test_that("Nested groups can also have relative urls", {
        vc <- varcat
        vc@orders$hier <- sub("hierarchical",
            "relative-and-nested-hierarchical", vc@orders$hier)
        expect_identical(entities(nested.ord),
            entities(VariableOrder(crGET(vc@orders$hier))))
    })
    test_that("Nested groups can serialize and deserialize", {
        vglist <- cereal(nested.ord)
        expect_identical(vglist, list(graph=list(
            list(`Group 1`=list(
                    ent.urls[1],
                    list(`Nested`=as.list(ent.urls[2:3])),
                    ent.urls[4]
                    )
                ),
            list(`Group 2`=as.list(ent.urls[5:6]))
        )))
    })

    ng <- list(ent.urls[1],
                VariableGroup(name="Nested", entities=ent.urls[2:3]),
                ent.urls[4])
    test_that("can assign nested groups in entities", {
        to <- test.ord
        try(entities(to) <- ng)
        expect_identical(entities(to), entities(ng))
        expect_identical(urls(to), ent.urls[1:4])
        expect_identical(to[[2]],
            VariableGroup(name="Nested", entities=ent.urls[2:3]))
        expect_identical(entities(to[[2]]), as.list(ent.urls[2:3]))
    })
    test_that("can assign group into order", {
        to <- test.ord
        try(to[[1]] <- VariableGroup(name="[[<-", entities=ng))
        expect_identical(entities(to[[1]]), ng)
        expect_identical(name(to[[1]]), "[[<-")
        expect_identical(urls(to[[1]]), ent.urls[1:4])
        expect_identical(to[[1]][[2]],
            VariableGroup(name="Nested", entities=ent.urls[2:3]))
    })
    test_that("can assign NULL into order to remove a group", {
        no <- no2 <- no3 <- nested.ord
        no[[2]] <- NULL
        expect_identical(no, VariableOrder(
            VariableGroup(name="Group 1", entities=list(ent.urls[1],
                VariableGroup(name="Nested", entities=ent.urls[2:3]),
                ent.urls[4]))))
        no2[["Group 2"]] <- NULL
        expect_identical(no2, VariableOrder(
            VariableGroup(name="Group 1", entities=list(ent.urls[1],
                VariableGroup(name="Nested", entities=ent.urls[2:3]),
                ent.urls[4]))))
        no3$`Group 2` <- NULL
        expect_identical(no3, VariableOrder(
            VariableGroup(name="Group 1", entities=list(ent.urls[1],
                VariableGroup(name="Nested", entities=ent.urls[2:3]),
                ent.urls[4]))))
    })
    test_that("can assign group into group by index", {
        to <- test.ord
        try(to[[1]] <- VariableGroup(name="[[<-", entities=ng))
        expect_identical(to[[1]][[1]], ent.urls[1])
        try(to[[1]][[1]] <- VariableGroup(name="Nest2",
            entities=to[[1]][[1]]))
        expect_identical(entities(to[[1]]),
            list(VariableGroup(name="Nest2", entities=ent.urls[1]),
                VariableGroup(name="Nested", entities=ent.urls[2:3]),
                ent.urls[4]))
        expect_identical(urls(to[[1]]), ent.urls[1:4])
    })
    test_that("can assign into a nested group", {
        to <- test.ord
        try(to[[1]] <- VariableGroup(name="[[<-", entities=ng))
        try(entities(to[[1]][[2]]) <- rev(entities(to[[1]][[2]])))
        expect_identical(entities(to[[1]]),
            list(ent.urls[1],
                VariableGroup(name="Nested", entities=ent.urls[3:2]),
                ent.urls[4]))
        expect_identical(urls(to[[1]]), ent.urls[c(1,3,2,4)])
        expect_identical(name(to[[1]]), "[[<-")
        try(name(to[[1]]) <- "Something better")
        expect_identical(name(to[[1]]), "Something better")
    })

    test_that("Assignment by new group name", {
        nested.o <- nested.ord
        nested.o[["Group 3"]] <- test.ds["starttime"]
        expect_identical(names(nested.o), c("Group 1", "Group 2", "Group 3"))
        expect_identical(entities(nested.o[["Group 3"]]),
            list(self(test.ds$starttime)))
        ## Test the "duplicates option": starttime should have been removed from
        ## Group 2
        expect_identical(entities(nested.o[["Group 2"]]),
            list(self(test.ds$catarray)))
    })

    test_that("Assignment by new group name with duplicates", {
        nested.o <- nested.ord
        duplicates(nested.o) <- TRUE
        expect_true(duplicates(nested.o))
        nested.o[["Group 3"]] <- test.ds["starttime"]
        expect_identical(names(nested.o), c("Group 1", "Group 2", "Group 3"))
        expect_identical(entities(nested.o[["Group 3"]]),
            list(self(test.ds$starttime)))
        ## Test the "duplicates option": starttime should not have been removed
        ## from Group 2
        expect_identical(entities(nested.o[["Group 2"]]),
            list(self(test.ds$starttime), self(test.ds$catarray)))
    })

    test_that("Update group with Dataset", {
        nested.o <- nested.ord
        nested.o[["Group 2"]] <- test.ds[c("gender", "starttime")]
        expect_identical(entities(nested.o[["Group 2"]]),
            lapply(test.ds[c("gender", "starttime")], self))
    })

    test_that("Assignment by new nested group name", {
        nested.o <- nested.ord
        nested.o[["Group 1"]][[2]][["More nesting"]] <- self(test.ds$gender)
        expect_identical(entities(nested.o[["Group 1"]]$Nested[["More nesting"]]),
            list(self(test.ds$gender)))
        ## Test duplicates option: gender should only be in "More nesting"
        expect_identical(nested.o[["Group 1"]]$Nested[[1]],
            self(test.ds$mymrset))
    })

    test_that("Assignment by new nested group name with duplicates", {
        nested.o <- nested.ord
        duplicates(nested.o) <- TRUE
        expect_true(duplicates(nested.o[["Group 1"]]))
        expect_true(duplicates(nested.o[["Group 1"]][[2]]))
        nested.o[["Group 1"]][[2]][["More nesting"]] <- self(test.ds$gender)
        expect_identical(entities(nested.o[["Group 1"]]$Nested[["More nesting"]]),
            list(self(test.ds$gender)))
        ## Test duplicates option: gender should still be in "More nesting"
        expect_identical(nested.o[["Group 1"]]$Nested[[1]],
            self(test.ds$gender))
        ## Test that duplicates option passes to new group
        expect_true(duplicates(nested.o[["Group 1"]][[2]][["More nesting"]]))
    })

    test_that("c() Order/Group", {
        skip("TODO")
    })

    test_that("Update group with URLs", {
        skip("TODO")
        to <- test.ord
        expect_true(inherits(to, "VariableOrder"))
        expect_error(to[[1]] <- ent.urls,
            "Correct error expectation here")
        ## Now try where [[1]] is a Group
        nested.o <- nested.ord
        try(nested.o[[1]] <- ent.urls)
        expect_identical(entities(nested.o[[1]]), as.list(ent.urls))
    })

    test_that("VariableOrder/Group show methods", {
        expect_identical(showVariableOrder(nested.ord, vars=variables(test.ds)),
            c("[+] Group 1",
              "    Birth Year",
              "    [+] Nested",
              "        Gender",
              "        mymrset",
              "    Text variable ftw",
              "[+] Group 2",
              "    starttime",
              "    Cat Array"))
        no <- nested.ord
        no[[3]] <- VariableGroup("Group 3", entities=list())
        expect_identical(showVariableOrder(no, vars=variables(test.ds)),
            c("[+] Group 1",
              "    Birth Year",
              "    [+] Nested",
              "        Gender",
              "        mymrset",
              "    Text variable ftw",
              "[+] Group 2",
              "    starttime",
              "    Cat Array",
              "[+] Group 3",
              "    (Empty group)"))
    })

    ord <- test.ord
    test_that("Composing a VariableOrder step by step: setup", {
        expect_identical(showVariableOrder(ord, vars=variables(test.ds)),
            c("Birth Year",
              "Gender",
              "mymrset",
              "Text variable ftw",
              "starttime",
              "Cat Array"))
    })
    test_that("Composing a VariableOrder step by step: group 1 by dataset", {
        ord$Demos <<- test.ds[c("gender", "birthyr")]
        expect_identical(showVariableOrder(ord, vars=variables(test.ds)),
            c("mymrset",
              "Text variable ftw",
              "starttime",
              "Cat Array",
              "[+] Demos",
              "    Gender",
              "    Birth Year"))
    })
    test_that("Composing a VariableOrder step by step: group by Order subset", {
        ord$Arrays <<- ord[c(1, 4)] #test.ds[c("mymrset", "catarray")]
        expect_identical(showVariableOrder(ord, vars=variables(test.ds)),
            c("Text variable ftw",
              "starttime",
              "[+] Demos",
              "    Gender",
              "    Birth Year",
              "[+] Arrays",
              "    mymrset",
              "    Cat Array"))
    })
    test_that("Composing a VariableOrder step by step: nested group by dataset", {
        ord$Demos[["Others"]] <<- test.ds[c("birthyr", "textVar")]
        expect_identical(showVariableOrder(ord, vars=variables(test.ds)),
            c("starttime",
              "[+] Demos",
              "    Gender",
              "    [+] Others",
              "        Birth Year",
              "        Text variable ftw",
              "[+] Arrays",
              "    mymrset",
              "    Cat Array"))
    })
    test_that("Composing a VariableOrder step by step: reorder group", {
        ord$Demos <<- ord$Demos[2:1]
        expect_identical(showVariableOrder(ord, vars=variables(test.ds)),
            c("starttime",
              "[+] Demos",
              "    [+] Others",
              "        Birth Year",
              "        Text variable ftw",
              "    Gender",
              "[+] Arrays",
              "    mymrset",
              "    Cat Array"))
    })
    test_that("Composing a VariableOrder step by step: reorder order", {
        ord <<- ord[3:1]
        expect_identical(showVariableOrder(ord, vars=variables(test.ds)),
            c("[+] Arrays",
              "    mymrset",
              "    Cat Array",
              "[+] Demos",
              "    [+] Others",
              "        Birth Year",
              "        Text variable ftw",
              "    Gender",
              "starttime"))
    })
    test_that("Composing a VariableOrder step by step: nested group by Group", {
        ord$Arrays$MR <<- ord$Arrays[1]
        expect_identical(showVariableOrder(ord, vars=variables(test.ds)),
            c("[+] Arrays",
              "    Cat Array",
              "    [+] MR",
              "        mymrset",
              "[+] Demos",
              "    [+] Others",
              "        Birth Year",
              "        Text variable ftw",
              "    Gender",
              "starttime"))
    })
})


if (run.integration.tests) {
    with(test.authentication, {
        with(test.dataset(df), {
            test_that("Can get VariableOrder from dataset", {
                expect_true(setequal(unlist(entities(ordering(ds))),
                    urls(allVariables(ds))))
            })
            test_that("Can construct VariableOrder from variables", {
                vg <- VariableOrder(
                    VariableGroup(name="Group 1",
                        variables=ds[c("v1", "v3", "v5")]),
                    VariableGroup(name="Group 2.5", entities=ds["v4"]),
                    VariableGroup(name="Group 2",
                        entities=ds[c("v6", "v2")]))
                vglist <- cereal(vg)
                expect_identical(vglist, list(graph=list(
                    list(`Group 1`=list(self(ds$v1), self(ds$v3), self(ds$v5))),
                    list(`Group 2.5`=list(self(ds$v4))),
                    list(`Group 2`=list(self(ds$v6), self(ds$v2)))
                )))
            })
            starting.vg <- vg <- VariableOrder(
                VariableGroup(name="Group 1",
                    entities=ds[c("v1", "v3", "v5")]),
                VariableGroup(name="Group 2.5", variables=ds["v4"]),
                VariableGroup(name="Group 2",
                    entities=ds[c("v6", "v2")]),
                duplicates=TRUE)

            test_that("Get urls from VariableOrder and Group", {
                expect_identical(urls(vg[[1]]),
                    c(self(ds$v1), self(ds$v3), self(ds$v5)))
                expect_identical(urls(vg),
                    c(self(ds$v1), self(ds$v3), self(ds$v5), self(ds$v4),
                    self(ds$v6), self(ds$v2)))
            })

            try(entities(vg[[2]]) <- self(ds$v2))
            test_that("Set URLs -> entities on VariableGroup", {
                expect_identical(urls(vg[[2]]), self(ds$v2))
                expect_identical(urls(vg),
                    c(self(ds$v1), self(ds$v3), self(ds$v5), self(ds$v2),
                    self(ds$v6)))
            })
            try(entities(vg[[2]]) <- list(ds$v3))
            test_that("Set variables -> entities on VariableGroup", {
                expect_identical(urls(vg[[2]]), self(ds$v3))
            })

            try(name(vg[[2]]) <- "Group 3")
            test_that("Set name on VariableGroup", {
                expect_identical(names(vg), c("Group 1", "Group 3", "Group 2"))
            })
            try(names(vg) <- c("G3", "G1", "G2"))
            test_that("Set names on VariableOrder", {
                expect_identical(names(vg), c("G3", "G1", "G2"))
            })

            try(vglist <- cereal(vg))
            test_that("VariableOrder to/fromJSON", {
                expect_identical(vglist, list(graph=list(
                    list(`G3`=list(self(ds$v1), self(ds$v3), self(ds$v5))),
                    list(`G1`=list(self(ds$v3))),
                    list(`G2`=list(self(ds$v6), self(ds$v2)))
                )))

                vg[1:2] <- vg[c(2,1)]
                expect_identical(cereal(vg), list(graph=list(
                    list(`G1`=list(self(ds$v3))),
                    list(`G3`=list(self(ds$v1), self(ds$v3), self(ds$v5))),
                    list(`G2`=list(self(ds$v6), self(ds$v2)))
                )))
            })

            original.order <- ordering(ds)
            test_that("Can set VariableOrder on dataset", {
                expect_false(identical(starting.vg, original.order))
                ordering(ds) <- starting.vg
                expect_identical(entities(grouped(ordering(ds))),
                    entities(starting.vg))
                expect_identical(entities(grouped(ordering(refresh(ds)))),
                    entities(starting.vg))
                expect_true(inherits(ungrouped(ordering(ds)), "VariableGroup"))
                expect_true(inherits(ungrouped(ordering(refresh(ds))),
                    "VariableGroup"))
                expect_identical(names(ordering(ds)),
                    c("Group 1", "Group 2.5", "Group 2"))

                ## Test that can reorder groups
                ordering(ds) <- starting.vg[c(2,1,3)]
                expect_identical(entities(grouped(ordering(ds))),
                    entities(starting.vg[c(2,1,3)]))
                expect_identical(names(ordering(ds)),
                    c("Group 2.5", "Group 1", "Group 2"))
                expect_identical(names(ordering(refresh(ds))),
                    c("Group 2.5", "Group 1", "Group 2"))

                ds <- refresh(ds)
                expect_false(identical(entities(ordering(variables(ds))),
                    entities(original.order)))
                ordering(variables(ds)) <- original.order
                expect_identical(entities(ordering(variables(ds))),
                    entities(original.order))
                expect_identical(entities(ordering(variables(refresh(ds)))),
                    entities(original.order))
            })

            test_that("A partial order results in 'ungrouped' variables", {
                ordering(ds) <- starting.vg[1:2]
                expect_true(inherits(grouped(ordering(ds)), "VariableOrder"))
                expect_identical(entities(grouped(ordering(ds))),
                    entities(starting.vg[1:2]))
                expect_true(inherits(ungrouped(ordering(ds)), "VariableGroup"))
                expect_true(setequal(unlist(entities(ungrouped(ordering(ds)))),
                    c(self(ds$v6), self(ds$v2))))
            })

            test_that("grouped and ungrouped within a group", {
                nesting <- VariableGroup("Nest", self(ds$v3))
                ordering(ds) <- starting.vg
                ordering(ds)[["Group 1"]][[2]] <- nesting
                ## Update fixture with duplicates=TRUE, as it should be found
                ## after setting on a duplicates=TRUE order
                duplicates(nesting) <- TRUE
                expect_identical(grouped(ordering(ds)[["Group 1"]]),
                    VariableGroup("Group 1", list(nesting), duplicates=TRUE))
                expect_identical(ungrouped(ordering(ds)[["Group 1"]]),
                    VariableGroup("ungrouped", list(self(ds$v1), self(ds$v5))))
            })

            test_that("Can manipulate VariableOrder that's part of a dataset", {
                ordering(ds) <- starting.vg
                expect_identical(names(ordering(ds)),
                    c("Group 1", "Group 2.5", "Group 2"))
                names(ordering(ds))[3] <- "Three"
                expect_identical(names(ordering(ds)),
                    c("Group 1", "Group 2.5", "Three"))
                expect_identical(names(grouped(ordering(ds))),
                    c("Group 1", "Group 2.5", "Three"))
            })

            test_that("duplicates property: setup", {
                expect_false(duplicates(ordering(ds)))
            })
            test_that("duplicates property set on order in dataset", {
                duplicates(ordering(ds)) <<- TRUE
                expect_true(duplicates(ordering(ds)))
            })
            test_that("duplicates property persists on refreshing dataset", {
                expect_true(duplicates(ordering(refresh(ds))))
            })
            ord <- ordering(ds)
            test_that("duplicates property persists on extracting order", {
                expect_true(duplicates(ord))
            })
            test_that("duplicates property persists on refreshing order", {
                skip("refresh method for VariableOrder not implemented")
                expect_true(duplicates(refresh(ord)))
            })
            test_that("duplicates property from order is set on assign to ds", {
                duplicates(ordering(ds)) <<- FALSE
                expect_false(duplicates(ordering(ds)))
                ordering(ds) <<- ord
                expect_true(duplicates(ordering(ds)))
            })

            test_that("ordering<- validation", {
                bad.vg <- starting.vg
                entities(bad.vg[[1]]) <- c(entities(bad.vg[[1]])[-2],
                    "/not/a/variable")
                expect_error(ordering(ds) <- bad.vg,
                    "Variable URL referenced in Order not present in catalog: /not/a/variable")
            })

            test_that("Creating VariableOrder with named list doesn't break", {
                bad.vg <- do.call(VariableOrder, c(sapply(names(starting.vg),
                    function (i) starting.vg[[i]], simplify=FALSE),
                    duplicates=TRUE))
                ## The list of entities is named because sapply default is
                ## USE.NAMES=TRUE, but the VariableOrder constructor should
                ## handle this
                ordering(ds) <- bad.vg
                expect_identical(ordering(ds)@graph,
                    starting.vg@graph)
            })
        })
    })
}
