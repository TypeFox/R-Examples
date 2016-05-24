library(testthat)
library(yummlyr)

context("Search Recipes")

test_that("Search with actual http requests", {
    app_id <- Sys.getenv("YUMMLY_APP_ID")
    app_key <- Sys.getenv("YUMMLY_APP_KEY")
    save_yummly_credentials(app_id, app_key)
    if (Sys.getenv("TRAVIS") != "" || Sys.getenv("APPVEYOR") != "") {
        result <- search_recipes("bacon")
        expect_is(result, "list")
        expect_equal(length(result), 5)
        result <- search_recipes("bacon", allowed_ingredient = "asparagus")
        expect_true(all(unlist(lapply(result$matches$ingredients, function(x) any(sapply(x, grepl, pattern="asparagus"))))))
        expect_error(search_recipes("bacon", allowed_ingredient = "asparagus2"))
        expect_warning(search_recipes("bacon", allowed_ingredient = "fried"))
        result <- search_recipes("bacon", excluded_ingredient = "asparagus")
        expect_false(any(unlist(lapply(result$matches$ingredients, function(x) any(sapply(x[[1]], grepl, pattern="asparagus"))))))
        expect_error(search_recipes("bacon", allowed_ingredient = "asparagus2"))
        expect_warning(search_recipes("bacon", allowed_ingredient = "fried"))
    }
})

save_yummly_credentials("APP_ID", "APP_KEY")
test_that("Search. search_words argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon"),
        expect_true(grepl("q=bacon$", result)),
        result <- search_recipes(c("onion", "bacon")),
        expect_true(grepl("q=onion%20bacon$", result))
    )
})

test_that("Search. require_pictures argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", require_picture = TRUE),
        expect_true(grepl("q=bacon&requirePictures=true$", result))
    )
})

test_that("Search. allowed_ingredient argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", allowed_ingredient="garlic"),
        expect_true(grepl("allowedIngredient\\[\\]=garlic$", result)),
        result <- search_recipes("bacon", allowed_ingredient=c("garlic", "asparagus")),
        expect_true(grepl("allowedIngredient\\[\\]=garlic&allowedIngredient\\[\\]=asparagus", result)),
        expect_error(search_recipes("bacon", allowed_ingredient = "asparagus2")),
        expect_warning(search_recipes("bacon", allowed_ingredient = "fried"))
    )
})

test_that("Search. excluded_ingredient argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", excluded_ingredient = c("garlic","asparagus")),
        expect_true(grepl("excludedIngredient\\[\\]=garlic&excludedIngredient\\[\\]=asparagus$", result)),
        expect_error(search_recipes("bacon", excluded_ingredient = "asparagus2")),
        expect_warning(search_recipes("bacon", excluded_ingredient = "fried")),
        result <- search_recipes("bacon", excluded_ingredient=c("onion soup mix", "asparagus")),
        expect_true(grepl("excludedIngredient\\[\\]=onion%20soup%20mix&excludedIngredient\\[\\]=asparagus$", result))
    )
})

test_that("Search. allowed_allergy argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", allowed_allergy =c("Dairy-Free", "Gluten-Free")),
        expect_true(grepl("allowedAllergy\\[\\]=396%5[Ee]Dairy-Free&allowedAllergy\\[\\]=393%5[Ee]Gluten-Free$", result)),
        expect_error(search_recipes("bacon", allowed_allergy = "Dairy2")),
        expect_warning(search_recipes("bacon", allowed_allergy = "Dairy"))
    )
})

test_that("Search. allowed_diet argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", allowed_diet =c("Pescetarian", "Lacto vegetarian")),
        expect_true(grepl("allowedDiet\\[\\]=390%5[Ee]Pescetarian&allowedDiet\\[\\]=388%5[Ee]Lacto%20vegetarian$", result)),
        expect_error(search_recipes("bacon", allowed_diet = "Lacto2")),
        expect_warning(search_recipes("bacon", allowed_diet = "Lacto"))
    )
})

test_that("Search. allowed_cuisine argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", allowed_cuisine =c("American")),
        expect_true(grepl("allowedCuisine\\[\\]=cuisine%5[Ee]cuisine-american", result)),
        expect_error(search_recipes("bacon", allowed_cuisine = "American2")),
        expect_warning(search_recipes("bacon", allowed_cuisine = "Ame"))
    )
})

test_that("Search. excluded_cuisine argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", excluded_cuisine =c("American")),
        expect_true(grepl("excludedCuisine\\[\\]=cuisine%5[Ee]cuisine-american", result)),
        expect_error(search_recipes("bacon", excluded_cuisine = "American2")),
        expect_warning(search_recipes("bacon", excluded_cuisine = "Ame"))
    )
})

test_that("Search. allowed_course argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", allowed_course =c("Appetizers")),
        expect_true(grepl("allowedCourse\\[\\]=course%5[Ee]course-Appetizers", result)),
        expect_error(search_recipes("bacon", allowed_course = "Appetizers2")),
        expect_warning(search_recipes("bacon", allowed_course = "Appeti"))
    )
})

test_that("Search. excluded_course argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", excluded_course =c("Appetizers")),
        expect_true(grepl("excludedCourse\\[\\]=course%5[Ee]course-Appetizers", result)),
        expect_error(search_recipes("bacon", excluded_course = "Appetizers2")),
        expect_warning(search_recipes("bacon", excluded_course = "Appeti"))
    )
})

test_that("Search. allowed_holiday argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", allowed_holiday =c("Thanksgiving")),
        expect_true(grepl("allowedHoliday\\[\\]=holiday%5[Ee]holiday-thanksgiving", result)),
        expect_error(search_recipes("bacon", allowed_holiday = "Thanksgiving2")),
        expect_warning(search_recipes("bacon", allowed_holiday = "Thanks"))
    )
})

test_that("Search. excluded_holiday argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", excluded_holiday =c("Thanksgiving")),
        expect_true(grepl("excludedHoliday\\[\\]=holiday%5[Ee]holiday-thanksgiving", result)),
        expect_error(search_recipes("bacon", excluded_holiday = "Thanksgiving2")),
        expect_warning(search_recipes("bacon", excluded_holiday = "Thanks"))
    )
})

test_that("Search. max_total_time argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", max_total_time = 5400),
        expect_true(grepl("&maxTotalTimeInSeconds=5400", result)),
        result <- search_recipes("bacon", max_total_time = "a"),
        expect_true(!grepl("maxTotalTimeInSeconds", result))
    )
})

test_that("Search. max_retults argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", max_results = 20),
        expect_true(grepl("&maxResult=20", result)),
        result <- search_recipes("bacon", max_results = "a"),
        expect_true(!grepl("maxResult", result))
    )
})

test_that("Search. start argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", start = 20),
        expect_true(grepl("&start=20", result)),
        result <- search_recipes("bacon", start = "a"),
        expect_true(!grepl("start", result))
    )
})

test_that("Search. nutrition argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", nutrition = list(Calcium=list(min=3, max=3.5))),
        expect_true(grepl("&nutrition.CA.min=3&nutrition.CA.max=3.5", result)),
        result <- search_recipes("bacon", nutrition = list(Calcium=list(max=3.5))),
        expect_true(grepl("&nutrition.CA.max=3.5", result)),
        result <- search_recipes("bacon", nutrition = list(Calcium=list(max=3.5), Cholesterol=list(min=2))),
        expect_true(grepl("&nutrition.CA.max=3.5&nutrition.CHOLE.min=2", result)),
        expect_error(search_recipes("bacon", nutrition = list(Calcium2=list(min=3, max=3.5)))),
        expect_error(search_recipes("bacon", nutrition = list(Calcium=list(a=3, max=3.5))))
    )
})

test_that("Search. flavor argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", flavor = list(sweet=list(min=0.1, max=1))),
        expect_true(grepl("&flavor.sweet.min=0.1&flavor.sweet.max=1", result)),
        result <- search_recipes("bacon", flavor = list(sweet=list(min=0.1))),
        expect_true(grepl("&flavor.sweet.min=0.1", result)),
        result <- search_recipes("bacon", flavor = list(sweet=list(min=0.1), sour=list(max=0.8))),
        expect_true(grepl("&flavor.sweet.min=0.1&flavor.sour.max=0.8", result)),
        expect_error(search_recipes("bacon", flavor = list(Calcium2=list(min=3, max=3.5)))),
        expect_error(search_recipes("bacon", flavor = list(sweet=list(min=3, max=3.5)))),
        expect_error(search_recipes("bacon", flavor = list(Calcium=list(a=3, max=3.5))))
    )
})

test_that("Search. facetField argument works (with mocks)", {
    with_mock(
        `yummlyr::perform_query` = function(query) query,
        `jsonlite::fromJSON` = function(content) content,
        result <- search_recipes("bacon", facet_field = "diet"),
        expect_true(grepl("&facetField\\[\\]=diet", result)),
        result <- search_recipes("bacon", facet_field = c("ingredient", "diet")),
        expect_true(grepl("&facetField\\[\\]=ingredient&facetField\\[\\]=diet", result)),
        expect_error(search_recipes("bacon", facet_field = "test"))
    )
})
