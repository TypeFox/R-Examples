library("testthat")
library("rfigshare")

# This loads the rOpenSci figshare sandbox credentials, so that the example 
# can run automatically during check and install.  Unlike normal figshare accounts,
# data loaded to this testing sandbox is periodically purged.  

authenticated <- fs_auth(token = "xdBjcKOiunwjiovwkfTF2QjGhROeLMw0y0nSCSgvg3YQxdBjcKOiunwjiovwkfTF2Q", token_secret = "4mdM3pfekNGO16X4hsvZdg")

if(!is(authenticated, "Token1.0")){
  warning("Could not authenticate the test server, skipping online tests")
} else {

  context("Authentication, creating authors, and content")

  test_that("We are able to create objects on figshare", {

# fs_create
      new_object <- fs_create("My Title", "A description of the object", "dataset")
      expect_is(new_object, "numeric")
      fs_delete(new_object)	  
# fs_create_author
      expect_error(fs_create_author()) # Can't create an empty author
  })

# ------------------------------------------------------------------


  context("Author creation")
  test_that("we can add new authors without accounts", {

    # We need a unique name...
    if(require(uuid)){
      full_name = paste("John", UUIDgenerate())

      author_id <- fs_create_author(full_name)
      expect_is(author_id, "numeric")

      # If we try to create the same author, 
      expect_warning(id <- fs_create_author(full_name))
     # expect to get back the original id  ## NOT implemented in the figshare API yet! Emailed Mark 
     # expect_equal(id, author_id)

      MrX <- paste("John", UUIDgenerate())
      id2 <- fs_new_article(title="title", description="description", 
                     authors = c("Karthik Ram", MrX))
      d <- fs_details(id2)

      # Expect that Mr X is now the third author FIXME we seem to confuse the API sometimes here
#  expect_match(sapply(d$authors, `[[`, "first_name")[[3]], strsplit(MrX, ' ')[[1]][1])
    }


  })



  context("Data Retrieval")

  test_that("Downloads work correctly", {

    # Downloading a csv file
    url <- fs_download(90818)
    data <- read.csv(url)
    expect_that(data, is_a("data.frame"))

    #downloading an Excel file
    library(gdata)
    excel_url <- fs_download(894949)
    excel_data <- gdata::read.xls(excel_url)
    expect_that(excel_data, is_a("data.frame"))
  })




# ----------------------------------------------------------------
  context("Metadata to/from existing objects")



# fs_update
  test_that("We are able to add metadata to existing objects", {
   data <- head(iris)
   new_fs_obj <- fs_create("Dummy data", "Fisher's iris", "dataset")
# fs_add_authors	  
    fs_add_authors(new_fs_obj, "Scott Chamberlain")

    fs_add_tags(new_fs_obj, "a random tag")
        
    fs_add_categories(new_fs_obj, "Ecology")
    
    fs_add_links(new_fs_obj, "http://ropensci.org")

    d <- fs_details(new_fs_obj)

    expect_equal(d$tags[[1]]$name, "a random tag")
    expect_equal(d$categories[[1]]$name, "Ecology")
    expect_equal(d$links[[1]]$link, "http://ropensci.org")
    expect_equal(d$title, "Dummy data")
     expect_equal(d$authors[[2]]$last_name, "Chamberlain")
   
    fs_update(new_fs_obj, title = "This is now the new title")
    d <- fs_details(new_fs_obj)
    expect_equal(d$title, "This is now the new title")

    # A blank category is undefined, but a blank tag is just ommitted 
    expect_error(fs_add_categories(new_fs_obj, ""))

    # A missing argument is an error
    expect_error(fs_add_tags(new_fs_obj))





  })

  test_that("articles are tagged as 'Published using rfigshare'", {


   write.csv(mtcars, "mtcars.csv")
   id <- fs_new_article(title="A Test of rfigshare",
                           description="This is a test",
                           type="dataset",
                           authors=c("Karthik Ram", "Scott Chamberlain"),
                           tags=c("ecology", "openscience"),
                           categories="Ecology",
                           links="http://ropensci.org",
                           files="mtcars.csv",
                           visibility="private")

    details <- fs_details(id, mine = TRUE)

    expect_match("Published using rfigshare", sapply(details$tags, `[[`, "name")[[1]])

    fs_delete(id)
    unlink("mtcars.csv")
  })

  test_that("We can browse our own files", {
    expect_is(fs_browse(mine = TRUE), "list")
  })

# ----------------------------------------------------------------
  context("Searching figshare")
  test_that("fs_category_list works and debugs", {
          expect_is(fs_category_list(), "table")

          expect_is(fs_category_list(debug=TRUE), "response")

            })
#  test_that("We are able to perform browsing correctly", {
#            expect_is(fs_browse()[[1]], "fs_object")

#  })

  test_that("We are able to perform search correctly", {
 #           expect_is(fs_search("Boettiger")[[1]], "fs_object")

            ids <- fs_ids(fs_search("Boettiger", category="Ecology"))
            expect_is(ids, "numeric")

  })

  test_that("We are able to perform author search correctly", {
       auth <-  fs_author_search("Karthik Ram")
       expect_equal(auth[[1]]$id, 97306)
      

  })


# ----------------------------------------------------------------
  context("Manipulating objects on figshare")
# fs_make_public
# fs_make_private
# fs_delete (delete second)
  test_that("We are able to change visibility of objects", {
      data <- head(iris)
      new_fs_obj <- fs_create("Dummy data", "Fisher's iris", "dataset")
      write.csv(data, file = "iris_data.csv")
      fs_upload(new_fs_obj, "iris_data.csv")

      response <- fs_make_private(new_fs_obj)
      expect_is(response, "response")
      expect_equal(response$status_code, 400)
      fs_add_categories(new_fs_obj, "Ecology")
      expect_equal(fs_make_private(new_fs_obj)$status_code, 200)
      fs_delete(new_fs_obj)
      unlink("iris_data.csv")
  })

  test_that("We can make articles public", {
      data <- head(iris)
      new_fs_obj <- fs_create("Dummy data", "Fisher's iris", "dataset")
      write.csv(data, file = "iris_data.csv")
      fs_upload(new_fs_obj, "iris_data.csv")
      expect_equal(fs_make_public(new_fs_obj)$status_code, 400)
      fs_add_categories(new_fs_obj, "Ecology")
      expect_equal(fs_make_public(new_fs_obj)$status_code, 200)
      fs_delete(new_fs_obj)
      unlink("iris_data.csv")

  })
# ----------------------------------------------------------------
  context("We can obtain summaries")
# plot_to_filename
# summary_fs_details
  test_that("We can obtain summaries", {
#    fs_auth() # needs interactive environment
      my_article <- fs_details(138)
      summary_fs_details(my_article)
      expect_that(summary_fs_details(my_article), prints_text("Article ID"))
     

  })


# ----------------------------------------------------------------
  context("Testing miscellaneous functions")
  test_that("We can save plots to disk while uploading", {
      library(ggplot2)
      p <- qplot(mpg, wt, data = mtcars)
      plot_to_filename(p, "myfilename", ".")
      expect_true(file.exists("myfilename.png"))
      unlink("myfilename.png")

  })


# miscellaneous functions
# fs_embed

  test_that("we can get the url of an image", {
    path <- fs_image_url(138)
    library(httr)
    resp <- GET(path)
    expect_equal(resp$headers$`content-type`, "image/png")
  })

# fs_ids
# Can't test because it's broken
# also missing tests for fs_embed
}
