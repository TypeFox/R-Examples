## ---- eval = FALSE-------------------------------------------------------
#  b2AuthorizeAccount(url = "https://api.backblaze.com/b2api/v1/b2_authorize_account",
#                     accountId = "yourAccountId", authorizationKey = "yourAuthorisationKey")

## ---- eval = FALSE-------------------------------------------------------
#  b2CreateBucket(bucketName = "a-unique-bucket-name", bucketType = "allPrivate")

## ---- eval = FALSE-------------------------------------------------------
#  b2ListBuckets()

## ---- eval = FALSE-------------------------------------------------------
#  b2DeleteBucket(bucketId = "aUniqueBucketId")

## ---- eval = FALSE-------------------------------------------------------
#  b2UpdateBucket(bucketId = "aUniqueBucketId", bucketType = "allPrivate-or-allPublic")

## ---- eval = FALSE-------------------------------------------------------
#  uploadUrlReturn <- b2GetUploadUrl(bucketId = "ff062d0e23056cb55226081c")
#  uploadUrl <- uploadUrlReturn$uploadUrl
#  authToken <- uploadUrlReturn$authorizationToken

## ---- eval = FALSE-------------------------------------------------------
#  b2UploadFile(authToken, uploadUrl, fileName = "blah.txt")

## ---- eval = TRUE--------------------------------------------------------
b2FileTypesLocation <- system.file("extdata", "b2FileTypes.rds", package = "backblazer")
b2FileTypes <- readRDS(b2FileTypesLocation)

## ---- eval = FALSE-------------------------------------------------------
#  b2ListFileNames(bucketId = "aUniqueBucketId")

## ---- eval = FALSE-------------------------------------------------------
#  b2ListFileVersions(bucketId = "aUniqueBucketId")

## ---- eval = FALSE-------------------------------------------------------
#  b2GetFileInfo(fileId = "a_unique_file_id")

## ---- eval = FALSE-------------------------------------------------------
#  b2HideFile(bucketId = "aUniqueBucketId", fileName = "blah.txt")

## ---- eval = FALSE-------------------------------------------------------
#  # Download file by name
#  b2DownloadFileByName(bucketName = "aUniqueBucketName", fileName = "blah.txt", overwrite = TRUE)
#  
#  # Download file by ID
#  b2DownloadFileById(fileId = "a_unique_file_id", overwrite = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  b2DeleteFileVersion(fileName = "blah.txt", fileId = "a_unique_file_id")

