context("zip utility functions")

test_that('can get information about sample zip file',{
  
  if(nchar(unzip())==0)
    skip("No unzip tool available!")
  
  zipfile=file.path('testdata','sample.zip')
  r=zipinfo(zipfile)
  # only check fields which will not give trouble on machines with timezone
  # differences i.e. drop Date/Time columns
  chkfields=intersect(names(r), c("Length", "Method", "Size", "Ratio", "CRC.32", "Name"))
  
  baseline=structure(list(Length = c(0L, 25L, 0L, 0L), Method = structure(c(1L, 
                  1L, 1L, 1L), .Label = "Stored", class = "factor"), Size = c(0L, 
              25L, 0L, 0L), Ratio = c("0%", "0%", "0%", "0%"), Date = c("06-17-13", 
              "06-17-13", "06-17-13", "06-17-13"), Time = c("03:48", "03:48", 
              "03:48", "03:48"), CRC.32 = c("00000000", "5e61c365", "00000000", 
              "00000000"), Name = c("sample/", "sample/README", "sample/somedir/", 
              "sample/somedir/empty")), .Names = c("Length", "Method", "Size", 
          "Ratio", "Date", "Time", "CRC.32", "Name"), class = "data.frame", row.names = c(NA, 
          -4L))
  expect_equal(r[chkfields], baseline[chkfields])
  
  zipfile_with_spaces=file.path('testdata','sample-spaces.zip')
  r2=zipinfo(zipfile_with_spaces)
  
  baseline2=structure(list(Length = c(0L, 25L), Method = structure(c(1L, 
                  1L), .Label = "Stored", class = "factor"), Size = c(0L, 25L), 
          Ratio = c("0%", "0%"), Date = c("06-17-13", "06-17-13"), 
          Time = c("11:05", "03:48"), CRC.32 = c("00000000", "5e61c365"
          ), Name = c("sample/", "sample/file with spaces in its name"
          )), .Names = c("Length", "Method", "Size", "Ratio", "Date", 
          "Time", "CRC.32", "Name"), row.names = c(NA, -2L), class = "data.frame")
  expect_equal(r2[chkfields], baseline2[chkfields],
               label='parse zip file with paths containing spaces')
  
  zipfile_onefile=file.path('testdata','sample-onefile.zip')
  r3=zipinfo(zipfile_onefile)
  
  baseline3=structure(list(Length = 25L, Method = structure(1L, .Label = "Stored", class = "factor"), 
          Size = 25L, Ratio = "0%", Date = "06-17-13", Time = "03:48", 
          CRC.32 = "5e61c365", Name = "sample/README"), .Names = c("Length", 
          "Method", "Size", "Ratio", "Date", "Time", "CRC.32", "Name"), row.names = c(NA, 
          -1L), class = "data.frame")
  
  expect_equal(r3[chkfields], baseline3[chkfields],
               label='parse zip file containing only one entry')

})

test_that('can get test integrity of sample zip file',{
  if(nchar(unzip())==0)
    skip("No unzip tool available!")
  
  zipfile=file.path('testdata','sample.zip')
  notzipfile='test-ziputils.r'
  expect_true(zipok(zipfile))
  expect_false(zipok(notzipfile))
})
