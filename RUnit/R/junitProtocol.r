##  RUnit : A unit test framework for the R programming language
##  Copyright (C) 2003-2009  Thomas Koenig, Matthias Burger, Klaus Juenemann
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; version 2 of the License.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


printJUnitProtocol <- function(testData,
                              fileName = "") {

  ##@bdescr
  ##  Report generator
  ##  Extracts the log information stored in the 'RUnitTestData' test run object
  ##  and generates a JUnit-style formated XML output.
  ##@edescr
  ##
  ##@in  testData            : [RUnitTestData] S3 class object
  ##@in  fileName            : [character] string, full path + file name to be written to
  ##@ret                     : [logical] TRUE if execution completed without error
  ##
  ##@codestatus : testing

  ##  preconditions
  if (!is(testData, "RUnitTestData")) {
    stop("Argument 'testData' must be of class 'RUnitTestData'.")
  }

  if (!is.character(fileName)) {
    stop("Argument 'fileName' has to be of type character.")
  }
  if (length(fileName) != 1) {
    stop("Argument 'fileName' must contain exactly one element.")
  }

  errInfo <- getErrors(testData)
  # Create entry fro all test suites
  xml.testsuites <- XML::newXMLNode("testsuites", 
                             attrs = c(
                               errors=errInfo$nErr, 
                               failures=errInfo$nFail, 
                               tests=errInfo$nTestFunc)
  )

  for (tsName in names(testData)) {
    # Create entry for test suite
    xml.testsuite <- XML::newXMLNode("testsuite", 
                                attrs = c(
                                  errors = testData[[tsName]]$nErr,
                                  failures = testData[[tsName]]$nFail,
                                  name = tsName,
                                  tests = testData[[tsName]]$nTestFunc
                                ))
    XML::addChildren(xml.testsuites, kids=c(xml.testsuite))    

    if (testData[[tsName]]$nErr + testData[[tsName]]$nFail >= 0) {
      srcFileRes <- testData[[tsName]][["sourceFileResults"]]
      for (i in seq_along(srcFileRes)) {
        testFuncNames <- names(srcFileRes[[i]])
        for (j in seq_along(testFuncNames)) {
          funcList <- srcFileRes[[i]][[testFuncNames[j]]]
          # Each tested function gets a testcase
          xml.testcase <- XML::newXMLNode("testcase", attrs=c(name=testFuncNames[j], time=funcList$time[['elapsed']]))
          XML::addChildren(xml.testsuite, kids=c(xml.testcase))
          if (funcList$kind == "success") {
          } else if (funcList$kind == "error") {
            xml.error <- XML::newXMLNode("error", attrs=c(
              "message"=funcList$msg,
              "type"="ERROR"))
            XML::addChildren(xml.testcase, kids=c(xml.error))
          }
          else if (funcList$kind == "failure") {
            xml.error <- XML::newXMLNode("failure", attrs=c(
              "message"=funcList$msg,
              "type"="FAILURE"))
            XML::addChildren(xml.testcase, kids=c(xml.error))          
          }          
          else if (funcList$kind == "deactivated") {
            xml.skipped <- XML::newXMLNode("skipped")
            XML::addChildren(xml.testcase, kids=c(xml.skipped))                      
          }
        }
      }
    }      
  }
  xml <- XML::saveXML(xml.testsuites)
  if(fileName=="") {
    write(xml, stdout())
  } else {
    dir.create(dirname(fileName), showWarnings=FALSE, recursive=TRUE)
    fileConn <- file(fileName)
    write(xml, fileConn)
    close(fileConn)
  }
  return(invisible(TRUE))
}
