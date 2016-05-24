
context('preprocess_response')

XML_FILES <- c(
    'GetChart.xml',  'GetComps.xml',  'GetDeepComps.xml',  'GetDeepSearchResults.xml',
    'GetMonthlyPayments.xml',  'GetRateSummary.xml',  'GetSearchResults.xml',
    'GetUpdatedPropertyDetails.xml',  'GetZestimate.xml'
)

for (xml_file in XML_FILES) {
    test_that(sprintf("%s response is coerced to an R-friendly structure", xml_file), {
        response <- readLines(system.file('xml_response_examples', xml_file, package = 'ZillowR'))

        expect_true(is.list(preprocess_response(response)))
        expect_identical(length(preprocess_response(response)), 3L)
        expect_identical(names(preprocess_response(response)), c('request', 'message', 'response'))
        if (xml_file == 'GetRateSummary.xml') {
            expect_true(methods::is(preprocess_response(response)$request, 'NULL'))
        } else {
            expect_true(methods::is(preprocess_response(response)$request, 'list'))
        }
        expect_true(methods::is(preprocess_response(response)$message, 'list'))
        expect_true(methods::is(preprocess_response(response)$response, 'XMLNode'))
    })
}

rm(xml_file, XML_FILES)
