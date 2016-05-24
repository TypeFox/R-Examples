## ----sample--------------------------------------------------------------
library(curlconverter)
library(jsonlite)

curl_command <- "curl 'http://financials.morningstar.com/ajax/ReportProcess4HtmlAjax.html?&t=XNAS:MSFT&region=usa&culture=en-US&cur=&reportType=is&period=12&dataType=A&order=asc&columnYear=5&curYearPart=1st5year&rounding=3&view=raw&r=973302&callback=jsonp1454021128757&_=1454021129337' -H 'Cookie: JSESSIONID=5E43C98903E865D72AA3C2DCEF317848; sfhabit=asc%7Craw%7C3%7C12%7CA%7C5%7Cv0.14; ScrollY=0' -H 'DNT: 1' -H 'Accept-Encoding: gzip, deflate, sdch' -H 'Accept-Language: en-US,en;q=0.8' -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/47.0.2526.111 Safari/537.36' -H 'Accept: text/javascript, application/javascript, */*' -H 'Referer: http://financials.morningstar.com/income-statement/is.html?t=MSFT&region=usa&culture=en-US' -H 'X-Requested-With: XMLHttpRequest' -H 'Connection: keep-alive' -H 'Cache-Control: max-age=0' --compressed"

toJSON((cURL <- straighten(curl_command)), pretty=TRUE)

## ----req1, echo=FALSE----------------------------------------------------
req <- make_req(cURL, add_clip=FALSE)


## ----req2, eval=FALSE----------------------------------------------------
#  req <- make_req(cURL)

## ----funct, eval=FALSE---------------------------------------------------
#  httr::VERB(verb = "GET", url = "http://financials.morningstar.com/ajax/ReportProcess4HtmlAjax.html?&t=XNAS:MSFT&region=usa&culture=en-US&cur=&reportType=is&period=12&dataType=A&order=asc&columnYear=5&curYearPart=1st5year&rounding=3&view=raw&r=973302&callback=jsonp1454021128757&_=1454021129337",
#             httr::add_headers(DNT = "1",
#                               `Accept-Encoding` = "gzip, deflate, sdch",
#                               `Accept-Language` = "en-US,en;q=0.8",
#                               `User-Agent` = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_3) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/47.0.2526.111 Safari/537.36",
#                               Accept = "text/javascript, application/javascript, */*",
#                               Referer = "http://financials.morningstar.com/income-statement/is.html?t=MSFT&region=usa&culture=en-US",
#                               `X-Requested-With` = "XMLHttpRequest",
#                               Connection = "keep-alive",
#                               `Cache-Control` = "max-age=0"),
#             httr::set_cookies(JSESSIONID = "5E43C98903E865D72AA3C2DCEF317848",
#                               sfhabit = "asc%7Craw%7C3%7C12%7CA%7C5%7Cv0.14",
#                               ScrollY = "0"))

## ----calling_req---------------------------------------------------------
res <- req[[1]]()

## ----muck----------------------------------------------------------------
library(V8)

ctx <- v8()

ctx$eval(sub("\\)$", ";", sub("^jsonp1454021128757\\(", "var dat=", httr::content(res, as="text"))))
names(ctx$get("dat"))

## ----pipe, eval=FALSE----------------------------------------------------
#  straighten() %>%
#    make_req() -> req

## ----bracket, eval=FALSE-------------------------------------------------
#  req <- make_req(straighten())[[1]]

