subset.rtconnect <-
function(x, sku=NULL, developer=NULL, title=NULL,
                             version=NULL, version.lte=NULL, version.gte=NULL,
                             product.type.identifier=NULL,
                             units=NULL, units.lte=NULL, units.gte=NULL,
                             developer.proceeds=NULL,
                             developer.proceeds.lte=NULL,
                             developer.proceeds.gte=NULL,
                             date=NULL, date.lte=NULL, date.gte=NULL,
                             customer.currency=NULL, country.code=NULL,
                             currency.of.proceeds=NULL, apple.identifier=NULL,
                             customer.price=NULL, promo.code=NULL,
                             parent.identifier=NULL, subscription=NULL,
                             period=NULL, ...) {
  rtc <- x
  class(rtc) <- "data.frame"

  # For R CMD check
  SKU <- NULL
  Developer <- NULL
  Title <- NULL
  Version <- NULL
  Product.Type.Identifier <- NULL
  Units <- NULL
  Developer.Proceeds <- NULL
  Begin.Date <- NULL
  End.Date <- NULL
  Customer.Currency <- NULL
  Country.Code <- NULL
  Currency.of.Proceeds <- NULL
  Apple.Identifier <- NULL
  Customer.Price <- NULL
  Promo.Code <- NULL
  Parent.Identifier <- NULL
  Subscription <- NULL
  Period <- NULL

  if (!(is.null(sku))) {
    rtc <- subset(rtc, SKU==sku)
  }

  if (!(is.null(developer))) {
    rtc <- subset(rtc, Developer==developer)
  }

  if (!(is.null(title))) {
    rtc <- subset(rtc, Title==title)
  }

  if (!(is.null(version))) {
    rtc <- subset(rtc, Version==version)
  }
  if (!(is.null(version.lte))) {
    rtc <- subset(rtc, Version<=version.lte)
  }
  if (!(is.null(version.gte))) {
    rtc <- subset(rtc, Version>=version.gte)
  }

  if (!(is.null(product.type.identifier))) {
    e <- NULL
    for (pti in product.type.identifier) {
      e <- rbind(e, subset(rtc, Product.Type.Identifier==pti))
    }
    rtc <- e
    e <- NULL
  }

  if (!(is.null(units))) {
    rtc <- subset(rtc, Units==units)
  }
  if (!(is.null(units.lte))) {
    rtc <- subset(rtc, Units<=units.lte)
  }
  if (!(is.null(units.gte))) {
    rtc <- subset(rtc, Units>=units.gte)
  }

  if (!(is.null(developer.proceeds))) {
    rtc <- subset(rtc, Developer.Proceeds==developer.proceeds)
  }
  if (!(is.null(developer.proceeds.lte))) {
    rtc <- subset(rtc, Developer.Proceeds<=developer.proceeds.lte)
  }
  if (!(is.null(developer.proceeds.gte))) {
    rtc <- subset(rtc, Developer.Proceeds>=developer.proceeds.gte)
  }

  if (!(is.null(date))) {
    rtc <- subset(rtc, Begin.Date<=date)
    rtc <- subset(rtc, End.Date>=date)
  }
  if (!(is.null(date.lte))) {
    rtc <- subset(rtc, End.Date<=date.lte)
  }
  if (!(is.null(date.gte))) {
    rtc <- subset(rtc, Begin.Date>=date.gte)
  }

  if (!(is.null(customer.currency))) {
    rtc <- subset(rtc, Customer.Currency==customer.currency)
  }

  if (!(is.null(country.code))) {
    rtc <- subset(rtc, Country.Code==country.code)
  }

  if (!(is.null(currency.of.proceeds))) {
    rtc <- subset(rtc, Currency.of.Proceeds==currency.of.proceeds)
  }

  if (!(is.null(apple.identifier))) {
    rtc <- subset(rtc, Apple.Identifier==apple.identifier)
  }

  if (!(is.null(customer.price))) {
    rtc <- subset(rtc, Customer.Price==customer.price)
  }

  if (!(is.null(promo.code))) {
    rtc <- subset(rtc, Promo.Code==promo.code)
  }

  if (!(is.null(parent.identifier))) {
    rtc <- subset(rtc, Parent.Identifier==parent.identifier)
  }

  if (!(is.null(subscription))) {
    rtc <- subset(rtc, Subscription==subscription)
  }

  if (!(is.null(period))) {
    rtc <- subset(rtc, Period==period)
  }

  class(rtc) <- c("rtconnect", "data.frame")

  return (rtc)
}
