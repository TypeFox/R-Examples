gedestp <-
function(x, method = c("inverse", "direct"))
{
    x <- as.vector(na.omit(x))
    p <- 2 #starting value
    mu <- 0
    if (!is.numeric(x))
        stop(" Non-numeric argument to mathematical function")
    method <- match.arg(method)
    ssp <- sum(abs(x - mu)^p)/length(x)
    sp <- ssp^(1/p)
    xstand <- (x - mu)/sp
    sa <- sum(abs(xstand))
    sb <- sum(xstand * xstand)
    vi <- sqrt(length(x) * sb)/sa
    vi <- vi + ((vi - 1)/length(x)) * 5
    if (method == "inverse") {
        if (vi < 1.1547005)
            pz <- 11.5
        else {
            zi <- (1.4142135 - vi)/0.259513
            if (zi < 0)
                pz <- 1
            else {
                if (zi < 0.6200052) {
                  tz <- zi/0.6200052
                  yy <- ((0.4738581 * tz - 0.4966873) * tz +
                    1.0532646) * tz + 1.2246159
                  pz <- 1 + tz^yy
                }
                else {
                  if (zi < 0.7914632) {
                    tz <- (zi - 0.6200052)/0.1714592
                    yy <- ((0.5246979 * tz - 0.8167733) * tz +
                      0.8805483) * tz + 1.0859246
                    pz <- 2 + tz^yy
                  }
                  else {
                    if (zi < 0.8670333) {
                      tz <- (zi - 0.7914632)/0.0755701
                      yy <- ((0.0743092 * tz - 0.1269859) * tz +
                        0.3588207) * tz + 1.1227837
                      pz <- 3 + tz^yy
                    }
                    else {
                      if (zi < 0.9072536) {
                        tz <- (zi - 0.8670333)/0.0402203
                        yy <- ((0.1097723 * tz - 0.2127039) *
                          tz + 0.3529203) * tz + 1.0761256
                        pz <- 4 + tz^yy
                      }
                      else {
                        if (zi < 0.9314555) {
                          tz <- (zi - 0.9072536)/0.0242019
                          yy <- ((0.0955441 * tz - 0.1891569) *
                            tz + 0.2961275) * tz + 1.0631784
                          pz <- 5 + tz^yy
                        }
                        else {
                          if (zi < 0.9472072) {
                            tz <- (zi - 0.9314555)/0.0157518
                            yy <- ((0.0862627 * tz - 0.1725326) *
                              tz + 0.256885) * tz + 1.0540746
                            pz <- 6 + tz^yy
                          }
                          else {
                            if (zi < 0.9580557) {
                              tz <- (zi - 0.9472072)/0.0108484
                              yy <- ((0.078785 * tz - 0.1581388) *
                                tz + 0.227011) * tz + 1.04735
                              pz <- 7 + tz^yy
                            }
                            else {
                              if (zi < 0.9658545) {
                                tz <- (zi - 0.9580557)/0.0077988
                                yy <- ((0.0663921 * tz - 0.1380841) *
                                  tz + 0.2010053) * tz + 1.0422984
                                pz <- 8 + tz^yy
                              }
                              else {
                                if (zi < 0.9716534) {
                                  tz <- (zi - 0.9658545)/0.0057989
                                  yy <- ((0.0557199 * tz - 0.1184033) *
                                    tz + 0.178176) * tz + 1.038582
                                  pz <- 9 + tz^yy
                                }
                                else pz <- 10 + (zi - 0.9716534)/0.0283466
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
            }
        }
    }
    if (method == "direct") {
        fvi <- function(p) (vi - sqrt(gamma(1/p) * gamma(3/p))/gamma(2/p)) *
            (vi - sqrt(gamma(1/p) * gamma(3/p))/gamma(2/p))
        pz <- optim(p, fvi, method = "BFGS")$par
        if (pz < 1)
            pz <- 1
        if (pz > 10)
            pz <- 11.5
    }
    pp <- pz
    pp
}
