angle <-
function(x, y) {
        if (x > 0) {
            atan(y/x)
        }
        else {
            if (x < 0 & y != 0) {
                atan(y/x) + sign(y) * pi
            }
            else {
                if (x < 0 & y == 0) {
                  pi
                }
                else {
                  if (y != 0) {
                    (sign(y) * pi)/2
                  }
                  else {
                    NA
                  }
                }
            }
        }
    }
