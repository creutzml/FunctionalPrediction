#' Meanfunction (polynomial)
#'
#' @param x function argument
#' @param params Parameters: params=c(shift,scale).
#' @examples
#' curve(meanf_poly(x, c(0, 2)),
#'   from = 0, to = 1,
#'   main = "Meanfct Poly", ylab = "", xlab = ""
#' )
#' curve(meanf_poly(x, c(0, 1)),
#'   from = 0, to = 1,
#'   lty = 2, add = TRUE
#' )
#' @export
meanf_poly <- function(x, params = c(0, 1)) {
  f <- params[1] + params[2] * (10 * x^3 - 15 * x^4 + 6 * x^5)
  names(f) <- x
  return(f)
}


#' Meanfunction (polynomial with simple shifting)
#'
#' @param x function argument
#' @param delta shifting parameter
#' @examples
#' curve(meanf_shift(x, 0), from = 0, to = 1, lty = 2)
#' curve(meanf_shift(x, .1), from = 0, to = 1, lty = 1, add = TRUE)
#' @export
meanf_shift <- function(x, delta = 0) {
  meanf_poly(x, c(delta, 1))
}


#' Meanfunction (polynomial with simple scaling)
#'
#' @param x function argument
#' @param delta scaling parameter
#' @examples
#' curve(meanf_scale(x, 1), from = 0, to = 1)
#' curve(meanf_scale(x, 0), from = 0, to = 1, lty = 2, add = TRUE)
#' @export
meanf_scale <- function(x, delta = 0) {
  meanf_poly(x, c(0, 1 + delta))
}


#' Meanfunction with local peak (polynomial with simple local peak)
#'
#' @param x function argument
#' @param delta shifting parameter
#' @examples
#' curve(meanf_localshift(x, 0), from = 0, to = 1, lty = 2)
#' curve(meanf_localshift(x, .1), from = 0, to = 1, lty = 1, add = TRUE)
#' @export
meanf_localshift <- function(x, delta = 0) {
  meanf_peak(x, c(0, 1, 1 + delta, 10, 1))
}


#' Meanfunction with local rectangles (polynomial with simple local rectangles)
#'
#' @param x function argument
#' @param delta hight of rectangle
#' @param rect location of the rectangle
#' @examples
#' curve(meanf_rect(x, delta = 1 / 5),
#'   from = 0, to = 1,
#'   main = "Meanfct Rect", ylab = "", xlab = ""
#' )
#' curve(meanf_rect(x, delta = 1 / 10),
#'   from = 0, to = 1,
#'   lty = 2, add = TRUE
#' )
#' @export
meanf_rect <- function(x, delta = 1 / 5, rect = c(0, 0.125)) {
  tmp <- numeric(length(x))
  tmp[rect[1] < x & x < rect[2]] <- delta
  result <- c(tmp + meanf_poly(x))
  return(result)
}


#' Meanfunction with local ellipse
#'
#' @param x function argument
#' @param params parameters params=c(start, end, height)
#' @examples
#' curve(meanf_ellipse(x, c(0, 0.25, 0.1)),
#'   from = 0, to = 1,
#'   main = "Meanfct Ellipse", ylab = "", xlab = ""
#' )
#' curve(meanf_ellipse(x, c(0.50, 0.75, 0.1)),
#'   from = 0, to = 1,
#'   lty = 2, add = TRUE
#' )
#' @export
meanf_ellipse <- function(x, params = c(0, 1 / 4, 1 / 5)) {
  ##
  a <- (params[2] - params[1]) / 2
  b <- params[3]
  x_s <- x - a - params[1]
  if (-a <= x_s & x_s <= a) {
    y <- sqrt(b^2 * (1 - x_s^2 / a^2))
  } else {
    y <- 0
  }
  ##
  return(y)
}
meanf_ellipse <- Vectorize(meanf_ellipse, vectorize.args = "x")



#' Meanfunction with local bump
#'
#' @param x function argument
#' @param quarter 1,2,3, or 4 to select domain-quarter.
#' @param height height of the bump function
#' @examples
#' curve(meanf_bump(x, 1),
#'   from = 0, to = 1,
#'   main = "Meanfct Ellipse", ylab = "", xlab = ""
#' )
#' curve(meanf_bump(x, 2),
#'   from = 0, to = 1,
#'   lty = 2, add = TRUE
#' )
#' @export
meanf_bump <- function(x, quarter = 1, height = .1) {
  ##
  x_s <- (x - (quarter - 1) * 0.25 - 0.125) / 0.125
  if (-1 <= x_s & x_s <= 1) {
    y <- height * exp(1 - 1 / (1 - x_s^2))
  } else {
    y <- 0
  }
  ##
  return(y)
}
meanf_bump <- Vectorize(meanf_bump, vectorize.args = "x")




#' Meanfunction with local peak
#'
#' @param x function argument
#' @param params Parameters
#' @examples
#' curve(meanf_peak(x, c(0, 1, 2, 16, 1)),
#'   from = 0, to = 1,
#'   main = "Meanfct Peak", ylab = "", xlab = ""
#' )
#' curve(meanf_peak(x, c(-1, 1, 2, 16, 1)),
#'   from = 0, to = 1,
#'   lty = 2, add = TRUE
#' )
#' curve(meanf_peak(x, c(0, 2, 2, 16, 1)),
#'   from = 0, to = 1,
#'   lty = 2, add = TRUE
#' )
#' curve(meanf_peak(x, c(0, 1, 1 + 1, 10, 1)),
#'   from = 0, to = 1,
#'   lty = 2, add = TRUE
#' )
#' @export
meanf_peak <- function(x, params = c(0, 1, 2, 16, 1)) {
  peak <- params[3] - params[4] * abs(x - 0.5)^(params[5])
  peak[peak < 0] <- 0 # params=c(shift,scale,peak-top,peak-sharpness,peak-power )
  f <- params[2] * (params[1] + peak)
  names(f) <- x
  return(f)
}
