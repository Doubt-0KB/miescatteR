#' @name spherical_Bessel
#'
#' @title
#' Spherical Bessel and Hankel functions.
#'
#' @description
#' Extend Bessel functions, modified Bessel functions and Hankel function to spherical functions.
#'
#' @details
#' The differential equation
#'
#' \deqn{x^2 \frac{d^2y}{dx^2} + 2x \frac{dy}{dx} + (x^2 - n(n+1))y = 0}
#'
#' has two linearly independent solutions which are called the spherical Bessel functions.
#'
#' \deqn{j_n(x) = \sqrt{\frac{\pi}{2x}} J_{n+\frac{1}{2}}(x)}
#' \deqn{y_n(x) = \sqrt{\frac{\pi}{2x}} Y_{n+\frac{1}{2}}(x)}
#'
#' These are the spherical Hankel function.
#'
#' \deqn{h^{(1)}_n(x) = j_n(x) + iy_n(x)}
#' \deqn{h^{(2)}_n(x) = j_n(x) - iy_n(x)}
#'
#' @param z A complex or numeric vector.
#' @param nu A numeric.
#' @param expon.scaled Logical. The result should be scaled by an exponential factor,
#'   typically to avoid under- or over-flow.
#' @param nSeq A positive integer. if >1, the result has \eqn{j_{nu}}, \eqn{j_{nu+1}}, ..., \eqn{j_{nu+nSeq-1}}
#' @param verbose A integer (default is 0). Indicating the level of verbosity notably from C code.
#' @param m 1 or 2. Indicate the kind of Hankel function.
#'
#' @examples
#' spherical_BesselJ(3+5i, 1)
#' spherical_BesselJ(3+5i, 1, nSeq = 4)
#' spherical_BesselJ(c(3+5i, 1+2i, 7+6i), 1)
#' spherical_BesselJ(c(3+5i, 1+2i, 7+6i), 1, nSeq = 4)
#'
NULL



#' @rdname spherical_Bessel
#' @importFrom Bessel BesselJ
#' @export
spherical_BesselJ <- function(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0){
  sqrt(pi/2/z) * Bessel::BesselJ(z, nu+1/2, expon.scaled, nSeq, verbose)
}



#' @rdname spherical_Bessel
#' @importFrom Bessel BesselY
#' @export
spherical_BesselY <- function(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0){
  sqrt(pi/2/z) * Bessel::BesselY(z, nu+1/2, expon.scaled, nSeq, verbose)
}



#' @rdname spherical_Bessel
#' @importFrom Bessel BesselI
#' @export
spherical_BesselI <- function(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0){
  sqrt(pi/2/z) * Bessel::BesselI(z, nu+1/2, expon.scaled, nSeq, verbose)
}



#' @rdname spherical_Bessel
#' @importFrom Bessel BesselK
#' @export
spherical_BesselK <- function(z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0){
  sqrt(2/pi/z) * Bessel::BesselK(z, nu+1/2, expon.scaled, nSeq, verbose)
}



#' @rdname spherical_Bessel
#' @importFrom Bessel BesselH
#' @export
spherical_BesselH <- function(m, z, nu, expon.scaled = FALSE, nSeq = 1, verbose = 0){
  sqrt(pi/2/z) * Bessel::BesselH(m, z, nu+1/2, expon.scaled, nSeq, verbose)
}
