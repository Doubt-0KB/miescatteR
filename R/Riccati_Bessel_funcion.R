#' Riccati-Bessel Funcions
#'
#' The special functions mostly used in calculating mie scattering.
#'
#' @details
#' \eqn{\psi_n(z)} and \eqn{\xi_n(z)}, called Riccati-Bessel functions, are the special functions which slightly differ from spherical Bessel functions.
#'
#' \deqn{\psi_n(z) = zj_n(z)}
#' \deqn{\xi_n(z) = zh^{(1)}_n(z)}
#'
#' These functions satisfy the differential equation;
#' \deqn{x^2\frac{d^2y}{dx^2} + (x^2 - n(n+1))y = 0}
#'
#' \eqn{\psi'_n(z)} and \eqn{\xi'_n(z)} are written like these by using recurrence formula;
#'
#' \deqn{\psi'_n(z) = j_{n-1}(z) - \frac{n}{z}j_n(z)}
#' \deqn{\xi'_n(z) = h^{(1)}_{n-1}(z) - \frac{n}{z}h^{(1)}_n(z)}
#'
#' @examples
#' psinz(1, 1, 10)
#' xinz(1, 1, 10)
#' psindz(1, 1, 10)
#' xindz(1, 1, 10)
#'
#' @inheritParams spherical_Bessel
#'
#' @export
psinz <- function(z, nu, nSeq = 1){
  z * spherical_BesselJ(z, nu, nSeq = nSeq)
}

#' @rdname psinz
#' @export
xinz <- function(z, nu, nSeq = 1){
  z * spherical_BesselH(m=1, z, nu, nSeq = nSeq)
}

#' @rdname psinz
#' @export
psindz <- function(z, nu, nSeq = 1){
  psinz(z, nu-1, nSeq) - t(sapply(z, \(x) nu:nSeq/x)) * psinz(z, nu, nSeq)
}

#' @rdname psinz
#' @export
xindz <- function(z, nu, nSeq = 1){
  xinz(z, nu-1, nSeq) - t(sapply(z, \(x) nu:nSeq/x)) * xinz(z, nu, nSeq)
}
