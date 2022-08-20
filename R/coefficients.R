#' The coeffiicients of scattered fields
#'
#' @details
#' \eqn{a_n} and \eqn{b_n} are the coefficients
#' used for calculating scattered light intensity and the Mie total scattering cross section.
#'
#' \deqn{
#' a_n =
#' \frac{m\psi_n(m\alpha){\psi_n}'(\alpha)-\psi_n(\alpha){\psi_n}'(m\alpha)}
#'   {m\psi_n(m\alpha){\xi_n}'(\alpha)-\xi_n(\alpha){\psi_n}'(m\alpha)}
#' }
#' \deqn{
#' b_n =
#' \frac{\psi_n(m\alpha){\psi_n}'(\alpha)-m\psi_n(\alpha){\psi_n}'(m\alpha)}
#'   {\psi_n(m\alpha){\xi_n}'(\alpha)-m\xi_n(\alpha){\psi_n}'(m\alpha)}
#' }
#'
#' @inheritParams spherical_Bessel
#' @param m A numeric or complex.
#'     This is the ratio of the refractive index of the particle to that of the surrounding medium.
#'
#' @examples
#' coef_an(2*pi*0.05/0.55, 1, 1.33, 10)
#' coef_bn(2*pi*0.05/0.55, 1, 1.33, 10)
#'
#' @export
coef_an <- function(z, nu, m, nSeq = 1){
  (m * psinz(m*z, nu, nSeq) * psindz(z, nu, nSeq) - psinz(z, nu, nSeq) * psindz(m*z, nu, nSeq)) /
    (m * psinz(m*z, nu, nSeq) * xindz(z, nu, nSeq) - xinz(z, nu, nSeq) * psindz(m*z, nu, nSeq))
}

#' @rdname coef_an
#' @export
coef_bn <- function(z, nu, m, nSeq = 1){
  (psinz(m*z, nu, nSeq) * psindz(z, nu, nSeq) - m * psinz(z, nu, nSeq) * psindz(m*z, nu, nSeq)) /
    (psinz(m*z, nu, nSeq) * xindz(z, nu, nSeq) - m * xinz(z, nu, nSeq) * psindz(m*z, nu, nSeq))
}
