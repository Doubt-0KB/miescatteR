% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coefficients.R
\name{coef_an}
\alias{coef_an}
\alias{coef_bn}
\title{The coeffiicients of scattered fields}
\usage{
coef_an(z, nu, m, nSeq = 1)

coef_bn(z, nu, m, nSeq = 1)
}
\arguments{
\item{z}{A complex or numeric vector.}

\item{nu}{A numeric.}

\item{m}{A numeric or complex.
This is the ratio of the refractive index of the particle to that of the surrounding medium.}

\item{nSeq}{A positive integer. if >1, the result has \eqn{j_{nu}}, \eqn{j_{nu+1}}, ..., \eqn{j_{nu+nSeq-1}}}
}
\description{
The coeffiicients of scattered fields
}
\details{
\eqn{a_n} and \eqn{b_n} are the coefficients
used for calculating scattered light intensity and the Mie total scattering cross section.

\deqn{
a_n =
\frac{m\psi_n(m\alpha){\psi_n}'(\alpha)-\psi_n(\alpha){\psi_n}'(m\alpha)}
  {m\psi_n(m\alpha){\xi_n}'(\alpha)-\xi_n(\alpha){\psi_n}'(m\alpha)}
}
\deqn{
b_n =
\frac{\psi_n(m\alpha){\psi_n}'(\alpha)-m\psi_n(\alpha){\psi_n}'(m\alpha)}
  {\psi_n(m\alpha){\xi_n}'(\alpha)-m\xi_n(\alpha){\psi_n}'(m\alpha)}
}
}
\examples{
coef_an(2*pi*0.05/0.55, 1, 1.33, 10)
coef_bn(2*pi*0.05/0.55, 1, 1.33, 10)

}
