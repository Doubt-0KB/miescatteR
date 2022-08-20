Csca <- function(radius, lambda, np, nm){

  m <- np/nm
  alpha <- 2 * pi * radius / lambda * nm
  nstop <- alpha + 4 * alpha^(1/3) + 2 + 15
  qq <- as.integer(nstop) - 1
  n <- 1:qq

  an <- as.complex(coef_an(alpha, 1, m, qq))
  bn <- as.complex(coef_bn(alpha, 1, m, qq))

  2 * pi * lambda^2 * sum( (2*n + 1) * (abs(an)^2 + abs(bn)^2) )

}



Cext <- function(radius, lambda, np, nm){

  m <- np/nm
  alpha <- 2 * pi * radius / lambda * nm
  nstop <- alpha + 4 * alpha^(1/3) + 2 + 15
  qq <- as.integer(nstop) - 1
  n <- 1:qq

  an <- as.complex(coef_an(alpha, 1, m, qq))
  bn <- as.complex(coef_bn(alpha, 1, m, qq))

  2 * pi * lambda^2 * sum( (2*n + 1) * Re(an + bn) )

}



#' Scattering cross-section
#'
#' Calculate the scattering cross-section.
#'
#' @inheritParams miescatteR
#'
#' @importFrom dplyr mutate
#' @importFrom purrr map2_dbl
#' @importFrom tidyr expand_grid
#'
#' @examples
#' cross_section_sca(50*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#' cross_section_sca(c(50, 100, 500)*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#'
#' @export
cross_section_sca <- function(radius, lambda, np, nm){

  tidyr::expand_grid(
    radius = radius,
    lambda = lambda
  ) |>
    dplyr::mutate(
      Csca = purrr::map2_dbl(!!as.name("radius"), !!as.name("lambda"), Csca, np=np, nm=nm)
    )

}



#' Extinction cross-section
#'
#' Calculate the extinction cross-section.
#'
#' @inheritParams miescatteR
#'
#' @importFrom dplyr mutate
#' @importFrom purrr map2_dbl
#' @importFrom tidyr expand_grid
#'
#' @examples
#' cross_section_ext(50*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#' cross_section_ext(c(50, 100, 500)*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#'
#' @export
cross_section_ext <- function(radius, lambda, np, nm){

  tidyr::expand_grid(
    radius = radius,
    lambda = lambda
  ) |>
    dplyr::mutate(
      Cext = purrr::map2_dbl(!!as.name("radius"), !!as.name("lambda"), Cext, np=np, nm=nm)
    )

}



#' Scattering and extinction cross-sections
#'
#' @description
#' Calculate scattering and extinction cross-section.
#' * Csca: scattering cross-section
#' * Cext: extinction cross-section
#'
#' @inheritParams miescatteR
#'
#' @importFrom dplyr inner_join
#'
#' @examples
#' cross_section(50*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#' cross_section(c(50, 100, 500)*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#'
#' @export
cross_section <- function(radius, lambda, np, nm){

  dplyr::inner_join(
    cross_section_sca(radius, lambda, np, nm),
    cross_section_ext(radius, lambda, np, nm),
    by = c("radius", "lambda")
  )

}
