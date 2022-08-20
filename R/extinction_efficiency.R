Qext <- function(radius, lambda, np, nm){

  Cext(radius, lambda, np, nm) / (pi * radius^2)

}

#' Extinction efficiency
#'
#' Calculate the extinction efficiency.
#'
#' @inheritParams miescatteR
#'
#' @importFrom dplyr mutate
#' @importFrom purrr map2_dbl
#' @importFrom tidyr expand_grid
#'
#' @examples
#' ext_efficiency(50*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#' ext_efficiency(c(50, 100, 500)*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33)
#'
#' @export
ext_efficiency <- function(radius, lambda, np, nm){

  tidyr::expand_grid(
    radius = radius,
    lambda = lambda
  ) |>
    dplyr::mutate(
      alpha = purrr::map2_dbl(!!as.name("radius"), !!as.name("lambda"), \(x,y) 2*pi*x/y*nm),
      Qext = purrr::map2_dbl(!!as.name("radius"), !!as.name("lambda"), Qext, np=np, nm=nm)
    )

}
