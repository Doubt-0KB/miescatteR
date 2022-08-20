#' Calculate the intensity of scattered light
#'
#' @details
#' \deqn{I(\theta,\varphi) = cos^2\varphi I_\parallel (\theta) + sin^2\varphi I_\perp(\theta)}
#'
#' \deqn{
#' I_\parallel (\theta)
#'  = (\frac{\lambda}{2 \pi r})^2|\sum^\infty_{n=1} \frac{2n+1}{n(n+1)}(a_n \pi_n(cos\theta) + b_n \tau_n(cos\theta))|^2
#' }
#' \deqn{
#' I_\perp (\theta)
#'  = (\frac{\lambda}{2 \pi r})^2|\sum^\infty_{n=1} \frac{2n+1}{n(n+1)}(a_n \tau_n(cos\theta) + b_n \pi_n(cos\theta))|^2
#' }
#'
#' @param radius A numeric of the particle radius.
#' @param lambda A numeric of the light's wavelength.
#' @param np A numeric or complex of the refractive indix of the particle.
#' @param nm A numeric or complex of the refractive indix of the medium.
#' @param nang A integer of the angle division number.
#'
#' @importFrom dplyr
#'   mutate
#'   summarize
#'   group_by
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom data.table data.table
#'
#' @examples
#' miescatteR(0.26*10^-6, 0.55*10^-6, 1.33+(1i*10^-8), 1.0, 360)
#' miescatteR( 50*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33, 360)
#' miescatteR(100*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33, 360)
#' miescatteR(500*10^-9, 0.55*10^-6, 0.57+2.45i, 1.33, 360)
#'
#' @export
miescatteR <- function(radius, lambda, np, nm, nang=360){

  alpha <- 2*pi*radius/lambda*nm
  if(length(alpha)!=1) stop("\u8907\u6570\u306e\u5024\u304c\u5165\u529b\u3055\u308c\u3066\u3044\u307e\u3059")
  m <- np/nm

  qq <- 50 - 1
  pii <- rep(0, qq)
  tau <- rep(0, qq)
  n <- 1:qq

  an <- as.complex(coef_an(alpha, 1, m, qq))
  bn <- as.complex(coef_bn(alpha, 1, m, qq))

  data.table::data.table(
    theta = 1:nang - 1
  ) |>
    dplyr::mutate(
      costh = cospi(2.0/(nang-1) * !!as.name("theta"))
    ) |>
    dplyr::mutate(
      dat = purrr::map(
        !!as.name("costh"),
        function(x){
          p0 <- 0
          p1 <- 1
          #pii[0] <- 0
          #tau[0] <- 0
          pii[1] <- 1
          tau[1] <- x

          for(n in 2:qq){
            pn <- (2*n-1)/(n-1)*x*p1 - n/(n-1)*p0
            tn <- n*x*pn - (n+1)*p1
            p0 <- p1
            p1 <- pn
            pii[n] <- pn
            tau[n] <- tn
          }

          data.table::data.table(
            pn = pii,
            tn = tau
          )
        }
      )
    ) |>
    tidyr::unnest(
      !!as.name("dat")
    ) |>
    dplyr::group_by(
      !!as.name("theta")
    ) |>
    dplyr::summarize(
      I_prep = sum((2*n+1)/(n*(n+1)) * (an * !!as.name("pn") + bn * !!as.name("tn"))),
      I_parallel = sum((2*n+1)/(n*(n+1)) * (an * !!as.name("tn") + bn * !!as.name("pn")))
    )

}
