#' Use simulations to compare the observed distribution with the modelled
#' distribution
#'
setGeneric(
  name = "fast_distribution_check",
  def = function(object, nsim = 1000) {
    standardGeneric("fast_distribution_check") # nocov
  }
)

#' library(INLA)
#' set.seed(20181202)
#' model <- inla(
#'   poisson ~ 1,
#'   family = "poisson",
#'   data = data.frame(
#'     poisson = rpois(20, lambda = 10),
#'     base = 1
#'   ),
#'   control.predictor = list(compute = TRUE)
#' )
#' fast_distribution_check(model)
setMethod(
  f = "fast_distribution_check",
  signature = signature(object = "inla"),
  definition = function(object, nsim = 1000) {
    assertthat::assert_that(assertthat::is.count(nsim))
    
    if (length(object$.args$family) > 1) {
      stop("Only single responses are handled")
    }
    
   # observed <- get_observed(object)
    observed=Cs.Copies
    mu <- fitted(object)[!is.na(observed)]
    observed <- observed[!is.na(observed)]
    n_mu <- length(mu)
    n_sampled <- switch(
      object$.args$family,
      poisson = {
        data.frame(
          run = rep(seq_len(nsim), each = n_mu),
          x = rpois(n = n_mu * nsim, lambda = mu)
        ) %>%
          count(.data$run, .data$x)
      },
      nbinomial = {
        relevant <- grep("overdispersion", rownames(object$summary.hyperpar))
        size <- object$summary.hyperpar[relevant, "mean"]
        data.frame(
          run = rep(seq_len(nsim), each = n_mu),
          x = rnbinom(n = n_mu * nsim, mu = mu, size = size)
        ) %>%
          count(.data$run, .data$x)
      },
      gpoisson = {
        relevant <- grep("Overdispersion", rownames(object$summary.hyperpar))
        phi <- object$summary.hyperpar[relevant, "mean"]
        data.frame(
          run = rep(seq_len(nsim), n_mu),
          x = as.vector(rgpoisson(n = nsim, mu = mu, phi = phi))
        ) %>%
          count(.data$run, .data$x)
      },
      zeroinflatedpoisson1 = {
        relevant <- grep("zero-probability", rownames(object$summary.hyperpar))
        zero <- object$summary.hyperpar[relevant, "mean"]
        data.frame(
          run = rep(seq_len(nsim), each = n_mu),
          x = rpois(n = n_mu * nsim, lambda = mu) *
            rbinom(n = n_mu * nsim, size = 1, prob = 1 - zero)
        ) %>%
          count(.data$run, .data$x)
      },
      zeroinflatedpoisson0 = {
        relevant <- grep("zero-probability", rownames(object$summary.hyperpar))
        if (length(relevant) == 1) {
          zero <- object$summary.hyperpar[relevant, "mean"]
        } else {
          zero <- object$all.hyper$family[[1]]$hyper$theta$from.theta(
            object$all.hyper$family[[1]]$hyper$theta$initial
          )
        }
        data.frame(
          run = rep(seq_len(nsim), each = n_mu),
          x = rtpois(n = n_mu * nsim, lambda = mu) *
            rbinom(n = n_mu * nsim, size = 1, prob = 1 - zero)
        ) %>%
          count(.data$run, .data$x)
      },
      zeroinflatednbinomial1 = {
        relevant <- grep("zero-probability", rownames(object$summary.hyperpar))
        zero <- object$summary.hyperpar[relevant, "mean"]
        relevant <- grep(
          "size for nbinomial",
          rownames(object$summary.hyperpar)
        )
        size <- object$summary.hyperpar[relevant, "mean"]
        data.frame(
          run = rep(seq_len(nsim), each = n_mu),
          x = rnbinom(n = n_mu * nsim, mu = mu, size = size) *
            rbinom(n = n_mu * nsim, size = 1, prob = 1 - zero)
        ) %>%
          count(.data$run, .data$x)
      },
      stop(object$.args$family, " is not yet handled")
    )
    
    data.frame(x = observed) %>%
      count(.data$x) -> n_observed
    n_count <- unique(c(n_observed$x, n_sampled$x))
    n_sampled %>%
      complete(run = .data$run, x = n_count, fill = list(n = 0)) %>%
      group_by(.data$run) %>%
      arrange(.data$x) %>%
      mutate(ecdf = cumsum(.data$n) / sum(.data$n)) %>%
      group_by(.data$x) %>%
      summarise(
        median = quantile(.data$ecdf, probs = 0.5),
        lcl = quantile(.data$ecdf, probs = 0.025),
        ucl = quantile(.data$ecdf, probs = 0.975)
      ) %>%
      inner_join(
        n_observed %>%
          complete(x = n_count, fill = list(n = 0)) %>%
          arrange(.data$x) %>%
          mutate(ecdf = cumsum(.data$n) / sum(.data$n))
        ,
        by = "x"
      ) -> ecdf
    class(ecdf) <- c("distribution_check", class(ecdf))
    return(ecdf)
  }
)

#' @rdname fast_distribution_check
#' @importFrom methods setMethod new
#' @importFrom purrr map map2
setMethod(
  f = "fast_distribution_check",
  signature = signature(object = "list"),
  definition = function(object, nsim = 1000) {
    ecdf <- map(object, fast_distribution_check)
    if (is.null(names(object))) {
      ecdf <- map2(ecdf, seq_along(object), ~mutate(.x, model = .y))
    } else {
      ecdf <- map2(ecdf, names(object), ~mutate(.x, model = .y))
    }
    ecdf <- bind_rows(ecdf)
    class(ecdf) <- c("distribution_check", class(ecdf))
    return(ecdf)
  }
)
number
dgpoisson <- function(y, mu, phi) {
  assert_that(
    is.integer(y),
    all(y >= 0)
  )
  assert_that(
    is.numeric(mu),
    all(mu > 0)
  )
  assert_that(
    is.number(phi),
    phi > 0
  )
  
  a <- outer(phi * y, mu, "+")
  b <- 1 + phi
  exp(
    matrix(log(mu), nrow = length(y), ncol = length(mu), byrow = TRUE) +
      (y - 1) * log(a) - y * log(b) - lfactorial(y) - a / b
  )
}

#' @noRd
#' @inheritParams dgpoisson
#' @param n the number of simulated values
#' @importFrom assertthat assert_that is.number is.count
rgpoisson <- function(n, mu, phi) {
  assert_that(is.count(n))
  assert_that(
    is.numeric(mu),
    all(mu > 0)
  )
  assert_that(
    is.number(phi),
    phi > 0
  )
  
  s <- sqrt(max(mu) * (1 + phi) ^ 2)
  low <- as.integer(max(0, min(mu) - 20 * s))
  high <- as.integer(max(mu) + 20 * s)
  prob <- dgpoisson(y = low:high, mu, phi)
  y <- apply(prob, 2, sample, x = low:high, replace = TRUE, size = n)
  return(y)
}

#' @importFrom assertthat assert_that is.count
rtpois <- function(n, lambda) {
 assertthat::assert_that(assertthat::is.count(n))
  assertthat::assert_that(inherits(lambda, "numeric"))
  assertthat::assert_that(all(lambda >= 0))
  if (length(lambda) < n) {
    lambda <- head(rep(lambda, ceiling(n / length(lambda))), n)
  }
  y <- rpois(n = n, lambda = lambda)
  while (any(y < 1)) {
    y[y < 1] <- rpois(sum(y < 1), lambda = lambda[y < 1])
  }
  return(y)
}