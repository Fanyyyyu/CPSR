####Data Generating Process############

dgp = function(T, P, p1, q, meanz, corZ, vare, eta_value, alpha1, alpha2, dgp_order){

  library(MASS)
  diag(corz) <- 1
  z <- array(1, dim = c(T, P))
  z[, 2:P] <- mvrnorm(T, meanz, corz)

  beta <- as.matrix(c(rep(alpha1, T/2), rep(alpha2, T/2)))
  # beta <- as.matrix(c(rep(alpha1, 3*T/10), rep(alpha2, 2*T/5), rep(alpha3, 3*T/10)))
  eta <- as.matrix( c(rep(eta_value, p1), rep(0, P - p1)) )

  if (dgp_order == 1) {
    x <- diag(arima.sim(n = T, list(ar = 0.5)), nrow = T, ncol = T)
    e <- as.matrix(rnorm(T, 0, vare))
  }else if  (dgp_order == 2) {
    x <- diag(rnorm(T, 0, 1), nrow = T, ncol = T )
    e <- vare * arima.sim(n = T, list(ar = 0.5), sd = sqrt(0.65))
  }else if  (dgp_order == 3) {
    x <- diag(arima.sim(n = T, list(ar = 0.45)), nrow = T, ncol = T)
    e <- as.matrix(c(rnorm(T/2, 0, vare), rnorm(T/2, 0, 0.8)))
  }else if  (dgp_order == 4) {
    x <- diag(arima.sim(n = T, list(ar = 0.5)), nrow = T, ncol = T)
    et <- rnorm(T, 0, 1)
    e <- vare * garch.sim(alpha = c(0.05, 0.05), beta = 0.9, n = T) * et
  }else if  (dgp_order == 5) {
    x <- diag(arima.sim(n = T, list(ar = 0.5)), nrow = T, ncol = T)
    e <- 0.5 * arima.sim(n = T, list(ma = 0.5), sd = sqrt(0.8))
  }

  y <- x %*% beta + z %*% eta + e

  dat_all = list(z = z, x = x, y = y)
  return(dat_all)
}
