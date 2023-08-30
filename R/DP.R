
DP <- function(T, P, p1, q, meanz, corZ, vare, eta_value, alpha1, alpha2, dgp_order, iteration){

  for (ii in 1:iteration) {

    data <- dgp(T, P, p1, q, meanz, corZ, vare, eta_value, alpha1, alpha2, dgp_order)
    y <- data$y
    z <- data$z
    x <- data$x

    L <- diag(1, T, T)
    for (i in 1:(T - 1)) {
      L[i + 1, i] <- -1
    }

    AM <- kronecker(L, diag(1, q, q))

    Q_z <- diag(1, T, T) - z %*% solve(t(z) %*% z + 10 * diag(1, P, P)) %*% t(z)
    beta_initial <- solve(t(x) %*% Q_z %*% x + 10 * t(AM) %*% AM) %*% t(x) %*% Q_z %*% y
    eta_initial <- solve(t(z) %*% z + 10 * diag(1, P, P)) %*% t(z) %*% (y - x %*% beta_initial)

    beta_weight <- 1/abs(c(beta_initial[1], diff(beta_initial)))
    eta_weight <- 1/abs(eta_initial)

    beta_count <- array(0, dim = c(T, Lambda_len))
    eta_count <- array(0, dim = c(P, Lambda_len, Gamma_len))
    eta_count_ols <- array(0, dim = c(P, Lambda_len, Gamma_len))
    break_count <- array(0, dim = c(1, Lambda_len))
    breakpos_count <- array(0, dim = c(T, Lambda_len))
    IC <- matrix(0, Gamma_len, Lambda_len)

    M <- list(sense = "min")
    b1 <- as.matrix.csr(L)
    B1 <- cbind(b1, -as(T, "matrix.diag.csr"), as(T, "matrix.diag.csr"),
                as.matrix.csr(0, T, 2*T), as.matrix.csr(0, T, 3*P), as.matrix.csr(0, T, T + 3) )
    B2 <- cbind(as(T, "matrix.diag.csr"),  as.matrix.csr(0, T, 2*T),
                -as(T, "matrix.diag.csr"), as(T, "matrix.diag.csr"),
                as.matrix.csr(0, T, 3*P), as.matrix.csr(0, T, T + 3))
    B3 <- cbind(as.matrix.csr(0, P, 5*T),  as(P, "matrix.diag.csr"),
                -as(P, "matrix.diag.csr"),
                as(P, "matrix.diag.csr"), as.matrix.csr(0, P, T + 3))

    a1 <- as.matrix.csr(x)
    a2 <- as.matrix.csr(z)
    A1 <- cbind(a1, as.matrix.csr(0, T, T*4), a2, as.matrix.csr(0, T, P*2),
                as(T, "matrix.diag.csr"), as.matrix.csr(0, T, 3))
    A2 <- cbind(as.matrix.csr(0, 1, T*6 + 3*P), as.matrix.csr(c(-.5, 1, 0), 1, 3) )
    A3 <- cbind(as.matrix.csr(0, 1, T*6 + 3*P), as.matrix.csr(c(-.5, 0, 1), 1, 3) )
    A <- rbind(B1, B2, B3, A1, A2, A3)
    M$A <- as(A,"CsparseMatrix")

    M$bc <- rbind(c(rep(0, T),rep(0, T),rep(0, P), y, -0.5, 0.5),
                  c(rep(0, T),rep(0, T),rep(0, P), y, -0.5, 0.5)   )

    M$bx <- rbind(c(rep(-Inf, T), rep(0, 2*T),rep(0, 2*T),
                    rep(-Inf, P), rep(0, 2*P),
                    rep(-Inf, T),
                    rep(0, 3) ),
                  c(rep(Inf, 6*T + 3*P + 3))        )
    M$cones <- matrix(list("QUAD", c(6*T + 3 + 3*P, (5*T + 3*P + 1):(5*T + 3*P + T), 6*T + 2 + 3*P)), 2, 1)
    rownames(M$cones) <- c("type", "sub")

    for (l in 1:Lambda_len) {
      for (g in 1:Gamma_len) {

        # M$c <- c(rep(0, T), rep(Lambda[g], 2*T), rep(0, 2*T), rep(0, P),
        #          rep(Gamma[g], 2*P), rep(0, T), 1/T, 0, 0)

        M$c <- c(rep(0, T), Lambda[l]*beta_weight, Lambda[l]*beta_weight, rep(0, 2*T), rep(0, P),
                 Gamma[g] * eta_weight, Gamma[g] * eta_weight, rep(0, T), 1/T, 0, 0)

        verb <- 2
        result <- mosek(M, opts = list(verbose = verb))
        xx <- result$sol$itr$xx

        beta_count[, l]  <- xx[1:T]
        beta_iteration_test[, l, ii] <- xx[1:T]
        eta_count[, l, g] <- xx[(5*T + 1):(5*T + P)]

        breakpoint <- 0
        breakpos <- 0
        for (ll in 1:(T - 1)) {
          if (abs(beta_count[ll + 1, l] - beta_count[ll, l]) >= 0.1) {
            breakpoint <- breakpoint + 1
            breakpos <- c(breakpos, ll)
          }
        }
        break_count[, l] <- breakpoint

        breakpos_count[, l] <- array(0, dim = c(T, 1))
        for (pp in 1:length(breakpos)) {
          breakpos_count[pp, l] <- breakpos[pp]
        }

        nonzero_pos <- which(abs(eta_count[, l, g]) >= (eta_value/10))
        eta_count_ols[, l, g] <- array(0, dim = c(P, 1))
        for (lp in 1:length(nonzero_pos)) {
          eta_count_ols[, l, g][nonzero_pos[lp]] <- coef(lm(y ~ diag(x) + z[, nonzero_pos] - 1))[lp + 1]
        }

        pos <- breakpos_count[, l][which(breakpos_count[, l] != 0)]
        pos <- c(0, pos, T)

        for (ls in 2:length(pos)) {
          y_ols <- y[(pos[ls - 1] + 1):pos[ls]]
          x_ols <- diag(x)[(pos[ls - 1] + 1):pos[ls]]
          z_ols <- z[(pos[ls - 1] + 1):pos[ls], nonzero_pos]
          if (length(y_ols) <= 1) {
            beta_count[(pos[ls - 1] + 1):pos[ls], l] <- as.matrix(beta_count[(pos[ls - 1] + 1):pos[ls], l],
                                                                  nrow = (pos[ls] - pos[ls - 1]), ncol = 1)
          }else{
            beta_ols <- coef(lm(y_ols ~ x_ols + z_ols - 1))[1]
            beta_count[(pos[ls - 1] + 1):pos[ls], l] <- as.matrix(beta_ols,
                                                                  nrow = (pos[ls] - pos[ls - 1]), ncol = 1)
          }
        }
      }

      for (g in 1:Gamma_len) {
        sum_all <- mean((y - z %*% eta_count_ols[, l, g] - x %*% beta_count[, l])^2)
        IC[g, l] <- log(sum_all) + ((break_count[, l] + 1) * q + P) * log(T * q + P) * sqrt(1/T)
      }
    }

    opt <- which(IC == min(IC), arr.ind = TRUE)
    if ((dim(opt)[1]) >= 2) {
      opt[1] <- opt[dim(opt)[1], 1]
      opt[2] <- opt[dim(opt)[1], 2]
    }

    beta_iteration_count[, ii] <- beta_count[, opt[2]]
    eta_iteration_count[, ii] <- eta_count[, opt[2], opt[1]]
    eta_iteration_count_ols[, ii] <- eta_count_ols[, opt[2], opt[1]]
    break_iteration_count[, ii] <- break_count[, opt[2]]
    breakpos_iteration_count[, ii] <- breakpos_count[, opt[2]]

    if (break_iteration_count[, ii] == 1) {
      hd_iteration_count[, ii] <- hausdorff_dist(breakpos_iteration_count[, ii][2], T/2)
    }
  }

  output <- list(
    beta_iteration_count <- beta_iteration_count,
    eta_iteration_count <- eta_iteration_count,
    eta_iteration_count_ols <- eta_iteration_count_ols,
    break_iteration_count <- break_iteration_count,
    breakpos_iteration_count <- breakpos_iteration_count,
    hd_iteration_count <- hd_iteration_count
  )

  return(output)
}
