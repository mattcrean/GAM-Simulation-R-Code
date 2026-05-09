# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# SIMULATION SCENARIO 4
{
  
library(mgcv)
twopi <- 2*pi

# MCSE bands for all plots (given m=5000)
{
  m <- 5000
  p_mcse_plot <- seq(0,0.15,length=200)
  mcse <- sqrt(p_mcse_plot*(1-p_mcse_plot)/m)
  upper <- p_mcse_plot + 1.96*mcse
  lower <- p_mcse_plot - 1.96*mcse
}

results_list_n <- vector("list", 4)
results_list_rho <- vector("list", 4)

# GAUSSIAN CASE ----------------------------------------------------------------
{
## 4.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(131)
  
  # fixed parameters
  m <- 5000
  sigma <- 0.5
  n <- 500
  
  # varying parameter
  rho_values <- c(0, 0.3, 0.6, 0.9)
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  eps_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  
  # prepare results vector
  results_list <- vector("list", length(rho_values))
  
  # loop over values of rho
  for (i in 1:length(rho_values))
  {
    # update rho for current loop
    rho <- rho_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # construct z and y for current rho
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    y_mat <- z_mat + eps_mat
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k=10, m=c(2,0)), 
                           family = gaussian(link = "identity"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      Response      = "Gaussian",
      rho           = rho,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.1

results_list_rho[[1]] <- results_list
  
## 4.2 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(132)
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  sigma <- 0.5
  
  # varying parameter
  n_values <- c(100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # DGP for each n
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
    eps_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
    g_mat <- scale(sin(twopi * x_mat))
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    y_mat <- z_mat + eps_mat
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k=10, m=c(2,0)), 
                           family = gaussian(link = "identity"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      Response      = "Gaussian",
      n             = n,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.2
  
results_list_n[[1]] <- results_list
  
} # END OF GAUSSIAN SIMULATION



# BINOMIAL CASE ----------------------------------------------------------------
{
## 4.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(231)
  
  # fixed parameters
  m <- 5000
  n <- 500
  alpha <- 0.5
  
  # varying parameter
  rho_values <- c(0, 0.3, 0.6, 0.9)
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  
  # prepare results vector
  results_list <- vector("list", length(rho_values))
  
  # loop over values of rho
  for (i in 1:length(rho_values))
  {
    # update rho for current loop
    rho <- rho_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # construct z and y for current rho
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- alpha * z_mat
    success_mat <- plogis(eta_mat)
    y_mat <- matrix(rbinom(n * m, 1, success_mat), nrow = n)
    
    # print data check (ensures 'well behaved' data)
    cat("\nrho =", rho,
        "\n p in [", round(min(success_mat), 3), ",", round(max(success_mat), 3), 
        "] \n",
        " event rate median =", round(median(colMeans(y_mat)), 3), "\n")
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10)+ s(z, k = 10, m = c(2, 0)), 
                           family = binomial(link = "logit"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      Response      = "Binomial",
      rho           = rho,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.1

results_list_rho[[2]] <- results_list
  
## 4.2 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(232)
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  alpha <- 0.5
  
  # varying sample size
  n_values <- c(100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
    g_mat <- scale(sin(twopi * x_mat))
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- alpha * z_mat
    success_mat <- plogis(eta_mat)
    y_mat <- matrix(rbinom(n * m, 1, success_mat), nrow = n)
    
    # print data check (ensures 'well behaved' data)
    cat("\nn =", n,
        "\n p in [", round(min(success_mat), 3), ",", round(max(success_mat), 3), 
        "] \n",
        " event rate median =", round(median(colMeans(y_mat)), 3), "\n")
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2, 0)), 
                           family = binomial(link = "logit"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      Response      = "Binomial",
      n             = n,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.2

results_list_n[[2]] <- results_list
  
## 4.3 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(233)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim4_binomial_method.pdf", width = 8, height = 4)
  par(mfrow=c(1,2))
  
  # fixed parameters
  m <- 5000
  n <- 500
  alpha <- 0.5
  rho <- 0.6
  
  # varying parameter
  sp_method <- c("GCV.Cp", "REML")
  
  # prepare results vector
  results_list <- vector("list", length(sp_method))
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
  eta_mat <- alpha * z_mat
  success_mat <- plogis(eta_mat)
  y_mat <- matrix(rbinom(n * m, 1, success_mat), nrow = n)
  
  # print data check (ensures 'well behaved' data)
  cat("\n p in [", round(min(success_mat), 3), ",", round(max(success_mat), 3), 
      "] \n",
      " event rate median =", round(median(colMeans(y_mat)), 3), "\n")
  
  # loop over methods
  for (i in 1:length(sp_method))
  {
    # update smoothing parameter method for current loop
    lambda <- sp_method[i]
    
    # create empty vectors for the current method
    p_values_z <- numeric(m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2, 0)), 
                           family = binomial(link = "logit"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_z[j] <- s_tab["s(z)", "p-value"]
    }
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = lambda,
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  
  dev.off()
  
} # END OF SIMULATION 4.3
  
  
  
} # END OF BINOMIAL SIMULATION



# POISSON CASE -----------------------------------------------------------------
{
## 4.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(331)
  
  # fixed parameters
  m <- 5000
  n <- 500
  alpha <- 0.5
  
  # varying parameter
  rho_values <- c(0, 0.3, 0.6, 0.9)
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  
  # prepare results vector
  results_list <- vector("list", length(rho_values))
  
  # loop over values of rho
  for (i in 1:length(rho_values))
  {
    # update rho for current loop
    rho <- rho_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # construct z and y for current rho
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- log(5) + alpha * z_mat
    lambda_mat  <- exp(eta_mat)
    y_mat   <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
    
    # print data check
    cat("\nrho =", rho,
        "\nlambda in [", round(min(lambda_mat), 3), ",", round(max(lambda_mat), 3), "]",
        "\nmedian mean(lambda) across reps =", round(median(colMeans(lambda_mat)), 3),
        "\nmedian mean(y) across reps =", round(median(colMeans(y_mat)), 3),
        "\nmax response observed =", max(y_mat), "\n")
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10)+ s(z, k = 10, m = c(2, 0)), 
                           family = poisson(link = "log"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      Response      = "Poisson",
      rho           = rho,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.1

results_list_rho[[3]] <- results_list
  
## 4.2 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(332)
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  alpha <- 0.5
  
  # varying parameter
  n_values <- c(100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors for the current sample size
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
    g_mat <- scale(sin(twopi * x_mat))
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- log(5) + alpha * z_mat
    lambda_mat  <- exp(eta_mat)
    y_mat   <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
    
    # print data check
    cat("\nrho =", rho,
        "\nlambda in [", round(min(lambda_mat), 3), ",", round(max(lambda_mat), 3), "]",
        "\nmedian mean(lambda) across reps =", round(median(colMeans(lambda_mat)), 3),
        "\nmedian mean(y) across reps =", round(median(colMeans(y_mat)), 3),
        "\nmax response observed =", max(y_mat), "\n")
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2, 0)), 
                           family = poisson(link = "log"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      Response      = "Poisson",
      n             = n,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.2
  
results_list_n[[3]] <- results_list
  
} # END OF POISSON SIMULATION




# GAMMA CASE -------------------------------------------------------------------
{
## 4.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(431)
  
  # fixed parameters
  m <- 5000
  n <- 500
  alpha <- 0.5
  kappa <- 5
  
  # varying correlation parameter
  rho_values <- c(0, 0.3, 0.6, 0.9)
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  
  # prepare results vector
  results_list <- vector("list", length(rho_values))
  
  # loop over values of rho
  for (i in 1:length(rho_values))
  {
    # update rho for current loop
    rho <- rho_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # construct z and y for current rho
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- log(5) + alpha * z_mat
    mu_mat  <- exp(eta_mat)
    y_mat   <- matrix(rgamma(n * m, shape = kappa, scale = c(mu_mat) / kappa), 
                      nrow = n, ncol = m)
    
    # print data check
    cat("\nrho =", rho,
        "\nmu in [", round(min(mu_mat), 3), ",", round(max(mu_mat), 3), "]",
        "\nmedian mean(mu) across reps =", round(median(colMeans(mu_mat)), 3),
        "\nmedian mean(y) across reps =", round(median(colMeans(y_mat)), 3),
        "\nmax response observed =", max(y_mat), "\n")
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10)+ s(z, k = 10, m = c(2, 0)), 
                           family = Gamma(link = "log"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame for current sample size
    row_df <- data.frame(
      Response      = "Gamma",
      rho           = rho,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.1

results_list_rho[[4]] <- results_list
  
## 4.2 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(432)
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  alpha <- 0.5
  kappa <- 5
  
  # varying parameter
  n_values <- c(100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of rho
  for (i in 1:length(n_values))
  {
    # update rho for current loop
    n <- n_values[i]
    
    # create empty vectors for the current sample size
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    beta <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
    g_mat <- scale(sin(twopi * x_mat))
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- log(5) + alpha * z_mat
    mu_mat  <- exp(eta_mat)
    y_mat   <- matrix(rgamma(n * m, shape = kappa, scale = c(mu_mat) / kappa), 
                      nrow = n, ncol = m)
    
    # print data check
    cat("\nrho =", rho,
        "\nmu in [", round(min(mu_mat), 3), ",", round(max(mu_mat), 3), "]",
        "\nmedian mean(mu) across reps =", round(median(colMeans(mu_mat)), 3),
        "\nmedian mean(y) across reps =", round(median(colMeans(y_mat)), 3),
        "\nmax response observed =", max(y_mat), "\n")
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2, 0)), 
                           family = Gamma(link = "log"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
    }
    
    # write row of data frame for current sample size
    row_df <- data.frame(
      Response      = "Gamma",
      n             = n,
      Covariate     = c("x", "z"),
      F1            = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5            = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10           = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean      = c(mean(edf_x), mean(edf_z)),
      edf_sd        = c(sd(edf_x), sd(edf_z)),
      beta_mean     = c(NA, mean(beta)),
      beta_sd       = c(NA, sd(beta))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
  }
} # END OF SIMULATION 4.2
  
results_list_n[[4]] <- results_list
  
} # END OF GAMMA SIMULATION



# GLM COMPARISON ---------------------------------------------------------------
{
  
  # set seed for reproducibility
  set.seed(400)
  
  # define 2*pi
  twopi <- 2*pi
  
  # fixed parameters
  m <- 5000
  n <- 500
  rho <- 0.6
  sigma <- 0.5
  kappa <- 5
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
  
  # response families
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  
  results_list <- vector("list", length(response))
  
  for (i in 1:length(response))
  {
    # update family for current loop
    response_family <- response[i]
    
    # set true beta under current family
    if (response_family == "Gaussian") {
      beta_true <- 1
    } else {
      beta_true <- 0.5
    }
    
    # create empty vectors for outputs
    beta_gam <- numeric(m)
    se_gam <- numeric(m)
    
    beta_glm <- numeric(m)
    se_glm <- numeric(m)
    
    p_values_x <- numeric(m)
    p_values_z <- numeric(m)
    
    glm_p_values <- numeric(m)
    gam_p_values <- numeric(m)
    
    
    # generate y under current family
    if (response_family == "Gaussian")
    {
      eps_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
      y_mat <- beta_true * z_mat + eps_mat
    }
    
    if (response_family == "Binomial")
    {
      eta_mat <- beta_true * z_mat
      q_mat <- plogis(eta_mat)
      y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
    }
    
    if (response_family == "Poisson")
    {
      eta_mat <- log(5) + beta_true * z_mat
      lambda_mat <- exp(eta_mat)
      y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
    }
    
    if (response_family == "Gamma")
    {
      eta_mat <- log(5) + beta_true * z_mat
      mu_mat <- exp(eta_mat)
      y_mat <- matrix(rgamma(n * m, shape = kappa, scale = c(mu_mat) / kappa),
                      nrow = n, ncol = m)
    }
    
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      
      # fit GAM and GLM under current response family
      if (response_family == "Gaussian")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k=10, m=c(2,0)),
                             family = gaussian(link = "identity"),
                             method = "REML",
                             data = sim)
        
        sim.glm <- glm(y ~ z,
                       family = gaussian(link = "identity"),
                       data = sim)
      }
      
      if (response_family == "Binomial")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k=10, m=c(2,0)),
                             family = binomial(link = "logit"),
                             method = "REML",
                             data = sim)
        
        sim.glm <- glm(y ~ z,
                       family = binomial(link = "logit"),
                       data = sim)
      }
      
      if (response_family == "Poisson")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k=10, m=c(2,0)),
                             family = poisson(link = "log"),
                             method = "REML",
                             data = sim)
        
        sim.glm <- glm(y ~ z,
                       family = poisson(link = "log"),
                       data = sim)
      }
      
      if (response_family == "Gamma")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k=10, m=c(2,0)),
                             family = Gamma(link = "log"),
                             method = "REML",
                             data = sim)
        
        sim.glm <- glm(y ~ z,
                       family = Gamma(link = "log"),
                       data = sim)
      }
      
      # extract GAM smooth p-values
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      
      # extract GAM beta outputs
      p_tab_gam <- summary(sim.gam)$p.table
      beta_gam[j] <- p_tab_gam["z", "Estimate"]
      gam_p_values[j] <- p_tab_gam["z", ncol(p_tab_gam)]
      
      # extract GLM beta outputs
      p_tab_glm <- summary(sim.glm)$coefficients
      beta_glm[j] <- p_tab_glm["z", "Estimate"]
      glm_p_values[j] <- p_tab_glm["z", ncol(p_tab_glm)]
    }
    
    # write rows of data frame for current family
    row_df <- data.frame(
      Response      = response_family,
      Model         = c("GAM", "GLM"),
      beta_mean     = c(mean(beta_gam), mean(beta_glm)),
      beta_sd       = c(sd(beta_gam), sd(beta_glm)),
      F1x           = c(mean(p_values_x <= 0.01), NA),
      F5x           = c(mean(p_values_x <= 0.05), NA),
      F10x          = c(mean(p_values_x <= 0.1), NA),
      F1z           = c(mean(p_values_z <= 0.01), NA),
      F5z           = c(mean(p_values_z <= 0.05), NA),
      F10z          = c(mean(p_values_z <= 0.1), NA),
      F1betaz       = c(mean(gam_p_values <= 0.01), mean(glm_p_values <= 0.01)),
      F5betaz       = c(mean(gam_p_values <= 0.05), mean(glm_p_values <= 0.05)),
      F10betaz      = c(mean(gam_p_values <= 0.1), mean(glm_p_values <= 0.1))
    )
    
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }

} # END OF GLM COMPARISON


writeLines("\nEmpirical rejection rates and EDF summaries for the smooth term of each covariate.\n")
writeLines("For rows with Covariate = z, beta_mean and beta_sd refer to the separate linear z coefficient.\n\n")

results_rho <- do.call(
  rbind,
  lapply(results_list_rho, function(x) do.call(rbind, x))
)
rownames(results_rho) <- NULL
print(results_rho)

writeLines("\n")

results_n <- do.call(
  rbind,
  lapply(results_list_n, function(x) do.call(rbind, x))
)
rownames(results_n) <- NULL
print(results_n)

writeLines("\n")

results_glm <- do.call(rbind, results_list)
rownames(results_glm) <- NULL
print(results_glm)

}
