# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# SIMULATION SCENARIO 3

library(mgcv)
library(xtable)
twopi <- 2*pi

# MCSE bands for all plots (given m=5000)
{
  m <- 5000
  p_mcse_plot <- seq(0,0.15,length=200)
  mcse <- sqrt(p_mcse_plot*(1-p_mcse_plot)/m)
  upper <- p_mcse_plot + 1.96*mcse
  lower <- p_mcse_plot - 1.96*mcse
}



# GAUSSIAN CASE ----------------------------------------------------------------
{
## 3.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(121)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim31.pdf", width = 8, height = 8)
  par(mfrow=c(2,2))
  
  # fixed parameters
  m <- 5000
  n <- 500
  sigma <- 0.5
  
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
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
    # construct z and y for current rho
    z_mat <- rho * g_mat + sqrt(1 - rho^2) * u_mat
    y_mat <- z_mat + eps_mat
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10), 
                           family = gaussian(link = "identity"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      rho          = rho,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("rho =", rho),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.1

## 3.2 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(122)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim32.pdf", width = 12, height = 4)
  par(mfrow=c(1,3))
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  sigma <- 0.5
  
  # varying parameter
  n_values <- c(50, 100, 250)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
    # DGP for current n
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
    eps_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
    g_mat <- scale(sin(twopi * x_mat))
    z_mat <- rho * g_mat + sqrt(1 - rho^2) * u_mat
    y_mat <- z_mat + eps_mat
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10), 
                           family = gaussian(link = "identity"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      n            = n,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("n =", n),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.2
} # END OF GAUSSIAN SIMULATION



# BINOMIAL CASE ----------------------------------------------------------------
{
## 3.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(221)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim31.pdf", width = 8, height = 8)
  par(mfrow=c(2,2))
  
  # fixed parameters
  m <- 5000
  n <- 500
  alpha <- 0.5
  lambda = "REML"
  
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
    
    # create empty vectors for the current sample size
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
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
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10), 
                           family = binomial(link = "logit"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      rho          = rho,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("rho =", rho),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.1

## 3.2 EFFECT OF SAMPLE SIZE
{
  # set seed for reproducibility
  set.seed(222)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim32.pdf", width = 12, height = 4)
  par(mfrow=c(1,3))
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  alpha <- 0.5
  
  # varying parameter
  n_values <- c(50, 100, 250)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
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
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10), 
                           family = binomial(link = "logit"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      n            = n,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("n =", n),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.2
} # END OF BINOMIAL SIMULATION



# POISSON CASE -----------------------------------------------------------------
{
## 3.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(321)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim31.pdf", width = 8, height = 8)
  par(mfrow=c(2,2))
  
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
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
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
      sim.gam <- mgcv::gam(y ~ z + s(x, k = 10),
                           family = poisson(link = "log"),
                           method = "REML",
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      rho          = rho,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("rho =", rho),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.1

## 3.2 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(322)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim32.pdf", width = 12, height = 4)
  par(mfrow=c(1,3))
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  alpha <- 0.5
  
  # varying parameter
  n_values <- c(50, 100, 250)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
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
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10), 
                           family = poisson(link = "log"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      n            = n,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("n =", n),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.2
} # END OF POISSON SIMULATION



# GAMMA CASE -------------------------------------------------------------------
{
## 3.1 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  # set seed for reproducibility
  set.seed(421)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim31.pdf", width = 8, height = 8)
  par(mfrow=c(2,2))
  
  # fixed parameters
  m <- 5000
  n <- 500
  alpha <- 0.5
  kappa <- 5
  
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
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
    # construct z and y for current rho
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- log(5) + alpha * z_mat
    mu_mat  <- exp(eta_mat)
    y_mat   <- matrix(rgamma(n * m, shape=kappa, scale=c(mu_mat) / kappa),
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
      sim.gam <- mgcv::gam(y ~ z + s(x, k = 10),
                           family = Gamma(link = "log"),
                           method = "REML",
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      rho          = rho,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("rho =", rho),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.1

## 3.2 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(422)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim32.pdf", width = 12, height = 4)
  par(mfrow=c(1,3))
  
  # fixed parameters
  m <- 5000
  rho <- 0.6
  alpha <- 0.5
  kappa <- 5
  
  # varying sample size
  n_values <- c(50, 100, 250)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors
    p_values <- numeric(m)
    edf <- numeric(m)
    beta <- numeric(m)
    concurv <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
    g_mat <- scale(sin(twopi * x_mat))
    z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
    eta_mat <- log(5) + alpha * z_mat
    mu_mat  <- exp(eta_mat)
    y_mat   <- matrix(rgamma(n * m, shape=kappa, scale=c(mu_mat)/kappa), 
                      nrow = n, ncol = m)
    
    # print data check (ensures 'well behaved' data)
    cat("\nn =", n,
        "\nmu in [", round(min(mu_mat), 3), ",", round(max(mu_mat), 3), "]",
        "\nmedian mean(lambda) across reps =", round(median(colMeans(mu_mat)), 3),
        "\nmedian mean(y) across reps =", round(median(colMeans(y_mat)), 3),
        "\nmax response observed =", max(y_mat), "\n")
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ z + s(x, k=10), 
                           family = Gamma(link = "log"), 
                           method = "REML",
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab["s(x)", "p-value"]
      edf[j] <- s_tab["s(x)", "edf"]
      beta[j] <- coef(sim.gam)["z"]
      concurv[j] <- mgcv::concurvity(sim.gam, full = TRUE)["estimate", "s(x)"]
    }
    
    # write row of data frame
    row_df <- data.frame(
      n            = n,
      F1           = mean(p_values <= 0.01),
      F5           = mean(p_values <= 0.05),
      F10          = mean(p_values <= 0.1),
      edf_mean     = mean(edf),
      edf_sd       = sd(edf),
      beta_mean    = mean(beta),
      beta_sd      = sd(beta),
      concurv_mean = mean(concurv),
      concurv_sd   = sd(concurv)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("n =", n),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
} # END OF SIMULATION 3.2
} # END OF GAMMA SIMULATION


