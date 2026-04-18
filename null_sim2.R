# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# SIMULATION SCENARIO 2

library(mgcv)
library(xtable)

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
## 2.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(111)
  
  # fixed parameters
  m <- 5000
  sigma = 0.5
  lambda = "REML"
  
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
    
    # DGP
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    z_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    y_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10), 
                           family = gaussian(link = "identity"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim21_n", n_values[i], ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
    
  }
} # END OF SIMULATION 2.1

## 2.2 EFFECT OF ERROR TERM VARIANCE -------------------------------------------
{
  # set seed for reproducibility
  set.seed(112)
  
  # fixed parameters
  m <- 5000
  n <- 500
  lambda <- "REML"
  
  # varying parameter
  sigma_values <- c(0.25, 0.5, 1, 2, 5)
  
  # prepare results vector
  results_list <- vector("list", length(sigma_values))
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  z_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  
  # loop over values of n
  for (i in 1:length(sigma_values))
  {
    # update sigma
    sigma <- sigma_values[i]
    
    # generate error/y term with corresponding variance for current loop
    y_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10), 
                           family = gaussian(link = "identity"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame for current value of sigma
    row_df <- data.frame(
      sigma      = sigma,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim22_sigma", sigma_values[i], ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
  }
} # END OF SIMULATION 2.2
} # END OF GAUSSIAN SIMULATION



# BINOMIAL CASE ----------------------------------------------------------------
{
## 2.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(211)
  
  # fixed parameters
  m <- 5000
  lambda = "REML"
  q <- 0.5
  
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
    
    # DGP
    # pre-generate all data for this sample size
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    z_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    y_mat <- matrix(rbinom(n * m, size = 1, prob = q), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10), 
                           family = binomial(link = "logit"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim21_n", n_values[i], ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
    
  }
} # END OF SIMULATION 2.1

## 2.2 EFFECT OF SUCCESS PROBABILITY FOR Y -------------------------------------
{
  # set seed for reproducibility
  set.seed(212)
  
  # fixed parameters
  m <- 5000
  n <- 500
  lambda = "REML"
  
  # varying parameter
  success_values <- c(0.3, 0.5, 0.7)
  
  # prepare results vector
  results_list <- vector("list", length(success_values))
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  z_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  
  # loop over values of success probability
  for (i in 1:length(success_values))
  {
    # update success probability for current loop
    success <- success_values[i]
    
    # create empty vectors for the current q
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    
    # DGP for current probability of success
    y_mat <- matrix(rbinom(n * m, size = 1, prob = success), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10), 
                           family = binomial(link = "logit"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame for q
    row_df <- data.frame(
      success    = success,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim22_success", success, ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
  }
} # END OF SIMULATION 2.2
} # END OF BINOMIAL SIMULATION



# POISSON CASE -----------------------------------------------------------------
{
## 2.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(311)
  
  # fixed parameters
  m <- 5000
  lambda = "REML"
  poi_mean <- 5
  
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
    
    # DGP
    x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
    z_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
    y_mat <- matrix(rpois(n * m, lambda = poi_mean), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10) + s(z, k = 10),
                           family = poisson(link = "log"),
                           method = "REML",
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim21_n", n_values[i], ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
    
  }
} # END OF SIMULATION 2.1

## 2.2 EFFECT OF MEAN OF POISSON DISTRIBUTION ----------------------------------
{
  # set seed for reproducibility
  set.seed(312)
  
  # fixed parameters
  m <- 5000
  n <- 500
  lambda = "REML"
  
  # varying parameter
  lambda_values <- c(3, 5, 7, 10)
  
  # prepare results vector
  results_list <- vector("list", length(lambda_values))
  
  # DGP
  x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
  z_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
  
  # loop over values of poisson mean
  for (i in 1:length(lambda_values))
  {
    # update poisson mean for current loop
    lam <- lambda_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    
    # DGP
    y_mat <- matrix(rpois(n * m, lambda = lam), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10) + s(z, k = 10),
                           family = poisson(link = "log"),
                           method = "REML",
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame
    row_df <- data.frame(
      lam        = lam,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim22_lambda", lam, ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
  }
} # END OF SIMULATION 2.2
} # END OF POISSON SIMULATION



# GAMMA CASE -------------------------------------------------------------------
{
## 2.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(411)
  
  # fixed parameters
  m <- 5000
  mu <- 5
  kappa <- 5
  lambda = "REML"
  
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
    
    # DGP
    x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
    z_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
    y_mat <- matrix(rgamma(n * m, shape = kappa, scale = mu / kappa),
                    nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10) + s(z, k = 10),
                           family = Gamma(link = "log"),
                           method = lambda,
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim21_n", n_values[i], ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
  }
} # END OF SIMULATION 2.1

## 2.2 EFFECT OF MEAN (AND VARIANCE) OF GAMMA DISTRIBUTION ---------------------
{
  # set seed for reproducibility
  set.seed(412)
  
  # fixed parameters
  m <- 5000
  n <- 500
  kappa <- 5
  lambda = "REML"
  
  # varying parameter
  mu_values <- c(3, 5, 7, 10)
  
  # prepare results vector
  results_list <- vector("list", length(mu_values))
  
  # DGP
  x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
  z_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
  
  # loop over values of gamma mean
  for (i in 1:length(mu_values))
  {
    # update gamma mean for current loop
    mu <- mu_values[i]
    
    # create empty vectors
    p_values_x <- numeric(m)
    edf_x <- numeric(m)
    p_values_z <- numeric(m)
    edf_z <- numeric(m)
    
    # DGP for current probability of success
    y_mat <- matrix(rgamma(n * m, shape = kappa, scale = mu / kappa),
                    nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10) + s(z, k = 10),
                           family = Gamma(link = "log"),
                           method = "REML",
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values_x[j] <- s_tab["s(x)", "p-value"]
      edf_x[j] <- s_tab["s(x)", "edf"]
      p_values_z[j] <- s_tab["s(z)", "p-value"]
      edf_z[j] <- s_tab["s(z)", "edf"]
    }
    # write row of data frame
    row_df <- data.frame(
      mu         = mu,
      F1x        = mean(p_values_x <= 0.01),
      F5x        = mean(p_values_x <= 0.05),
      F10x       = mean(p_values_x <= 0.1),
      edf_meanx  = mean(edf_x),
      edf_sdx    = sd(edf_x),
      F1z        = mean(p_values_z <= 0.01),
      F5z        = mean(p_values_z <= 0.05),
      F10z       = mean(p_values_z <= 0.1),
      edf_meanz  = mean(edf_z),
      edf_sdz    = sd(edf_z)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # prepare plot grid and pdf file to be saved to
    title <- paste0("sim22_mu", mu, ".pdf")
    pdf(title, width = 8, height = 4)
    par(mfrow=c(1,2))
    
    # export plot
    plot(ecdf(p_values_x),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(x)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty=2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    # export plot
    plot(ecdf(p_values_z),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0("s(z)"),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
    dev.off()
  }
} # END OF SIMULATION 2.2
} # END OF GAMMA SIMULATION


