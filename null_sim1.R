# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# SIMULATION SCENARIO 1

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
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(101)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim11.pdf", width = 10, height = 6)
  par(mfrow=c(2,3))
  
  # constant parameter
  m <- 5000
  lambda <- "GCV.Cp"
  
  # varying parameter
  n_values <- c(25,50,100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors for the current sample size
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    y_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10), 
                           family = gaussian(link = "identity"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
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
  
} # END OF SIMULATION 1.1
  
## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # vary:
  # smoothing parameter estimation method (GCV, REML)
  # keep constant across simulations:
  # 5000 simulation replications
  # 1000 sample size n
  # basis dimension for smooth, k=10
  
  # set seed for reproducibility
  set.seed(102)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim12.pdf", width = 8, height = 4)
  par(mfrow=c(1,2))
  
  # constant parameters
  m <- 5000
  n <- 1000
  
  # varying parameter
  sp_method <- c("GCV.Cp", "REML")
  
  # prepare results vector
  results_list <- vector("list", length(sp_method))
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  y_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  
  # loop over values of n
  for (i in 1:length(sp_method))
  {
    # update smoothing parameter method for current loop
    lambda <- sp_method[i]
    
    # create empty vectors for the current method
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10), 
                           family = gaussian(link = "identity"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    # write row of data frame
    row_df <- data.frame(
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
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
  
} # END OF SIMULATION 1.2
  
## 1.3 EFFECT OF BASIS DIMENSION -----------------------------------------------
{
  # set seed for reproducibility
  set.seed(103)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim13.pdf", width = 8, height = 8)
  par(mfrow=c(2,2))
  
  # constant parameters
  m <- 5000
  n <- 250
  lambda <- "REML"
  
  # varying parameter
  k_values <- c(5,10,20,30)
  
  # prepare results vector
  results_list <- vector("list", length(k_values))
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  y_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  
  # loop over values of k
  for (i in 1:length(k_values))
  {
    # update k for current loop
    k <- k_values[i]
    
    # create empty vectors for the current basis dimension
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=k), 
                           family = gaussian(link = "identity"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    # write row of data frame for current basis dimension
    row_df <- data.frame(
      k          = k,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste("k =", k),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  }
  dev.off()
  
} # END OF SIMULATION 1.3
} # END OF GAUSSIAN CASE



# BINOMIAL CASE ----------------------------------------------------------------
{
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(201)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim11.pdf", width = 10, height = 6)
  par(mfrow=c(2,3))
  
  # fixed parameters
  m <- 5000
  lambda <- "GCV.Cp"
  
  # varying parameter
  n_values <- c(25,50,100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors for the current sample size
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    y_mat <- matrix(rbinom(n * m, size = 1, prob = 0.5), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10), 
                           family = binomial(link = "logit"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
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
  
} # END OF SIMULATION 1.1

## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(202)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim12.pdf", width = 8, height = 4)
  par(mfrow=c(1,2))
  
  # constant parameters
  m <- 5000
  n <- 1000
  
  # varying parameter
  sp_method <- c("GCV.Cp", "REML")
  
  # prepare results vector
  results_list <- vector("list", length(sp_method))
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  y_mat <- matrix(rbinom(n * m, size = 1, prob = 0.5), nrow = n, ncol = m)
  
  # loop over sp methods
  for (i in 1:length(sp_method))
  {
    # update smoothing parameter method for current loop
    lambda <- sp_method[i]
    
    # create empty vectors
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k=10), 
                           family = binomial(link = "logit"), 
                           method = lambda,
                           data=sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    # write row of data frame for current method
    row_df <- data.frame(
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
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

} # END OF SIMULATION 1.2
} # END OF BINOMIAL CASE



# POISSON CASE -----------------------------------------------------------------
{
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(301)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim11.pdf", width = 10, height = 6)
  par(mfrow=c(2,3))
  
  # constant parameters
  m <- 5000
  lambda <- "GCV.Cp"
  
  # varying parameter
  n_values <- c(25,50,100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors for the current sample size
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
    y_mat <- matrix(rpois(n * m, lambda = 5), nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      # reporting any warning messages if appear
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10),
                           family = poisson(link = "log"),
                           method = lambda,
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
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
  
} # END OF SIMULATION 1.1
  
## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(302)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim12.pdf", width = 8, height = 4)
  par(mfrow=c(1,2))
  
  # constant parameters
  m <- 5000
  n <- 1000
  
  # varying parameter
  sp_method <- c("GCV.Cp", "REML")
  
  # prepare results vector
  results_list <- vector("list", length(sp_method))
  
  # DGP
  x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
  y_mat <- matrix(rpois(n * m, lambda = 5), nrow = n, ncol = m)
  
  # loop over values of n
  for (i in 1:length(sp_method))
  {
    # update smoothing parameter method for current loop
    lambda <- sp_method[i]
    
    # create empty vectors for the current sample size
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10),
                           family = poisson(link = "log"),
                           method = lambda,
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    # write row of data frame for current sample size
    row_df <- data.frame(
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
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
  
} # END OF SIMULATION 1.2
} # END OF POISSON CASE



# GAMMA CASE -------------------------------------------------------------------
{
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(401)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim11.pdf", width = 10, height = 6)
  par(mfrow=c(2,3))
  
  # constant parameters
  m <- 5000
  lambda <- "GCV.Cp"
  mu <- 5
  kappa <- 5
  
  # varying parameter
  n_values <- c(25,50,100,250,500,1000)
  
  # prepare results vector
  results_list <- vector("list", length(n_values))
  
  # loop over values of n
  for (i in 1:length(n_values))
  {
    # update n for current loop
    n <- n_values[i]
    
    # create empty vectors for the current sample size
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # DGP
    x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
    y_mat <- matrix(rgamma(n * m, shape=kappa, scale=mu/kappa), 
                    nrow = n, ncol = m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10),
                           family = Gamma(link = "log"),
                           method = lambda,
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    
    # write row of data frame for current sample size
    row_df <- data.frame(
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
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
  
} # END OF SIMULATION 1.1

## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(402)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim12.pdf", width = 8, height = 4)
  par(mfrow=c(1,2))
  
  # constant parameters
  m <- 5000
  n <- 1000
  mu <- 5
  kappa <- 5
  
  # varying parameter
  sp_method <- c("GCV.Cp", "REML")
  
  # prepare results vector
  results_list <- vector("list", length(sp_method))
  
  # DGP
  x_mat <- matrix(runif(n * m, 0, 1), nrow = n, ncol = m)
  y_mat <- matrix(rgamma(n * m, shape=kappa, scale=mu/kappa), nrow = n, ncol = m)
  
  # loop over values of n
  for (i in 1:length(sp_method))
  {
    # update smoothing parameter method for current loop
    lambda <- sp_method[i]
    
    # create empty vectors
    p_values <- numeric(m)
    edf <- numeric(m)
    
    # loop over number of simulation replications
    for (j in 1:m)
    {
      # form data frame and fit GAM
      sim <- data.frame(x = x_mat[, j], y = y_mat[, j])
      sim.gam <- mgcv::gam(y ~ s(x, k = 10),
                           family = Gamma(link = "log"),
                           method = lambda,
                           data = sim)
      
      # extract current outputs for results
      s_tab <- summary(sim.gam)$s.table
      p_values[j] <- s_tab[1, "p-value"]
      edf[j] <- s_tab[1, "edf"]
    }
    # write row of data frame for current sample size
    row_df <- data.frame(
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
    # export plot
    plot(ecdf(p_values),
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

} # END OF SIMULATION 1.2
} # END OF GAMMA CASE
