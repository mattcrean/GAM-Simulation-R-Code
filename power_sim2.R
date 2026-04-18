# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# POWER SIMULATION SCENARIO 2

library(mgcv)
library(xtable)
twopi <- 2*pi

## 2.1 EFFECT OF SIGNAL STRENGTH -----------------------------------------------
{
  # parameters for current simulation
  set.seed(2100)
  m <- 5000
  n <- 500
  lambda <- "REML"
  twopi <- 2*pi
  rho <- 0.6
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  response_colour <- c("black", "red", "blue", "darkgreen")
  
  results_power21_list <- vector("list", 16)
  row_index <- 0
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power_21.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
  h_mat <- scale(sin(twopi * z_mat))
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    
    # generate y under current family
    if (response_family == "Gaussian")
    {
      delta_values <- c(0.01, 0.05, 0.10, 0.25)
      eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
    }
    
    if (response_family == "Binomial")
    {
      delta_values <- c(0.10, 0.25, 0.50, 1.00)
    }
    
    if (response_family == "Poisson")
    {
      delta_values <- c(0.01, 0.05, 0.10, 0.25)
    }
    
    if (response_family == "Gamma")
    {
      delta_values <- c(0.01, 0.05, 0.10, 0.25)
    }
    
    p_values_x_list <- vector("list", length(delta_values))
    p_values_z_list <- vector("list", length(delta_values))
    
    for (d in 1:length(delta_values))
    {
      row_index <- row_index + 1
      delta <- delta_values[d]
      p_values_x <- numeric(m)
      p_values_z <- numeric(m)
      
      if (response_family == "Gaussian")
      {
        y_mat <- z_mat + delta * h_mat + eps_mat
      }
      
      if (response_family == "Binomial")
      {
        eta_mat <- 0.5 * z_mat + delta * h_mat
        q_mat <- plogis(eta_mat)
        y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Poisson")
      {
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        lambda_mat <- exp(eta_mat)
        y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Gamma")
      {
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        mu_mat <- exp(eta_mat)
        y_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_mat) / 5),
                        nrow = n, ncol = m)
      }
      
      # loop over number of simulation replications
      for (k in 1:m)
      {
        # form data frame
        sim <- data.frame(x = x_mat[, k], y = y_mat[, k], z = z_mat[, k])
        
        # fit GAM and GLM under current response family
        if (response_family == "Gaussian")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = gaussian(link = "identity"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Binomial")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = binomial(link = "logit"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Poisson")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = poisson(link = "log"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Gamma")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = Gamma(link = "log"),
                               method = lambda,
                               data = sim)
        }
        
        s_tab <- summary(sim.gam)$s.table
        p_values_x[k] <- s_tab["s(x)", "p-value"]
        p_values_z[k] <- s_tab["s(z)", "p-value"]
        
        
      } # END OF m SIMULATIONS LOOP
      
      p_values_x_list[[d]] <- p_values_x
      p_values_z_list[[d]] <- p_values_z
      
      row_df <- data.frame(
        family        = response_family,
        delta         = delta,
        F1x           = mean(p_values_x <= 0.01),
        F5x           = mean(p_values_x <= 0.05),
        F10x          = mean(p_values_x <= 0.1)
      )
      
      results_power21_list[[row_index]] <- row_df
      
    } # END OF DELTA LOOP
    
    # export plot
    plot(ecdf(p_values_z_list[[1]]),
         col = "black",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = paste0(response_family),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    lines(ecdf(p_values_z_list[[2]]), col = "red")
    lines(ecdf(p_values_z_list[[3]]), col = "blue")
    lines(ecdf(p_values_z_list[[4]]), col = "darkgreen")
    
  } # END OF RESPONSE LOOP
  
  dev.off()
  
} # END OF SIMULATION 2.1



## 2.2 EFFECT OF CORRELATION PARAMETER -----------------------------------------
{
  
  # parameters for current simulation
  set.seed(2200)
  m <- 5000
  n <- 500
  lambda <- "REML"
  twopi <- 2*pi
  rho_values <- c(0,0.3,0.6,0.9)
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  response_colour <- c("black", "red", "blue", "darkgreen")
  
  results_power22_list <- vector("list", 16)
  row_index <- 0
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power_22.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    p_values_x_list <- vector("list", length(rho_values))
    p_values_z_list <- vector("list", length(rho_values))
    
    for (i in 1:length(rho_values))
    {
      rho <- rho_values[i]
      z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
      h_mat <- scale(sin(twopi * z_mat))
      row_index <- row_index + 1
      
      if (response_family == "Gaussian")
      {
        delta <- 0.05
        eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
        y_mat <- z_mat + delta * h_mat + eps_mat
      }
      
      if (response_family == "Binomial")
      {
        delta <- 0.25
        eta_mat <- 0.5 * z_mat + delta * h_mat
        q_mat <- plogis(eta_mat)
        y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Poisson")
      {
        delta <- 0.05
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        lambda_mat <- exp(eta_mat)
        y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Gamma")
      {
        delta <- 0.05
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        mu_mat <- exp(eta_mat)
        y_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_mat) / 5),
                        nrow = n, ncol = m)
      }
      
      p_values_x <- numeric(m)
      p_values_z <- numeric(m)
      
      # loop over number of simulation replications
      for (k in 1:m)
      {
        # form data frame
        sim <- data.frame(x = x_mat[, k], y = y_mat[, k], z = z_mat[, k])
        
        # fit GAM and GLM under current response family
        if (response_family == "Gaussian")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = gaussian(link = "identity"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Binomial")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = binomial(link = "logit"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Poisson")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = poisson(link = "log"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Gamma")
        {
          sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                               family = Gamma(link = "log"),
                               method = lambda,
                               data = sim)
        }
        
        s_tab <- summary(sim.gam)$s.table
        p_values_x[k] <- s_tab["s(x)", "p-value"]
        p_values_z[k] <- s_tab["s(z)", "p-value"]
        
      } # END OF m SIMULATIONS LOOP
      
      p_values_x_list[[i]] <- p_values_x
      p_values_z_list[[i]] <- p_values_z
      
      row_df <- data.frame(
        family        = response_family,
        rho           = rho,
        F1x           = mean(p_values_x <= 0.01),
        F5x           = mean(p_values_x <= 0.05),
        F10x          = mean(p_values_x <= 0.1)
      )
      
      results_power22_list[[row_index]] <- row_df
      
    } # END OF RHO LOOP
    
    # export plot
    plot(ecdf(p_values_z_list[[1]]),
         col = "black",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = paste0(response_family),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    lines(ecdf(p_values_z_list[[2]]), col = "red")
    lines(ecdf(p_values_z_list[[3]]), col = "blue")
    lines(ecdf(p_values_z_list[[4]]), col = "darkgreen")
    
  } # END OF RESPONSE LOOP
  
  dev.off()
  
} # END OF SIMULATION 2.2



## 2.3 COMPARISON TO NULL SMOOTH -----------------------------------------------
{
  
  # parameters for current simulation
  set.seed(2300)
  m <- 5000
  n <- 500
  lambda <- "REML"
  twopi <- 2*pi
  rho <- 0.6
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  # plot 3
  # response_delta <- c(0.5, 1.0, 0.5, 0.5)
  # plot 2
  response_delta <- c(0.1, 0.4, 0.1, 0.1)
  c("black", "red", "blue", "darkgreen")
  
  results_power23_list <- vector("list", 4)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power_23_plot2.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
  h_mat <- scale(sin(twopi * z_mat))
  eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    delta <- response_delta[j]
    p_values_x <- numeric(m)
    p_values_z <- numeric(m)
    
    if (response_family == "Gaussian")
    {
      y_mat <- z_mat + delta * h_mat + eps_mat
    }
    
    if (response_family == "Binomial")
    {
      eta_mat <- 0.5 * z_mat + delta * h_mat
      q_mat <- plogis(eta_mat)
      y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
    }
    
    if (response_family == "Poisson")
    {
      eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
      lambda_mat <- exp(eta_mat)
      y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
    }
    
    if (response_family == "Gamma")
    {
      eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
      mu_mat <- exp(eta_mat)
      y_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_mat) / 5),
                      nrow = n, ncol = m)
    }
    
    # loop over number of simulation replications
    for (k in 1:m)
    {
      # form data frame
      sim <- data.frame(x = x_mat[, k], y = y_mat[, k], z = z_mat[, k])
      
      # fit GAM and GLM under current response family
      if (response_family == "Gaussian")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                             family = gaussian(link = "identity"),
                             method = lambda,
                             data = sim)
      }
      
      if (response_family == "Binomial")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                             family = binomial(link = "logit"),
                             method = lambda,
                             data = sim)
      }
      
      if (response_family == "Poisson")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                             family = poisson(link = "log"),
                             method = lambda,
                             data = sim)
      }
      
      if (response_family == "Gamma")
      {
        sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k = 10, m = c(2,0)),
                             family = Gamma(link = "log"),
                             method = lambda,
                             data = sim)
      }
      
      s_tab <- summary(sim.gam)$s.table
      p_values_x[k] <- s_tab["s(x)", "p-value"]
      p_values_z[k] <- s_tab["s(z)", "p-value"]
      
      
    } # END OF m SIMULATIONS LOOP
    
    row_df <- data.frame(
      family        = response_family,
      F1x           = mean(p_values_x <= 0.01),
      F5x           = mean(p_values_x <= 0.05),
      F10x          = mean(p_values_x <= 0.1)
    )
    
    results_power23_list[[j]] <- row_df
    
    # export plot
    plot(ecdf(p_values_z),
         col = "black",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = paste0(response_family),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    lines(ecdf(p_values_x), col = "red")
    abline(0, 1, lty = 2)
    
    
  } # END OF RESPONSE LOOP
  
  dev.off()
  
} # END OF SIMULATION 2.3


