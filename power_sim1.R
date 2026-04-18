# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# POWER SIMULATION SCENARIO 1

library(mgcv)
library(xtable)
twopi <- 2*pi

## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # parameters for current simulation
  set.seed(1100)
  m <- 5000
  n_values <- c(50,250,1000)
  lambda <- "REML"
  twopi <- 2*pi
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power_11_plot3.pdf", width = 12, height = 4)
  par(mfrow=c(1,3))
  
  for (i in 1:length(n_values))
  {
    # update to current sample size
    n <- n_values[i]
    
    # DGP for all responses
    x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
    h_mat <- scale(sin(twopi * x_mat))
    
    for (j in 1:length(response))
    {
      response_family <- response[j]
      
      # generate y under current family
      if (response_family == "Gaussian")
      {
        delta <- 0.05
        p_values_gaus <- numeric(m)
        eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
        y_mat <- delta * h_mat + eps_mat
      }
      
      if (response_family == "Binomial")
      {
        delta <- 0.20
        p_values_bin <- numeric(m)
        eta_mat <- delta * h_mat
        q_mat <- plogis(eta_mat)
        y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Poisson")
      {
        delta <- 0.04
        p_values_poi <- numeric(m)
        eta_mat <- log(5) + delta * h_mat
        lambda_mat <- exp(eta_mat)
        y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Gamma")
      {
        delta <- 0.04
        p_values_gamma <- numeric(m)
        eta_mat <- log(5) + delta * h_mat
        mu_mat <- exp(eta_mat)
        y_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_mat) / 5),
                        nrow = n, ncol = m)
      }
      
      
      # loop over number of simulation replications
      for (k in 1:m)
      {
        # form data frame
        sim <- data.frame(x = x_mat[, k], y = y_mat[, k])
        
        # fit GAM and GLM under current response family
        if (response_family == "Gaussian")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = gaussian(link = "identity"),
                               method = lambda,
                               data = sim)
          s_tab <- summary(sim.gam)$s.table
          p_values_gaus[k] <- s_tab["s(x)", "p-value"]
        }
        
        if (response_family == "Binomial")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = binomial(link = "logit"),
                               method = lambda,
                               data = sim)
          s_tab <- summary(sim.gam)$s.table
          p_values_bin[k] <- s_tab["s(x)", "p-value"]
        }
        
        if (response_family == "Poisson")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = poisson(link = "log"),
                               method = lambda,
                               data = sim)
          s_tab <- summary(sim.gam)$s.table
          p_values_poi[k] <- s_tab["s(x)", "p-value"]
        }
        
        if (response_family == "Gamma")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = Gamma(link = "log"),
                               method = lambda,
                               data = sim)
          s_tab <- summary(sim.gam)$s.table
          p_values_gamma[k] <- s_tab["s(x)", "p-value"]
        }
        
        
      } # END OF m SIMULATIONS LOOP
      
    } # END OF RESPONSE LOOP
    
    # export plot
    plot(ecdf(p_values_gaus),
         col = "black",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = paste0("n = ", n),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    
    lines(ecdf(p_values_bin), col="red")
    lines(ecdf(p_values_poi), col="blue")
    lines(ecdf(p_values_gamma), col="darkgreen")
    
  } # END OF SAMPLE SIZE LOOP
  
  dev.off()
  
} # END OF SIMULATION 1.1



## 1.2 EFFECT OF TRUE SMOOTH FUNCTION ------------------------------------------
{
  # parameters for current simulation
  set.seed(1200)
  m <- 5000
  n <- 500
  lambda <- "REML"
  twopi <- 2*pi
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  response_colour <- c("black", "red", "blue", "darkgreen")
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power_12_plot2.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  h_sin_mat <- scale(sin(twopi * x_mat))
  h_cubic_mat <- scale((x_mat - 0.5)^3)
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    
    # generate y under current family
    if (response_family == "Gaussian")
    {
      delta <- 0.05
      p_values_sin_gaus <- numeric(m)
      p_values_cubic_gaus <- numeric(m)
      eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
      y_sin_mat <- delta * h_sin_mat + eps_mat
      y_cubic_mat <- delta * h_cubic_mat + eps_mat
    }
    
    if (response_family == "Binomial")
    {
      delta <- 0.20
      p_values_sin_bin <- numeric(m)
      p_values_cubic_bin <- numeric(m)
      eta_sin_mat <- delta * h_sin_mat
      eta_cubic_mat <- delta * h_cubic_mat
      q_sin_mat <- plogis(eta_sin_mat)
      q_cubic_mat <- plogis(eta_cubic_mat)
      y_sin_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_sin_mat)), nrow = n, ncol = m)
      y_cubic_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_cubic_mat)), nrow = n, ncol = m)
    }
    
    if (response_family == "Poisson")
    {
      delta <- 0.04
      p_values_sin_poi <- numeric(m)
      p_values_cubic_poi <- numeric(m)
      eta_sin_mat <- log(5) + delta * h_sin_mat
      eta_cubic_mat <- log(5) + delta * h_cubic_mat
      lambda_sin_mat <- exp(eta_sin_mat)
      lambda_cubic_mat <- exp(eta_cubic_mat)
      y_sin_mat <- matrix(rpois(n * m, lambda = c(lambda_sin_mat)), nrow = n, ncol = m)
      y_cubic_mat <- matrix(rpois(n * m, lambda = c(lambda_cubic_mat)), nrow = n, ncol = m)
    }
    
    if (response_family == "Gamma")
    {
      delta <- 0.04
      p_values_sin_gamma <- numeric(m)
      p_values_cubic_gamma <- numeric(m)
      eta_sin_mat <- log(5) + delta * h_sin_mat
      eta_cubic_mat <- log(5) + delta * h_cubic_mat
      mu_sin_mat <- exp(eta_sin_mat)
      mu_cubic_mat <- exp(eta_cubic_mat)
      y_sin_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_sin_mat) / 5),
                          nrow = n, ncol = m)
      y_cubic_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_cubic_mat) / 5),
                            nrow = n, ncol = m)
    }
    
    # loop over number of simulation replications
    for (k in 1:m)
    {
      # form data frame
      sim <- data.frame(x = x_mat[, k], y_sin = y_sin_mat[, k], 
                        y_cubic = y_cubic_mat[, k])
      
      # fit GAM and GLM under current response family
      if (response_family == "Gaussian")
      {
        sim.gam <- mgcv::gam(y_sin ~ s(x, k=10),
                             family = gaussian(link = "identity"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_sin_gaus[k] <- s_tab["s(x)", "p-value"]
        
        sim.gam <- mgcv::gam(y_cubic ~ s(x, k=10),
                             family = gaussian(link = "identity"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_cubic_gaus[k] <- s_tab["s(x)", "p-value"]
      }
      
      if (response_family == "Binomial")
      {
        sim.gam <- mgcv::gam(y_sin ~ s(x, k=10),
                             family = binomial(link = "logit"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_sin_bin[k] <- s_tab["s(x)", "p-value"]
        
        sim.gam <- mgcv::gam(y_cubic ~ s(x, k=10),
                             family = binomial(link = "logit"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_cubic_bin[k] <- s_tab["s(x)", "p-value"]
      }
      
      if (response_family == "Poisson")
      {
        sim.gam <- mgcv::gam(y_sin ~ s(x, k=10),
                             family = poisson(link = "log"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_sin_poi[k] <- s_tab["s(x)", "p-value"]
        
        sim.gam <- mgcv::gam(y_cubic ~ s(x, k=10),
                             family = poisson(link = "log"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_cubic_poi[k] <- s_tab["s(x)", "p-value"]
      }
      
      if (response_family == "Gamma")
      {
        sim.gam <- mgcv::gam(y_sin ~ s(x, k=10),
                             family = Gamma(link = "log"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_sin_gamma[k] <- s_tab["s(x)", "p-value"]
        
        sim.gam <- mgcv::gam(y_cubic ~ s(x, k=10),
                             family = Gamma(link = "log"),
                             method = lambda,
                             data = sim)
        s_tab <- summary(sim.gam)$s.table
        p_values_cubic_gamma[k] <- s_tab["s(x)", "p-value"]
      }
      
    } # END OF m SIMULATIONS LOOP
    
  } # END OF RESPONSE LOOP
  
  # export plot
  plot(ecdf(p_values_sin_gaus),
       col = "black",
       lty = 1,
       xlim = c(0, 1),
       ylim = c(0, 1),
       main = response[1],
       xlab = "p-value", 
       ylab = "Empirical CDF")
  lines(ecdf(p_values_cubic_gaus), col = "red")
  
  plot(ecdf(p_values_sin_bin),
       col = "black",
       lty = 1,
       xlim = c(0, 1),
       ylim = c(0, 1),
       main = response[2],
       xlab = "p-value", 
       ylab = "Empirical CDF")
  lines(ecdf(p_values_cubic_bin), col = "red")
  
  plot(ecdf(p_values_sin_poi),
       col = "black",
       lty = 1,
       xlim = c(0, 1),
       ylim = c(0, 1),
       main = response[3],
       xlab = "p-value", 
       ylab = "Empirical CDF")
  lines(ecdf(p_values_cubic_poi), col = "red")
  
  plot(ecdf(p_values_sin_gamma),
       col = "black",
       xlim = c(0, 1),
       ylim = c(0, 1),
       main = response[4],
       xlab = "p-value", 
       ylab = "Empirical CDF")
  lines(ecdf(p_values_cubic_gamma), col = "red")
  
  dev.off()
  
} # END OF 1.2



## 1.3 EFFECT OF SIGNAL STRENGTH -----------------------------------------------
{
  # parameters for current simulation
  set.seed(1300)
  m <- 5000
  n <- 500
  lambda <- "REML"
  twopi <- 2*pi
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  response_colour <- c("black", "red", "blue", "darkgreen")
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power_13_plot2.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  h_mat <- scale(sin(twopi * x_mat))
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    
    # generate y under current family
    if (response_family == "Gaussian")
    {
      delta_values <- c(0.01, 0.05, 0.10, 0.20)
      eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
    }
    
    if (response_family == "Binomial")
    {
      delta_values <- c(0.10, 0.25, 0.50, 1.00)
    }
    
    if (response_family == "Poisson")
    {
      delta_values <- c(0.01, 0.05, 0.10, 0.20)
    }
    
    if (response_family == "Gamma")
    {
      delta_values <- c(0.01, 0.05, 0.10, 0.20)
    }
    
    p_values_list <- vector("list", length(delta_values))
    
    for (d in 1:length(delta_values))
    {
      delta <- delta_values[d]
      p_values <- numeric(m)
      
      if (response_family == "Gaussian")
      {
        y_mat <- delta * h_mat + eps_mat
      }
      
      if (response_family == "Binomial")
      {
        eta_mat <- delta * h_mat
        q_mat <- plogis(eta_mat)
        y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Poisson")
      {
        eta_mat <- log(5) + delta * h_mat
        lambda_mat <- exp(eta_mat)
        y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
      }
      
      if (response_family == "Gamma")
      {
        eta_mat <- log(5) + delta * h_mat
        mu_mat <- exp(eta_mat)
        y_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_mat) / 5),
                        nrow = n, ncol = m)
      }
      
      # loop over number of simulation replications
      for (k in 1:m)
      {
        # form data frame
        sim <- data.frame(x = x_mat[, k], y = y_mat[, k])
        
        # fit GAM and GLM under current response family
        if (response_family == "Gaussian")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = gaussian(link = "identity"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Binomial")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = binomial(link = "logit"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Poisson")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = poisson(link = "log"),
                               method = lambda,
                               data = sim)
        }
        
        if (response_family == "Gamma")
        {
          sim.gam <- mgcv::gam(y ~ s(x, k=10),
                               family = Gamma(link = "log"),
                               method = lambda,
                               data = sim)
        }
        
        s_tab <- summary(sim.gam)$s.table
        p_values[k] <- s_tab["s(x)", "p-value"]
        
        
      } # END OF m SIMULATIONS LOOP
      
      p_values_list[[d]] <- p_values
      
    } # END OF DELTA LOOP
    
    # export plot
    plot(ecdf(p_values_list[[1]]),
         col = "black",
         xlim = c(0, 1),
         ylim = c(0, 1),
         main = paste0(response_family),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    lines(ecdf(p_values_list[[2]]), col = "red")
    lines(ecdf(p_values_list[[3]]), col = "blue")
    lines(ecdf(p_values_list[[4]]), col = "darkgreen")
    
  } # END OF RESPONSE LOOP
  
  dev.off()
  
} # END OF SIMULATION 1.3


