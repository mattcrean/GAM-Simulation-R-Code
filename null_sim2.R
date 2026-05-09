# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# SIMULATION SCENARIO 2
{
library(mgcv)

# MCSE bands for all plots (given m=5000)
{
  m <- 5000
  p_mcse_plot <- seq(0,0.15,length=200)
  mcse <- sqrt(p_mcse_plot*(1-p_mcse_plot)/m)
  upper <- p_mcse_plot + 1.96*mcse
  lower <- p_mcse_plot - 1.96*mcse
}

results_list_n <- vector("list", 4)

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
      Response   = "Gaussian",
      n          = n,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 2.1

results_list_n[[1]] <- results_list
  
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
  results_list_gaussian <- vector("list", length(sigma_values))
  
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
      Response   = "Gaussian",
      sigma      = sigma,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list_gaussian[[i]] <- row_df
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
      Response   = "Binomial",
      n          = n,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 2.1

results_list_n[[2]] <- results_list
  
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
  results_list_binomial <- vector("list", length(success_values))
  
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
      Response   = "Binomial",
      Success    = success,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list_binomial[[i]] <- row_df
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
      Response   = "Poisson",
      n          = n,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 2.1
  
results_list_n[[3]] <- results_list

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
  results_list_poisson <- vector("list", length(lambda_values))
  
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
      Response   = "Poisson",
      Mean       = lam,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list_poisson[[i]] <- row_df
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
      Response   = "Gamma",
      n          = n,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 2.1

results_list_n[[4]] <- results_list
  
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
  results_list_gamma <- vector("list", length(mu_values))
  
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
    
    # DGP for current Gamma mean
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
      Response   = "Gamma",
      Mean       = mu,
      Smooth     = c("f1(x)", "f2(z)"),
      F1         = c(mean(p_values_x <= 0.01), mean(p_values_z <= 0.01)),
      F5         = c(mean(p_values_x <= 0.05), mean(p_values_z <= 0.05)),
      F10        = c(mean(p_values_x <= 0.10), mean(p_values_z <= 0.10)),
      edf_mean   = c(mean(edf_x), mean(edf_z)),
      edf_sd     = c(sd(edf_x), sd(edf_z))
    )
    # write this row to corresponding element of list
    results_list_gamma[[i]] <- row_df
  }
} # END OF SIMULATION 2.2

} # END OF GAMMA SIMULATION



# ONE VS TWO FITTED NULL SMOOTHS -----------------------------------------------
{
  # parameters for current simulation
  set.seed(511)
  m <- 5000
  n <- 500
  lambda <- "REML"
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim2_null_comp.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  z_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    p_values_x_1 <- numeric(m)
    p_values_x_2 <- numeric(m)
    
    if (response_family == "Gaussian")
    {
      y_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
    }
    
    if (response_family == "Binomial")
    {
      y_mat <- matrix(rbinom(n * m, size = 1, prob = 0.5), nrow = n, ncol = m)
    }
    
    if (response_family == "Poisson")
    {
      y_mat <- matrix(rpois(n * m, lambda = 5), nrow = n, ncol = m)
    }
    
    if (response_family == "Gamma")
    {
      y_mat <- matrix(rgamma(n * m, shape = 5, scale = 5 / 5),
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
        sim.gam1 <- mgcv::gam(y ~ s(x, k=10),
                             family = gaussian(link = "identity"),
                             method = lambda,
                             data = sim)
        sim.gam2 <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10),
                             family = gaussian(link = "identity"),
                             method = lambda,
                             data = sim)
      }
      
      if (response_family == "Binomial")
      {
        sim.gam1 <- mgcv::gam(y ~ s(x, k=10),
                             family = binomial(link = "logit"),
                             method = lambda,
                             data = sim)
        sim.gam2 <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10),
                             family = binomial(link = "logit"),
                             method = lambda,
                             data = sim)
      }
      
      if (response_family == "Poisson")
      {
        sim.gam1 <- mgcv::gam(y ~ s(x, k=10),
                             family = poisson(link = "log"),
                             method = lambda,
                             data = sim)
        sim.gam2 <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10),
                             family = poisson(link = "log"),
                             method = lambda,
                             data = sim)
      }
      
      if (response_family == "Gamma")
      {
        sim.gam1 <- mgcv::gam(y ~ s(x, k=10),
                             family = Gamma(link = "log"),
                             method = lambda,
                             data = sim)
        sim.gam2 <- mgcv::gam(y ~ s(x, k=10) + s(z, k=10),
                             family = Gamma(link = "log"),
                             method = lambda,
                             data = sim)
      }
      
      p_values_x_1[k] <- summary(sim.gam1)$s.table["s(x)", "p-value"]
      p_values_x_2[k] <- summary(sim.gam2)$s.table["s(x)", "p-value"]
      
        
      } # END OF m SIMULATIONS LOOP
      
    # export plot
    plot(ecdf(p_values_x_1),
         col = "black",
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0(response_family),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    lines(ecdf(p_values_x_2), col = "red")
    
    # reference line of the CDF of U(0,1)
    abline(0, 1, lty = 2)
    
    # MCSE bands
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  } # END OF RESPONSE LOOP
  
  dev.off()
  
} # END OF ONE VS TWO NULL COMPARISON



# RESPONSE PARAMETER PLOT ------------------------------------------------------
{
  # 2x2 grid of plots for f1
  # Each panel shows empirical rejection rates at alpha = 0.05 and 0.10
  # MCSE bands
  m <- 5000
  alpha05 <- 0.05
  alpha10 <- 0.10
  
  mcse05 <- sqrt(alpha05 * (1 - alpha05) / m)
  mcse10 <- sqrt(alpha10 * (1 - alpha10) / m)
  
  upper05 <- alpha05 + 1.96 * mcse05
  lower05 <- alpha05 - 1.96 * mcse05
  
  upper10 <- alpha10 + 1.96 * mcse10
  lower10 <- alpha10 - 1.96 * mcse10
  
  # line colours
  col_05 <- "black"
  col_10 <- "red"
  
  pdf("sim2_varf1.pdf", width = 9, height = 7)
  
  par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 1, 0))
  
  # GAUSSIAN
  sigma_vals <- c(0.25, 0.50, 1.00, 2.00, 5.00)
  gauss_05 <- c(0.055, 0.059, 0.054, 0.051, 0.047)
  gauss_10 <- c(0.108, 0.118, 0.107, 0.105, 0.105)
  
  plot(sigma_vals, gauss_05,
       type = "l", pch = 16, lwd = 2,
       col = col_05,
       ylim = c(0.04, 0.125),
       xlab = expression("Error s.d., " * sigma),
       ylab = expression(hat(F)),
       main = "Gaussian")
  
  lines(sigma_vals, gauss_10, type = "l", pch = 17, lwd = 2, col = col_10)
  
  abline(h = alpha05, lty = 2, lwd = 1.2)
  abline(h = alpha10, lty = 2, lwd = 1.2)
  abline(h = upper05, lty = 3, col = "grey40")
  abline(h = lower05, lty = 3, col = "grey40")
  abline(h = upper10, lty = 3, col = "grey40")
  abline(h = lower10, lty = 3, col = "grey40")
  
  # BINOMIAL
  q_vals <- c(0.3, 0.5, 0.7)
  binom_05 <- c(0.055, 0.051, 0.051)
  binom_10 <- c(0.119, 0.110, 0.115)
  
  plot(q_vals, binom_05,
       type = "l", pch = 16, lwd = 2,
       col = col_05,
       ylim = c(0.04, 0.125),
       xlab = expression("Success probability, " * q),
       ylab = expression(hat(F)),
       main = "Binomial")
  
  lines(q_vals, binom_10, type = "l", pch = 17, lwd = 2, col = col_10)
  
  abline(h = alpha05, lty = 2, lwd = 1.2)
  abline(h = alpha10, lty = 2, lwd = 1.2)
  abline(h = upper05, lty = 3, col = "grey40")
  abline(h = lower05, lty = 3, col = "grey40")
  abline(h = upper10, lty = 3, col = "grey40")
  abline(h = lower10, lty = 3, col = "grey40")
  
  # POISSON
  lambda_vals <- c(3, 5, 7, 10)
  pois_05 <- c(0.052, 0.051, 0.058, 0.054)
  pois_10 <- c(0.114, 0.103, 0.117, 0.106)
  
  plot(lambda_vals, pois_05,
       type = "l", pch = 16, lwd = 2,
       col = col_05,
       ylim = c(0.04, 0.125),
       xlab = expression("Poisson mean, " * lambda),
       ylab = expression(hat(F)),
       main = "Poisson")
  
  lines(lambda_vals, pois_10, type = "l", pch = 17, lwd = 2, col = col_10)
  
  abline(h = alpha05, lty = 2, lwd = 1.2)
  abline(h = alpha10, lty = 2, lwd = 1.2)
  abline(h = upper05, lty = 3, col = "grey40")
  abline(h = lower05, lty = 3, col = "grey40")
  abline(h = upper10, lty = 3, col = "grey40")
  abline(h = lower10, lty = 3, col = "grey40")
  
  # GAMMA
  mu_vals <- c(3, 5, 7, 10)
  gamma_05 <- c(0.053, 0.052, 0.051, 0.055)
  gamma_10 <- c(0.111, 0.110, 0.105, 0.112)
  
  plot(mu_vals, gamma_05,
       type = "l", pch = 16, lwd = 2,
       col = col_05,
       ylim = c(0.04, 0.125),
       xlab = expression("Gamma mean, " * mu),
       ylab = expression(hat(F)),
       main = "Gamma")
  
  lines(mu_vals, gamma_10, type = "l", pch = 17, lwd = 2, col = col_10)
  
  abline(h = alpha05, lty = 2, lwd = 1.2)
  abline(h = alpha10, lty = 2, lwd = 1.2)
  abline(h = upper05, lty = 3, col = "grey40")
  abline(h = lower05, lty = 3, col = "grey40")
  abline(h = upper10, lty = 3, col = "grey40")
  abline(h = lower10, lty = 3, col = "grey40")
  
  dev.off()
}


results_n <- do.call(
  rbind,
  lapply(results_list_n, function(x) do.call(rbind, x))
)
rownames(results_n) <- NULL
print(results_n)

writeLines("\n")

results_gaussian <- do.call(rbind, results_list_gaussian)
rownames(results_gaussian) <- NULL
print(results_gaussian)

writeLines("\n")

results_binomial <- do.call(rbind, results_list_binomial)
rownames(results_binomial) <- NULL
print(results_binomial)

writeLines("\n")

results_poisson <- do.call(rbind, results_list_poisson)
rownames(results_poisson) <- NULL
print(results_poisson)

writeLines("\n")

results_gamma <- do.call(rbind, results_list_gamma)
rownames(results_gamma) <- NULL
print(results_gamma)

}
