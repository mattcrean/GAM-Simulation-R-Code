# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# SIMULATION SCENARIO 1
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
results_list_sp <- vector("list", 4)


# GAUSSIAN CASE ----------------------------------------------------------------
{
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(101)
  
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
      Response   = "Gaussian",
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 1.1

results_list_n[[1]] <- results_list

## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(102)
  
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
      Response   = "Gaussian",
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 1.2

results_list_sp[[1]] <- results_list
  
## 1.3 EFFECT OF BASIS DIMENSION -----------------------------------------------
{
  # set seed for reproducibility
  set.seed(103)
  
  # constant parameters
  m <- 5000
  n <- 250
  lambda <- "REML"
  
  # varying parameter
  k_values <- c(5,10,20,30)
  
  # prepare results vector
  results_list_k <- vector("list", length(k_values))
  
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
    results_list_k[[i]] <- row_df
  }
  
} # END OF SIMULATION 1.3
} # END OF GAUSSIAN CASE



# BINOMIAL CASE ----------------------------------------------------------------
{
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(201)
  
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
      Response   = "Binomial",
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 1.1
  
results_list_n[[2]] <- results_list
  
## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(202)
  
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
      Response   = "Binomial",
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }

} # END OF SIMULATION 1.2

results_list_sp[[2]] <- results_list

} # END OF BINOMIAL CASE



# POISSON CASE -----------------------------------------------------------------
{
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(301)
  
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
      Response   = "Poisson",
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 1.1
  
results_list_n[[3]] <- results_list
  
## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(302)
  
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
      Response   = "Poisson",
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 1.2

results_list_sp[[3]] <- results_list

} # END OF POISSON CASE



# GAMMA CASE -------------------------------------------------------------------
{
## 1.1 EFFECT OF SAMPLE SIZE ---------------------------------------------------
{
  # set seed for reproducibility
  set.seed(401)
  
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
      Response   = "Gamma",
      n          = n,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }
  
} # END OF SIMULATION 1.1
  
results_list_n[[4]] <- results_list
  
## 1.2 EFFECT OF SMOOTHING PARAMETER ESTIMATION METHOD -------------------------
{
  # set seed for reproducibility
  set.seed(402)
  
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
      Response   = "Gamma",
      method     = lambda,
      F1         = mean(p_values <= 0.01),
      F5         = mean(p_values <= 0.05),
      F10        = mean(p_values <= 0.1),
      edf_mean   = mean(edf),
      edf_sd     = sd(edf)
    )
    # write this row to corresponding element of list
    results_list[[i]] <- row_df
    
  }

} # END OF SIMULATION 1.2

results_list_sp[[4]] <- results_list

} # END OF GAMMA CASE



# GLM COMPARISON PLOT ----------------------------------------------------------
{
  # parameters for current simulation
  set.seed(501)
  m <- 5000
  lambda <- "REML"
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  n <- 500
  mu <- 5
  kappa <- 5
  poi_mean <- 5
  q <- 0.5
  sigma <- 0.5
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim1_glm_comp.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    
    p_values_gam <- numeric(m)
    p_values_glm <- numeric(m)
    
    # generate y under current family
    if (response_family == "Gaussian")
    {
      y_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
      fam <- gaussian(link = "identity")
    }
    
    if (response_family == "Binomial")
    {
      y_mat <- matrix(rbinom(n * m, size = 1, prob = q), nrow = n, ncol = m)
      fam <- binomial(link = "logit")
    }
    
    if (response_family == "Poisson")
    {
      y_mat <- matrix(rpois(n * m, lambda = poi_mean), nrow = n, ncol = m)
      fam <- poisson(link = "log")
    }
    
    if (response_family == "Gamma")
    {
      y_mat <- matrix(rgamma(n * m, shape = kappa, scale = mu / kappa),
                      nrow = n, ncol = m)
      fam <- Gamma(link = "log")
    }
    
    # loop over number of simulation replications
    for (k in 1:m)
    {
      # form data frame
      sim <- data.frame(x = x_mat[, k], y = y_mat[, k])
      
      sim.gam <- mgcv::gam(y ~ s(x, k=10),
                           family = fam,
                           method = lambda,
                           data = sim)
      sim.glm <- glm(y ~ x,
                     family = fam,
                     data = sim)
      
      p_values_gam[k] <- summary(sim.gam)$s.table["s(x)", "p-value"]
      p_values_glm[k] <- summary(sim.glm)$coefficients["x", 4]
      
    } # END OF m SIMULATIONS LOOP
    
    # export plot
    plot(ecdf(p_values_gam),
         col = "black",
         xlim = c(0, 0.1),
         ylim = c(0, 0.1),
         main = paste0(response_family),
         xlab = "p-value", 
         ylab = "Empirical CDF")
    lines(ecdf(p_values_glm), col="red")
    abline(0, 1, lty=2)
    lines(p_mcse_plot,upper,lty=3,col="grey40")
    lines(p_mcse_plot,lower,lty=3,col="grey40")
    
  } # END OF RESPONSE LOOP
  
  dev.off()
  
}  # END OF GLM COMPARISON



# SAMPLE SIZE PLOT -------------------------------------------------------------
{
  # line plot of empirical 5% rejection rate across response families
  # using summary values already obtained from the simulations
  
  # sample sizes and equally spaced plotting positions
  n_values <- c(25, 50, 100, 250, 500, 1000)
  x_pos <- 1:6
  
  # empirical rejection rates at alpha = 0.05
  gaussian <- c(0.075, 0.061, 0.060, 0.056, 0.056, 0.052)
  binomial <- c(0.019, 0.032, 0.046, 0.053, 0.052, 0.060)
  poisson  <- c(0.057, 0.051, 0.056, 0.058, 0.059, 0.051)
  gamma <- c(0.079, 0.067, 0.058, 0.061, 0.058, 0.055)
  
  # MCSE band for nominal 5% level
  m <- 5000
  alpha <- 0.05
  mcse <- sqrt(alpha * (1 - alpha) / m)
  upper <- alpha + 1.96 * mcse
  lower <- alpha - 1.96 * mcse
  
  # choose colours
  col_gaussian <- "black"
  col_binomial <- "red"
  col_poisson  <- "blue"
  col_gamma    <- "darkgreen"
  
  pdf("sim1_lineplot05.pdf", width = 8, height = 5)
  
  par(mar = c(5, 5, 3, 2))
  
  plot(x_pos, gaussian,
       type = "l",
       lwd = 2,
       xaxt = "n",
       col = col_gaussian,
       ylim = c(0.015, 0.08),
       xlab = expression("Sample size, " * n),
       ylab = expression(hat(F)[0.05]))
  axis(1, at = x_pos, labels = n_values)
  
  lines(x_pos, binomial, col=col_binomial, type = "l", lwd = 2)
  lines(x_pos, poisson, col=col_poisson, type = "l", lwd = 2)
  lines(x_pos, gamma, col=col_gamma, type = "l", lwd = 2)
  
  # nominal line
  abline(h = alpha, lty = 2, lwd = 1.2)
  
  # MCSE band
  abline(h = upper, lty = 3, lwd = 1, col = "grey40")
  abline(h = lower, lty = 3, lwd = 1, col = "grey40")
  
  legend("bottomright",
         legend = c("Gaussian", "Binomial", "Poisson", "Gamma"),
         lty = c(1, 1, 1, 1),
         lwd = c(2, 2, 2, 2),
         col = c(col_gaussian, col_binomial, col_poisson, col_gamma),
         bty = "n",
         cex = 0.85)
  
  dev.off()
} # END OF SAMPLE SIZE PLOT CODE


results_n <- do.call(
  rbind,
  lapply(results_list_n, function(x) do.call(rbind, x))
)
rownames(results_n) <- NULL
print(results_n)

writeLines("\n")

results_sp <- do.call(
  rbind,
  lapply(results_list_sp, function(x) do.call(rbind, x))
)
rownames(results_sp) <- NULL
print(results_sp)

writeLines("\n")

results_k <- do.call(rbind, results_list_k)
rownames(results_k) <- NULL
print(results_k)

}
