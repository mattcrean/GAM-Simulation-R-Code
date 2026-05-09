# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# POWER SIMULATION SCENARIO 2
{
library(mgcv)
twopi <- 2*pi


## 2.1 EFFECT OF SIGNAL STRENGTH RELATIVE TO NULL ------------------------------
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
  
  results_list <- vector("list", 16)
  row_index <- 0
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power2_delta.pdf", width = 8, height = 6.5)
  par(mfrow=c(2,2))
  
  # DGP for all responses
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
  h_mat <- scale(z_mat^2 - 1)
  eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
  
  for (j in 1:length(response))
  {
    response_family <- response[j]
    delta_values <- c(0.005, 0.01, 0.02, 0.04)
    
    if (response_family == "Binomial")
    {
      delta_values <- c(0.05, 0.10, 0.20, 0.30)
    }
    
    p_values_x_list <- vector("list", length(delta_values))
    p_values_z_list <- vector("list", length(delta_values))
    p_values_x_null_list <- vector("list", length(delta_values))
    p_values_z_null_list <- vector("list", length(delta_values))
    
    for (d in 1:length(delta_values))
    {
      row_index <- row_index + 1
      delta <- delta_values[d]
      p_values_x <- numeric(m)
      p_values_z <- numeric(m)
      p_values_x_null <- numeric(m)
      p_values_z_null <- numeric(m)
      
      if (response_family == "Gaussian")
      {
        y_mat <- z_mat + delta * h_mat + eps_mat
        y_mat_null <- z_mat + eps_mat
        
        fam <- gaussian(link = "identity")
      }
      
      if (response_family == "Binomial")
      {
        eta_mat <- 0.5 * z_mat + delta * h_mat
        q_mat <- plogis(eta_mat)
        y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
        y_mat_null <- matrix(rbinom(n * m, size = 1, prob = c(plogis(0.5 * z_mat))), nrow = n, ncol = m)
        
        fam <- binomial(link = "logit")
      }
      
      if (response_family == "Poisson")
      {
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        lambda_mat <- exp(eta_mat)
        y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
        y_mat_null <- matrix(rpois(n * m, lambda = c(exp(log(5) + 0.5 * z_mat))), nrow = n, ncol = m)
        
        fam <- poisson(link = "log")
      }
      
      if (response_family == "Gamma")
      {
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        mu_mat <- exp(eta_mat)
        y_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_mat) / 5),
                        nrow = n, ncol = m)
        y_mat_null <- matrix(rgamma(n * m, shape = 5, scale = c(exp(log(5) + 0.5 * z_mat)) / 5),
                        nrow = n, ncol = m)
        
        fam <- Gamma(link = "log")
      }
      
      # loop over number of simulation replications
      for (k in 1:m)
      {
        # form data frame
        sim <- data.frame(x = x_mat[, k], y = y_mat[, k], z = z_mat[, k], y_null = y_mat_null[, k])
        
        sim.gam <- mgcv::gam(y ~ z + s(x, k = 10) + s(z, k = 10, m = c(2, 0)),
                            family = fam,
                            method = lambda,
                            data = sim
        )
        
        sim.gam_null <- mgcv::gam(y_null ~ z + s(x, k = 10) + s(z, k = 10, m = c(2, 0)),
                             family = fam,
                             method = lambda,
                             data = sim
        )
        
        s_tab <- summary(sim.gam)$s.table
        p_values_x[k] <- s_tab["s(x)", "p-value"]
        p_values_z[k] <- s_tab["s(z)", "p-value"]
        p_values_x_null[k] <- summary(sim.gam_null)$s.table["s(x)", "p-value"]
        p_values_z_null[k] <- summary(sim.gam_null)$s.table["s(z)", "p-value"]
        
        
      } # END OF m SIMULATIONS LOOP
      
      p_values_x_list[[d]] <- p_values_x
      p_values_z_list[[d]] <- p_values_z
      p_values_x_null_list[[d]] <- p_values_x_null
      p_values_z_null_list[[d]] <- p_values_z_null
      
      row_df <- data.frame(
        family        = response_family,
        delta         = delta,
        type          = c("Null", "Alternative"),
        F1z           = c(mean(p_values_z_null <= 0.01), mean(p_values_z <= 0.01)),
        F5z           = c(mean(p_values_z_null <= 0.05), mean(p_values_z <= 0.05)),
        F10z          = c(mean(p_values_z_null <= 0.10), mean(p_values_z <= 0.10)),
        F1x           = c(mean(p_values_x_null <= 0.01), mean(p_values_x <= 0.01)),
        F5x           = c(mean(p_values_x_null <= 0.05), mean(p_values_x <= 0.05)),
        F10x          = c(mean(p_values_x_null <= 0.10), mean(p_values_x <= 0.1))
      )
      
      results_list[[row_index]] <- row_df
      
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
  rho_values <- c(0,0.3,0.6,0.9)
  response <- c("Gaussian", "Binomial", "Poisson", "Gamma")
  response_colour <- c("black", "red", "blue", "darkgreen")
  
  row_index <- 0
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim_power2_rho.pdf", width = 8, height = 6.5)
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
      h_mat <- scale(z_mat^2 - 1)
      row_index <- row_index + 1
      
      if (response_family == "Gaussian")
      {
        delta <- 0.02
        eps_mat <- matrix(rnorm(n * m, mean = 0, sd = 0.5), nrow = n, ncol = m)
        y_mat <- z_mat + delta * h_mat + eps_mat
        fam <- gaussian(link = "identity")
      }
      
      if (response_family == "Binomial")
      {
        delta <- 0.20
        eta_mat <- 0.5 * z_mat + delta * h_mat
        q_mat <- plogis(eta_mat)
        y_mat <- matrix(rbinom(n * m, size = 1, prob = c(q_mat)), nrow = n, ncol = m)
        fam <- binomial(link = "logit")
      }
      
      if (response_family == "Poisson")
      {
        delta <- 0.02
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        lambda_mat <- exp(eta_mat)
        y_mat <- matrix(rpois(n * m, lambda = c(lambda_mat)), nrow = n, ncol = m)
        fam <- poisson(link = "log")
      }
      
      if (response_family == "Gamma")
      {
        delta <- 0.02
        eta_mat <- log(5) + 0.5 * z_mat + delta * h_mat
        mu_mat <- exp(eta_mat)
        y_mat <- matrix(rgamma(n * m, shape = 5, scale = c(mu_mat) / 5),
                        nrow = n, ncol = m)
        fam <- Gamma(link = "log")
      }
      
      p_values_x <- numeric(m)
      p_values_z <- numeric(m)
      
      # loop over number of simulation replications
      for (k in 1:m)
      {
        # form data frame
        sim <- data.frame(x = x_mat[, k], y = y_mat[, k], z = z_mat[, k])
        
        # fit GAM under current response family
        sim.gam <- mgcv::gam(y ~ z + s(x, k = 10) + s(z, k = 10, m = c(2, 0)),
                             family = fam,
                             method = lambda,
                             data = sim
        )
        
        s_tab <- summary(sim.gam)$s.table
        p_values_x[k] <- s_tab["s(x)", "p-value"]
        p_values_z[k] <- s_tab["s(z)", "p-value"]
        
      } # END OF m SIMULATIONS LOOP
      
      p_values_x_list[[i]] <- p_values_x
      p_values_z_list[[i]] <- p_values_z
      
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


results_delta <- do.call(rbind, results_list)
rownames(results_delta) <- NULL
print(results_delta)

}
