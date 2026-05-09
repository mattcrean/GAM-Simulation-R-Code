# EVALUATING SIGNIFICANCE TESTS IN GENERALISED ADDITIVE MODELS -----------------
# SIMULATION SCENARIO 4 --------------------------------------------------------
# PERMUTATION TEST (GAUSSIAN ONLY)
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
  
  # set seed for reproducibility
  set.seed(401)
  
  # fixed parameters
  m <- 5000
  n <- 500
  rho <- 0.6
  sigma <- 0.5
  no_perm <- 499
  beta_true <- 1
  
  # DGP
  x_mat <- matrix(runif(n * m, min = 0, max = 1), nrow = n, ncol = m)
  u_mat <- matrix(rnorm(n * m, mean = 0, sd = 1), nrow = n, ncol = m)
  g_mat <- scale(sin(twopi * x_mat))
  z_mat <- scale(rho * g_mat + sqrt(1 - rho^2) * u_mat)
  eps_mat <- matrix(rnorm(n * m, mean = 0, sd = sigma), nrow = n, ncol = m)
  y_mat <- beta_true * z_mat + eps_mat
  
  # create empty vectors for outputs
  p_values_z <- numeric(m)
  p_perm_z <- numeric(m)
  
  # loop over number of simulation replications
  for (j in 1:m)
  {
    
    if (j %% 50 == 0) print(j)
    
    # form data frame
    sim <- data.frame(x = x_mat[, j], z = z_mat[, j], y = y_mat[, j])
    
    # GAM with null smooth term in z still present
    sim.gam <- mgcv::gam(y ~ z + s(x, k=10) + s(z, k=10, m=c(2,0)),
                         family = gaussian(link = "identity"),
                         method = "REML",
                         data = sim)
    
    # extract GAM smooth p-values
    s_tab <- summary(sim.gam)$s.table
    p_values_z[j] <- s_tab["s(z)", "p-value"]
    F_obs <- s_tab["s(z)", "F"]
    
    # GAM with the null smooth term in z removed
    sim.gam_null <- mgcv::gam(y ~ z + s(x, k=10),
                         family = gaussian(link = "identity"),
                         method = "REML",
                         data = sim)
    
    # fitted values and residuals under the null
    y_fitted_null <- fitted(sim.gam_null)
    residual_null <- residuals(sim.gam_null, type = "response")
    
    # prepare vector of permutation F test statistics
    F_perm <- numeric(no_perm)
    
    # conduct residual permutation test under fitted null model
    for (b in 1:no_perm)
    {
      sim_perm <- sim
      
      # permute the residuals under null
      sim_perm$y_perm <- y_fitted_null + sample(residual_null, replace=FALSE)
      
      # refit full model to permuted response
      perm.gam <- mgcv::gam(y_perm ~ z + s(x, k = 10) + s(z, k = 10, m = c(2, 0)),
                            family = gaussian(link = "identity"),
                            method = "REML",
                            data = sim_perm)
      
      s_tab <- summary(perm.gam)$s.table
      F_perm[b] <- s_tab["s(z)", "F"]
    }
    
    p_perm_z[j] <- (1 + sum(F_perm >= F_obs)) / (no_perm + 1)
    
  }
  
  # write rows of data frame for current family
  row_df <- data.frame(
    Model         = c("GAM", "PERM"),
    F1z           = c(mean(p_values_z <= 0.01), mean(p_perm_z <= 0.01)),
    F5z           = c(mean(p_values_z <= 0.05), mean(p_perm_z <= 0.05)),
    F10z          = c(mean(p_values_z <= 0.10), mean(p_perm_z <= 0.10))
  )
  
  print(row_df)
  
  # prepare plot grid and pdf file to be saved to
  pdf("sim4_residual_permutation.pdf", width = 8, height = 6.5)
  par(mfrow=c(1,1))
  
  # export plot
  plot(ecdf(p_perm_z),
       col = "black",
       xlim = c(0, 0.1),
       ylim = c(0, 0.1),
       do.points = FALSE,
       verticals = TRUE,
       main = paste0("Gaussian"),
       xlab = "p-value", 
       ylab = "Empirical CDF")
  lines(ecdf(p_values_z), col="red")
  abline(0, 1, lty=2)
  lines(p_mcse_plot,upper,lty=3,col="grey40")
  lines(p_mcse_plot,lower,lty=3,col="grey40")
  dev.off()
  
}
