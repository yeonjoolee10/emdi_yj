#' Two-fold Fay-Herriot Models for Indicators at a disaggregated regional level
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the linear mixed regression model with the
#' dependent variable on the left of a ~ operator and the explanatory
#' variables on the right, separated by + operators.
#' @param vardir a character string indicating the name of the variable
#' containing the domain-specific sampling variances of the direct estimates
#' that are included in \cr \code{combined_data}.
#' @param combined_data a data set containing all the input variables that are
#' needed for the estimation of the two-fold Fay-Herriot model: the direct estimates,
#' the sampling variances, the explanatory variables, the domains, the subdomains,
#' and the population size. In addition, the effective sample size needs to be
#' included, if the arcsin transformation is chosen.
#' @param domains a character string indicating the domain variable that is
#' included in \code{combined_data}.
#' @param subdomains a character string indicating the subdomain variable that is
#' included in \code{combined_data}.
#' @param subdomsize a character string indicating the variable which includes
#' the population size at subdomain level in \code{combined_data}.
#' @param transformation a character that determines the type of transformation
#' of the dependent variable and of the sampling variances. Methods that can be
#' chosen (i) no transformation ("\code{no}"), (ii) log transformation
#' ("\code{log}") of the dependent variable and of the sampling variances,
#' (iii) arcsin transformation ("\code{arcsin}") of the dependent variable and
#' of the sampling variances following. Defaults to "\code{no}".
#' @param eff_smpsize a character string indicating the name of the variable
#' containing the effective sample sizes that are included in
#' \code{combined_data}. Required argument when the arcsin transformation is
#' chosen. Defaults to \code{NULL}.
#' @param maxit a number determining the maximum number of iterations for the
#' variance components estimation of the random effects. Defaults to 50.
#' @param MSE if \code{TRUE}, MSE estimates are calculated based on parametric
#' bootstrap method. Defaults to \code{FALSE}.
#' @param B a single number which defines the number of bootstrap iterations
#' for the MSE estimation.
#' @param seed an integer to set the seed for the random number generator. For
#' the usage of random number generation see details. Defaults to \code{123}.
#' @return An object of class "fh_tf", "emdi" that provides estimators
#' for regional indicators at a disaggregated level like means and ratios and optionally
#' corresponding MSE estimates.
#' @references
#' @examples
#' @export
#' @importFrom magic adiag
#' @importFrom formula.tools get.vars
library(magic)
# library(MASS)
library(formula.tools)
################################################################################
fixed = fitopt_mean_posto$fixed
vardir = "Var_Mean_Pois"
domains = "region2"
subdomains = "ID_Posto"
transformation = "log"
subdomsize = "N_Post"
MSE = T
B = 10
combined_data = agg_data_posto_wo_bld
seed = 123
eff_smpsize = NA
maxit = 30
################################################################################

fh_tf <- function(fixed, vardir, combined_data, domains, subdomains,
                    subdomsize, transformation = "no", maxit = 50,
                    eff_smpsize = NA, MSE = FALSE, B = 50, seed = 123){

  point_estim <- point_fh_tf(fixed = fixed, vardir = vardir, domains = domains,
                          subdomains = subdomains, nunits = subdomsize,
                          trafo = transformation, maxit = maxit,
                          eff_n = eff_smpsize, data = combined_data)

  if(MSE == TRUE){
    # MSE estimation based on Gonzalez-Manteiga et al.
    smp_hat_varu <- point_estim$varu
    smp_hat_varv <- point_estim$varv
    smp_hat_vare <- point_estim$data[, c(domains, subdomains, "vare")]
    smp_beta <- point_estim$beta

    set.seed(seed)
    mse_estim <- mse_fh_tf(fixed = fixed, vardir = vardir, domains = domains,
                           subdomains = subdomains, nunits = subdomsize,
                           trafo = transformation, seed = seed, maxit = maxit,
                           eff_n = eff_smpsize, B = B,
                           smp_hat_varu = smp_hat_varu,
                           smp_hat_varv = smp_hat_varu,
                           smp_hat_vare = smp_hat_vare, smp_beta = smp_beta,
                           data = combined_data)
  }

  if (MSE == FALSE){
    return(list(EBLUP_Area = point_estim$EBLUP_Area,
                EBLUP_SubArea = point_estim$EBLUP_SubArea,
                fixed = fixed,
                domains = domains,
                subdomains = subdomains,
                vardir = vardir,
                varu = point_estim$varu,
                varv = point_estim$varv,
                beta = point_estim$beta,
                v_tilde = point_estim$v_tilde,
                u_tilde = point_estim$u_tilde,
                resid = point_estim$resid,
                R2 = point_estim$R2,
                data = point_estim$data))
  } else if(MSE == TRUE){
    return(list(EBLUP_Area = point_estim$EBLUP_Area,
                EBLUP_SubArea = point_estim$EBLUP_SubArea,
                MSE_Area = mse_estim$MSE_Area,
                MSE_SubArea = mse_estim$MSE_SubArea,
                fixed = fixed,
                domains = domains,
                subdomains = subdomains,
                vardir = vardir,
                varu = point_estim$varu,
                varv = point_estim$varv,
                beta = point_estim$beta,
                v_tilde = point_estim$v_tilde,
                u_tilde = point_estim$u_tilde,
                resid = point_estim$resid,
                R2 = point_estim$R2,
                data = point_estim$data))
  }
}
################################################################################
################################################################################
point_fh_tf <- function(fixed, vardir, domains, subdomains, nunits, trafo,
                        maxit, eff_n, data){
  # Framework-------------------------------------------------------------------
  # Number of units in each sub-area (in population)
  N_ij <- data[[nunits]]
  # Number of units in each area (in population)
  N_i <- aggregate(data[[nunits]], by = list(data[[domains]]), FUN = sum)$x
  # Total number of units (in population)
  N <- sum(N_ij)
  # Names of area
  DomName <- unique(data[[domains]])

  # data of sampled sub-areas-----------------------------------------------------
  # Observation of sub-area
  data$ObsSub <- ifelse(is.na(data[[fixed[[2]]]]), "no", "yes")

  # data of sampled sub-areas
  smp_data <- data[data$ObsSub == "yes", ]

  X <- model.matrix(object = fixed, data = smp_data)
  rownames(X) <- smp_data[[subdomains]]

  # Dir_y and vare at subarea level-----------------------------------------------
  # With arcsin trafo
  if(trafo == "arcsin"){
    y_ij_orig <- smp_data[[fixed[[2]]]]
    vare_orig <- smp_data[[vardir]]
    y_ij <- asin(sqrt(y_ij_orig))
    vare <- 1 / (4 * smp_data[[eff_n]])
    data$vare <- 1 / (4 * data[[eff_n]])
    smp_data[[fixed[[2]]]] <- y_ij
    smp_data$y_orig <- y_ij_orig
  } else if(trafo == "log") {
    y_ij_orig <- smp_data[[fixed[[2]]]]
    vare_orig <- smp_data[[vardir]]
    vare <- (y_ij_orig)^(-2) * vare_orig
    data$vare <- (1 / data[[fixed[[2]]]])^2 * data[[vardir]] # for bc back-trafo
    y_ij <- log(y_ij_orig)
    smp_data[[fixed[[2]]]] <- y_ij
    smp_data$y_orig <- y_ij_orig
  } else if(trafo == "no") {
    y_ij <- smp_data[[fixed[[2]]]]
    vare <- smp_data[[vardir]]
    data$vare <- data[[vardir]]
  }

  names(vare) <- smp_data[[subdomains]]
  names(y_ij) <- smp_data[[subdomains]]
  # Number of sub-areas in each area
  if(is.factor(smp_data[[domains]])){no <- as.vector(table(as.character(smp_data[[domains]])))}
  else {no <- as.vector(table(smp_data[[domains]]))}
  # Total number of areas in the sample
  m <- length(unique(smp_data[[domains]]))
  # Total number of areas in pop
  M <- length(unique(data[[domains]]))
  # Total number of subareas in the sample
  n <- length(smp_data[[subdomains]])
  # Number of x variables
  p <- ncol(X)

  # Data aggregation at area level (sampled data)
  agg_data <- aggregate(smp_data[, get.vars(fixed)],
                        by = list(smp_data[[domains]]), FUN = mean, na.rm=TRUE)
  y_bar <- agg_data[[fixed[[2]]]]
  names(y_bar) <- agg_data[, "Group.1"]
  X_bar <- model.matrix(object = fixed, data = agg_data)
  rownames(X_bar) <- agg_data[, "Group.1"]


  # Estimate var_v and var_u (based on sampled data)----------------------------
  varv0 <- var(smp_data[[fixed[[2]]]])
  varu0 <- var(smp_data[[fixed[[2]]]])

  dif.ml <- ml_NR_diff(varv = varv0, varu = varu0, vare = vare,
                       yranr = y_ij,
                       no = no, X = X, m = m)
  varv1 <- varv0
  varu1 <- varu0

  iter <- 0
  while((abs(dif.ml[1]) > 0.0001 | abs(dif.ml[2]) > 0.0001) & iter < maxit){
    iter <- iter + 1
    varv1 <- varv1+dif.ml[1]
    varu1 <- varu1+dif.ml[2]
    dif.ml <- ml_NR_diff(varv = varv1, varu = varu1, vare = vare,
                         yranr = y_ij,
                         no = no, X = X, m = m)
  }

  if(iter >= maxit) {
    warning(strwrap(prefix = " ", initial = "",
    paste0("There is a convergence problem in the variance components estimation process.
           You may try to increase the number of the maximum interation.")))}

  varv <- varv1
  varu <- varu1
  if(varv <= 0){
    varv <- 0.00001
    warning(strwrap(prefix = " ", initial = "",
                    paste0("The estimated variance of the domain level random
                           effect is negative. Therefore, it is replaced by 0.0001")))}

  if(varu <= 0){
    varu <- 0.00001
    warning(strwrap(prefix = " ", initial = "",
                    paste0("The estimated variance of the subdomain level random
                           effect is negative. Therefore, it is replaced by 0.00001")))
  }

  beta <- beta_tilde(varv = varv, varu = varu, vare = vare,
                     yranr = y_ij, no = no, m = m,
                     X = X)
  beta_tild <- beta$beta

  # Obtain EBLUP at subarea level
  gamma_ij <- varu/(varu + vare)
  gamma1_i <- numeric()
  y_bar_igamma_prep <- numeric()
  X_bar_igamma_prep <- matrix(ncol = p, nrow = m)

  for(i in 1:m){
    if (i == 1){
      a <- 1
      b <- no[1]
    } else {
      a <- a + no[i - 1]
      b <- b + no[i]
    }
    gamma1_i[i] <- sum(gamma_ij[a:b])
    y_bar_igamma_prep[i] <- sum((gamma_ij * y_ij)[a:b])
    if (p == 1) {X_bar_igamma_prep[i, ] <- sum((gamma_ij * X)[a:b, ])}
    else{
      if(b == a) {X_bar_igamma_prep[i, ] <- ((gamma_ij * X)[a:b, ])}
      else{X_bar_igamma_prep[i, ] <- colSums((gamma_ij * X)[a:b, ])}
    }
  }
  names(gamma1_i) <- unique(smp_data[[domains]])
  gamma2_i <- varv/(varv + (varu/gamma1_i))

  y_bar_igamma <- (1/gamma1_i)*y_bar_igamma_prep
  X_bar_igamma <- (1/gamma1_i)*X_bar_igamma_prep
  rownames(X_bar_igamma) <- names(gamma1_i)
  v_tilde_i <- gamma2_i * (y_bar_igamma - X_bar_igamma %*% beta_tild)

  v_tilde_i_sub <- rep(v_tilde_i, times = no)
  gamma2_i_sub <- rep(gamma2_i, times = no)
  y_bar_igamma_sub <- rep(y_bar_igamma, times = no)

  X_bar_igamma_sub <- do.call("rbind", replicate(no[1], X_bar_igamma[1, ], simplify = FALSE))
  for (i in 2:m) {
    X_bar_igamma_sub_i <- do.call("rbind", replicate(no[i], X_bar_igamma[i, ], simplify = FALSE))
    X_bar_igamma_sub <- rbind(X_bar_igamma_sub, X_bar_igamma_sub_i)
  }

  u_tilde_ij <- gamma_ij * (y_ij - X %*% beta_tild) - gamma2_i_sub * gamma_ij *
    (y_bar_igamma_sub - X_bar_igamma_sub %*% beta_tild)

  # EBLUP for mu_ij for sampled sub-areas
  mu_ij_smp <- X %*% beta_tild + v_tilde_i_sub + u_tilde_ij

  # EBLUP for mu_ij for non-sampled sub-areas & in sample areas
  no_smp_data <- data[data$ObsSub == "no", ]
  no_dom <- DomName[which(!is.element(DomName, unique(smp_data[[domains]])))]
  vec <- which(is.element(no_smp_data[[domains]], no_dom))

  if(length(no_dom) > 0){
    no_sub_no_area_data <- no_smp_data[vec, ]
    no_sub_in_area_data <- no_smp_data[-vec, ]
  } else{
    no_sub_in_area_data <- no_smp_data
  }

  x_var <- c(labels(terms(fixed)))
  X_no_sub_in_area <- data.frame(Intercept = rep(1, length(no_sub_in_area_data[[subdomains]])))
  X_no_sub_in_area <- cbind(X_no_sub_in_area, no_sub_in_area_data[, x_var])
  Dom_X_no_sub_in_area <- no_sub_in_area_data[[domains]]
  rownames(X_no_sub_in_area) <- no_sub_in_area_data[[subdomains]]

  v_tilde_no <- numeric()
  for (k in 1:length(no_sub_in_area_data[[subdomains]])) {
    v_tilde_no[k] <- v_tilde_i[which(rownames(v_tilde_i) == Dom_X_no_sub_in_area[k])]
  }

  mu_ij_no_sub_in_area <- as.matrix(X_no_sub_in_area) %*% beta_tild + v_tilde_no

  if(length(no_dom) > 0){
    X_no_sub_no_area <- data.frame(Intercept = rep(1, length(no_sub_no_area_data[[subdomains]])))
    X_no_sub_no_area <- cbind(X_no_sub_no_area, no_sub_no_area_data[, x_var])
    rownames(X_no_sub_no_area) <- no_sub_no_area_data[[subdomains]]

    mu_ij_no_sub_no_area <- as.matrix(X_no_sub_no_area) %*% beta_tild

    mu_ij_no_smp_dat <- data.frame(mu_ij = c(mu_ij_no_sub_no_area, mu_ij_no_sub_in_area),
                                   SubDomain = c(rownames(mu_ij_no_sub_no_area),
                                                 rownames(mu_ij_no_sub_in_area)))
  } else{
    mu_ij_no_smp_dat <- data.frame(mu_ij = c(mu_ij_no_sub_in_area),
                                   SubDomain = c(rownames(mu_ij_no_sub_in_area)))
  }


  mu_ij_smp_dat <- data.frame(mu_ij = mu_ij_smp, SubDomain = rownames(mu_ij_smp))

  mu_ij <- rbind(mu_ij_smp_dat, mu_ij_no_smp_dat)

  data <- merge(data, mu_ij, by.x = subdomains, by.y = "SubDomain")
  EBLUP_ij <- data.frame(EBLUP = data$mu_ij, SubArea = data[[subdomains]])
  ################################################################################
  # Preparation for back transformation
  gamma_ij <- data.frame(SubArea = names(gamma_ij),
                         gamma_ij = gamma_ij)
  gamma_i <- data.frame(Area = names(gamma2_i),
                        gamma_i = gamma2_i)
  gamma_i0 <- data.frame(Area = names(gamma1_i),
                         gamma_i0 = gamma1_i)
  y_bar_igamma <- data.frame(Area = names(y_bar_igamma),
                             y_bar_igamma = y_bar_igamma)
  X_bar_igamma_x_beta <- data.frame(Area = y_bar_igamma$Area,
                                    X_bar_igamma_x_beta = X_bar_igamma %*% beta_tild)

  data <- merge(data, gamma_ij, by.x = subdomains,
                by.y = "SubArea", all.x = TRUE, all.y = FALSE)
  data <- merge(data, gamma_i, by.x = domains,
                by.y = "Area", all.x = TRUE, all.y = FALSE)
  data <- merge(data, gamma_i0, by.x = domains,
                by.y = "Area", all.x = TRUE, all.y = FALSE)
  data <- merge(data,  y_bar_igamma, by.x = domains,
                by.y = "Area", all.x = TRUE, all.y = FALSE)
  data <- merge(data,  X_bar_igamma_x_beta, by.x = domains,
                by.y = "Area", all.x = TRUE, all.y = FALSE)

  data$gamma_ij <- ifelse(is.na(data$gamma_ij), 0, data$gamma_ij)
  data$gamma_i0 <- ifelse(is.na(data$gamma_i0) & data$gamma_ij == 0, 0, data$gamma_i0)
  data$gamma_i <- ifelse(is.na(data$gamma_i), 0, data$gamma_i)

  ################################################################################
  # Back transformation at sub-area level
  if(trafo == "arcsin"){
    #EBLUP_ij$EBLUP <- arcsin_backtransform_naive(eblup = EBLUP_ij$EBLUP)
    #data$mu_ij <- arcsin_backtransform_naive(eblup = data$mu_ij)

    bt_eblup <- arcsin_backtransform_bc(data = data, eblup = "mu_ij", subdomains = subdomains,
                                        varv = varv, varu = varu)
    bt_naive <- arcsin_backtransform_naive(data = data, eblup = "mu_ij", subdomains = subdomains)
    data <- merge(data, bt_eblup, by = subdomains, all.x = TRUE, all.y = FALSE)
    data <- merge(data, bt_naive[, c(subdomains, "bt_naive_mu_ij")], by = subdomains, all.x = TRUE, all.y = FALSE)
    EBLUP_ij <- data[, c("bt_mu_ij", "bt_naive_mu_ij", subdomains)]
    colnames(EBLUP_ij) <- c("EBLUP_BC", "EBLUP_Naive", "SubArea")
    #  if(all.equal(EBLUP_ij$SubArea, data[[subdomains]])){
    #    data$mu_ij <- EBLUP_ij$EBLUP
    #  }
  } else if(trafo == "log"){

    bt_eblup <- log_backtransform_bc(data = data, eblup = "mu_ij", subdomains = subdomains,
                                     varv = varv, varu = varu)
    bt_naive <- log_backtransform_naive(data = data, eblup = "mu_ij", subdomains = subdomains)
    data <- merge(data, bt_eblup, by = subdomains, all.x = TRUE, all.y = FALSE)
    data <- merge(data, bt_naive[, c(subdomains, "bt_naive_mu_ij")], by = subdomains, all.x = TRUE, all.y = FALSE)

    EBLUP_ij <- data[, c("bt_mu_ij", "bt_naive_mu_ij", subdomains)]
    colnames(EBLUP_ij) <- c("EBLUP_BC", "EBLUP_Naive", "SubArea")
  }
  ################################################################################
  # Obtain area level EBLUP by aggregating back-transformed at sub-area level

  if(trafo == "log" | trafo == "arcsin"){
    mu_i <- numeric()
    naive_mu_i <- numeric()
    eblup <- "bt_mu_ij"
    eblup_naive <- "bt_naive_mu_ij"
    for (i in 1:M) {
      dat <- data[data[[domains]] == DomName[i], c(eblup, eblup_naive, "ObsSub", nunits)]
      N_i <- sum(dat[[nunits]])
      mu_i[i] <- sum(dat[[eblup]] * dat[[nunits]])/N_i
      naive_mu_i[i] <- sum(dat[[eblup_naive]] * dat[[nunits]])/N_i
    }
    EBLUP_i <- data.frame(EBLUP_BC = mu_i, EBLUP_Naive = naive_mu_i, Area = DomName)
  }
  else {
    mu_i <- numeric()
    eblup <- "mu_ij"
    for (i in 1:M) {
      dat <- data[data[[domains]] == DomName[i], c(eblup, "ObsSub", nunits)]
      N_i <- sum(dat[[nunits]])
      mu_i[i] <- sum(dat[[eblup]] * dat[[nunits]])/N_i
    }
    EBLUP_i <- data.frame(EBLUP = mu_i, Area = DomName)
  }

  ##############################################################################
  # Get residuals and predicted random effects
  #all.equal(names(y_ij), rownames(X))
  residuals <- data.frame(SubArea = names(y_ij),
                          y_ij = y_ij,
                          xbeta = X %*% beta_tild)
  residuals <- merge(residuals, data[, c(domains, subdomains)],
                     by.x = "SubArea", by.y = subdomains, all.x = T, all.y = F)
  residuals <- merge(residuals, data.frame(Area = rownames(v_tilde_i),
                                           v_tilde = v_tilde_i),
                     by.x = domains, by.y = "Area", all.x = T, all.y = F)
  residuals <- merge(residuals, data.frame(SubArea = rownames(u_tilde_ij),
                                           u_tilde = u_tilde_ij),
                     by = "SubArea", all.x = T, all.y = F)
  mar_y_hat <- residuals$xbeta
  con_y_hat <- residuals$xbeta + residuals$v_tilde + residuals$u_tilde
  mar_r2 <- var(mar_y_hat)/var(y_ij)
  con_r2 <- var(con_y_hat)/var(y_ij)
  R2 <- c(mar_r2, con_r2)
  names(R2) <- c("Mar.R2", "Con.R2")
  resid <- residuals$y_ij - con_y_hat
  ##############################################################################

  return(list(EBLUP_Area = EBLUP_i,
              EBLUP_SubArea = EBLUP_ij,
              varu = varu,
              varv = varv,
              vare = vare,
              beta = beta_tild,
              v_tilde = v_tilde_i,
              u_tilde = u_tilde_ij,
              resid = resid,
              R2 = R2,
              data = data))
}
################################################################################
################################################################################
mse_fh_tf <- function(fixed, vardir, domains, subdomains, nunits, trafo,
                     eff_n, B, smp_hat_varu, smp_hat_varv, smp_hat_vare,
                     smp_beta, seed, maxit, data){
  # Data frame for bias
  bias_subarea <- data.frame(SubArea = unique(data[[subdomains]]))
  bias_area <- data.frame(Area = unique(data[[domains]]))
  if(trafo == "log" | trafo == "arcsin"){
    bias_subarea_naive <- data.frame(SubArea = unique(data[[subdomains]]))
    bias_area_naive <- data.frame(Area = unique(data[[domains]]))
  }

  #  for(b in 1:B) {
  b <- 1
  while(b <= B) {
    print(b)
    # Generate pseudo errors and pseudo true indicator (in transformed scale)
    current.na.action <- options('na.action')
    options(na.action = 'na.pass')
    X_pop <- model.matrix(object = fixed, data = data)
    rownames(X_pop) <- data[[subdomains]]
    options('na.action' = current.na.action$na.action)
    if(!is.na(eff_n)){
      new_data <- data[, c(subdomains, domains, vardir, nunits, eff_n, get.vars(fixed, data = data))]
    } else if(is.na(eff_n)){new_data <- data[, c(subdomains, domains, vardir, nunits, get.vars(fixed, data = data))]}


    xbeta <- X_pop %*% smp_beta
    rownames(xbeta)
    dat_xbeta <- data.frame(Subarea = rownames(xbeta), xbeta = xbeta)
    new_data <- merge(new_data, dat_xbeta, by.x = subdomains, by.y = "Subarea", all.x = TRUE,
                      all.y = FALSE)
    new_data$pseudo_u <- rnorm(n = length(unique(new_data[[subdomains]])),
                               mean = 0, sd = sqrt(smp_hat_varu))

    pseudo_v <- data.frame(Area = unique(data[[domains]]),
                           pseudo_v = rnorm(n = length(unique(new_data[[domains]])),
                                            mean = 0, sd = sqrt(smp_hat_varv)))
    new_data <- merge(new_data, pseudo_v, by.x = domains, by.y = "Area",
                      all.x = TRUE, all.y = FALSE)
    new_data <- merge(new_data, smp_hat_vare[, c(subdomains, "vare")], by = subdomains, all.x = TRUE,
                      all.y = FALSE)
    new_data$boot_true_ind_star <- new_data$xbeta + new_data$pseudo_v + new_data$pseudo_u
    for (i in 1:dim(new_data)[1]) {
      if(is.na(new_data$vare[i])){
        new_data$pseudo_e[i] <- NA
        new_data$boot_dir_ind_star[i] <- NA
      } else{
        new_data$pseudo_e[i] <- rnorm(1, mean = 0, sd = sqrt(new_data$vare[i]))
        new_data$boot_dir_ind_star[i] <- new_data$boot_true_ind_star[i] + new_data$pseudo_e[i]
        #new_data$boot_dir_ind <- new_data$xbeta + new_data$pseudo_v + new_data$pseudo_u
      }
    }
    # Get true indicator in original scale
    if(trafo == "log"){
      new_data$boot_true_ind <- exp(new_data$boot_true_ind_star)
      new_data$boot_dir_ind <- exp(new_data$boot_dir_ind_star)
      ##########################################################################
      # Adjust boot_var_dir-----------------------------------------------------
      new_data$boot_var_dir <- new_data$boot_dir_ind^2 * new_data$vare
      boot_vardir <- "boot_var_dir"
      ##########################################################################

    } else if(trafo == "arcsin"){
      # Why do we need to truncate? -> because of monotonity
      new_data$boot_true_ind_star <- ifelse(new_data$boot_true_ind_star < 0,
                                            0, new_data$boot_true_ind_star)

      new_data$boot_true_ind_star <- ifelse(new_data$boot_true_ind_star > (pi / 2),
                                            (pi / 2), new_data$boot_true_ind_star)

      new_data$boot_true_ind <- (sin(new_data$boot_true_ind_star))^2
      new_data$boot_dir_ind_star <- ifelse(new_data$boot_dir_ind_star < 0,
                                           0, new_data$boot_dir_ind_star)

      new_data$boot_dir_ind_star <- ifelse(new_data$boot_dir_ind_star > (pi / 2),
                                           (pi / 2), new_data$boot_dir_ind_star)

      new_data$boot_dir_ind <- (sin(new_data$boot_dir_ind_star))^2
    } else if(trafo == "no"){
      new_data$boot_true_ind <- new_data$boot_true_ind_star
      new_data$boot_dir_ind <- new_data$boot_dir_ind_star
    }
    # Get true indicator at area level
    DomName <- unique(new_data[[domains]])
    new_data_area <- data.frame(Area = DomName)
    for (m in 1:length(DomName)) {
      dat <- new_data[new_data[[domains]] == DomName[m], c("boot_true_ind", nunits)]
      N <- sum(dat[[nunits]])
      new_data_area$boot_true_ind[m] <- sum(dat[["boot_true_ind"]] * dat[[nunits]])/N
    }
    # Estimate boot TF indicator
    boot_fixed <- update(fixed, "boot_dir_ind ~ .")
    boot_point_tf <- NULL
    if (trafo == "log"){
      try(boot_point_tf <- point_fh_tf(fixed = boot_fixed, vardir = boot_vardir, domains = domains,
                                    subdomains = subdomains, nunits = nunits, trafo = trafo,
                                    maxit = maxit,
                                    eff_n = eff_n, data = new_data), silent = TRUE)
    } else if(trafo == "arcsin" | trafo == "no"){
      try(boot_point_tf <- point_fh_tf(fixed = boot_fixed, vardir = vardir, domains = domains,
                                    subdomains = subdomains, nunits = nunits, trafo = trafo,
                                    maxit = maxit,
                                    eff_n = eff_n, data = new_data), silent = TRUE)
    }

    if(!is.null(boot_point_tf)){
      point_subarea <- merge(new_data[, c(subdomains, "boot_true_ind")], boot_point_tf$EBLUP_SubArea,
                             by.x = subdomains, by.y = "SubArea", all.x = TRUE, all.y = FALSE)
      point_area <- merge(new_data_area, boot_point_tf$EBLUP_Area,
                          by = "Area", all.x = TRUE, all.y = FALSE)


      if(trafo == "log" | trafo == "arcsin"){
        point_subarea$bias <- point_subarea$EBLUP_BC - point_subarea$boot_true_ind
        point_area$bias <- point_area$EBLUP_BC - point_area$boot_true_ind

        bias_subarea <- merge(bias_subarea, point_subarea[, c(subdomains, "bias")],
                              by.x = "SubArea", by.y = subdomains, all.x = TRUE, all.y = FALSE)
        colnames(bias_subarea)[b + 1] <- paste0("b", b)

        bias_area <- merge(bias_area, point_area[, c("Area", "bias")],
                           by = "Area", all.x = TRUE, all.y = FALSE)
        colnames(bias_area)[b + 1] <- paste0("b", b)

        point_subarea$bias_naive <- point_subarea$EBLUP_Naive - point_subarea$boot_true_ind
        point_area$bias_naive <- point_area$EBLUP_Naive - point_area$boot_true_ind

        bias_subarea_naive <- merge(bias_subarea_naive, point_subarea[, c(subdomains, "bias_naive")],
                                    by.x = "SubArea", by.y = subdomains, all.x = TRUE, all.y = FALSE)
        colnames(bias_subarea_naive)[b + 1] <- paste0("b", b)

        bias_area_naive <- merge(bias_area_naive, point_area[, c("Area", "bias_naive")],
                                 by = "Area", all.x = TRUE, all.y = FALSE)
        colnames(bias_area_naive)[b + 1] <- paste0("b", b)
      } else if(trafo == "no"){
        point_subarea$bias <- point_subarea$EBLUP - point_subarea$boot_true_ind
        point_area$bias <- point_area$EBLUP - point_area$boot_true_ind
        bias_subarea <- merge(bias_subarea, point_subarea[, c(subdomains, "bias")],
                              by.x = "SubArea", by.y = subdomains, all.x = TRUE, all.y = FALSE)
        colnames(bias_subarea)[b + 1] <- paste0("b", b)

        bias_area <- merge(bias_area, point_area[, c("Area", "bias")],
                           by = "Area", all.x = TRUE, all.y = FALSE)
        colnames(bias_area)[b + 1] <- paste0("b", b)
      }
      b <- b + 1
    }
  }

  if(trafo == "no"){
    subarea_quality <- data.frame(SubArea = bias_subarea$SubArea,
                                  bias = rowMeans(bias_subarea[, 2:(B+1)]),
                                  MSE = rowMeans((bias_subarea[, 2:(B+1)])^2))
    area_quality <- data.frame(Area = bias_area$Area,
                               bias = rowMeans(bias_area[, 2:(B+1)]),
                               MSE = rowMeans((bias_area[, 2:(B+1)])^2))
  } else if(trafo == "log" | trafo == "arcsin"){
    subarea_quality <- data.frame(SubArea = bias_subarea$SubArea,
                                  bias = rowMeans(bias_subarea[, 2:(B+1)]),
                                  MSE = rowMeans((bias_subarea[, 2:(B+1)])^2),
                                  bias_naiveBT = rowMeans(bias_subarea_naive[, 2:(B+1)]),
                                  MSE_naiveBT = rowMeans((bias_subarea_naive[, 2:(B+1)])^2))
    area_quality <- data.frame(Area = bias_area$Area,
                               bias = rowMeans(bias_area[, 2:(B+1)]),
                               MSE = rowMeans((bias_area[, 2:(B+1)])^2),
                               bias_naiveBT = rowMeans(bias_area_naive[, 2:(B+1)]),
                               MSE_naiveBT = rowMeans((bias_area_naive[, 2:(B+1)])^2))
  }

  return(list(MSE_Area = area_quality,
              MSE_SubArea = subarea_quality))

}
################################################################################
ml_NR_diff <- function(varv, varu, vare, yranr, no, X, m){
  a <- 1
  b <- no[1]
  if (is.vector(X)) {X_i <- X[a:b]}
  else{X_i <- X[a:b,]}
  J <- matrix(1,ncol=no[1],nrow=no[1])
  V <- Var_y(varv = varv, varu = varu, vare = vare, no = no, m = m)
  V_i <- V[a:b, a:b]
  inv_Vi <- solve(V_i)
  beta_til <- beta_tilde(varv = varv, varu = varu, vare = vare,
                         yranr = yranr, no = no, m = m, X = X)$beta

  score_v_term1 <- sum(diag(inv_Vi %*% J))
  score_u_term1 <- sum(diag(inv_Vi))
  score_v_term2 <-t(yranr[a:b] - X_i %*% beta_til)%*%(-inv_Vi %*% J %*% inv_Vi)%*%(yranr[a:b] -X_i %*% beta_til)
  score_u_term2 <- t(yranr[a:b] - X_i %*% beta_til)%*%(-inv_Vi %*% inv_Vi)%*%(yranr[a:b] - X_i %*% beta_til)
  score_v <- - 0.5 * (score_v_term1 + score_v_term2)
  score_u <- - 0.5 * (score_u_term1 + score_u_term2)

  H_hatv <- 0.5 * sum(diag(inv_Vi %*% J %*% inv_Vi %*% J))
  H_hatu <- 0.5 * sum(diag(inv_Vi %*% inv_Vi))
  H_hatvu <- 0.5 * sum(diag(inv_Vi  %*% J %*% inv_Vi))
  for(i in 2:m){
    a <- a + no[i - 1]
    b <- b + no[i]
    if (is.vector(X)) {X_i <- X[a:b]}
    else{X_i <- X[a:b,]}
    J <- matrix(1, ncol = no[i], nrow = no[i])
    V_i <- V[a:b, a:b]
    inv_Vi <- solve(V_i)

    score_v_term1 <- sum(diag(inv_Vi %*% J))
    score_u_term1 <- sum(diag(inv_Vi))
    score_v_term2 <-t(yranr[a:b] - X_i %*% beta_til)%*%(-inv_Vi %*% J %*% inv_Vi)%*%(yranr[a:b] -X_i %*% beta_til)
    score_u_term2 <- t(yranr[a:b] - X_i %*% beta_til)%*%(-inv_Vi %*%inv_Vi)%*%(yranr[a:b] - X_i %*% beta_til)
    score_v <- score_v -0.5*(score_v_term1 + score_v_term2)
    score_u <- score_u -0.5*(score_u_term1 + score_u_term2)

    H_hatv <- H_hatv + 0.5*sum(diag(inv_Vi %*% J %*% inv_Vi %*% J))
    H_hatu <- H_hatu + 0.5*sum(diag(inv_Vi %*% inv_Vi))
    H_hatvu <- H_hatvu + 0.5*sum(diag(inv_Vi %*% J %*% inv_Vi))
  }
  # score elements (S(theta^(r)) for theta = (sigma_v, sigma_u)) in Eq. 6.17 (Morales et al., 2021)
  score <- c(score_v, score_u)

  # Fisher information matrix (F(\theta^(r))) = - (hessian matrix)^(-1)  (Eq. 6.21 Morales et al., 2021)
  Inv_F <- solve(matrix(c(H_hatv, H_hatvu,
                          H_hatvu, H_hatu),
                        nrow = 2, ncol = 2, byrow = T))
  NR_diff <- Inv_F %*% score

  return(NR_diff)
}

#-------------------------------------------------------------------------------
# Get variance of y: V(y)  with sub-area as unit (N x N matrix)
Var_y <- function(varv, varu, vare, no, m) {
  R <- diag(vare, ncol = sum(no))
  Z <- cbind(c(rep(1, no[1])), diag(1, ncol = no[1], nrow = no[1]))
  G <- diag(c(varv, rep(varu, no[1])))
  for (i in 2:m) {
    Znext <- cbind(c(rep(1,no[i])),diag(1,ncol=no[i], nrow=no[i]))
    Gnext <- diag(c(varv, rep(varu, no[i])))
    Z <- adiag(Z, Znext)
    G <- adiag(G, Gnext)
  }
  V <- R + Z %*% G %*% t(Z)
  return(V)
}

#
# Under/For Eq. 2.6: beta_tilde(delta) =
# (sum_i^m(t(X_i) %*% solve(V_i) %*% X_i))^-1 * sum(t(X_i) %*% solve(V_i) %*% y_i)
beta_tilde <- function(varv, varu, vare, yranr, no, m, X){
  V <- Var_y(varv = varv, varu = varu, vare = vare, no = no, m = m)
  inv_V <- solve(V)
  beta <- solve(t(X) %*% inv_V %*% X) %*% (t(X) %*% inv_V %*% yranr)
  return(list(beta = beta,
              V = V,
              inv_V = inv_V))
}

################################################################################
log_backtransform_naive <- function(data, eblup, subdomains){
  bt_naive <- data[, c(eblup, subdomains)]
  bt_naive$bt_naive_mu_ij <- exp(bt_naive[[eblup]])
  return(bt_naive)
}

arcsin_backtransform_naive <- function(data, eblup, subdomains) {

  bt_naive <- data[, c(eblup, subdomains)]

  # Truncate the estimates on the transformed scale
  bt_naive[[eblup]] <- ifelse(bt_naive[[eblup]] < 0,
                              0, bt_naive[[eblup]])

  bt_naive[[eblup]] <- ifelse(bt_naive[[eblup]] > (pi / 2),
                              (pi / 2), bt_naive[[eblup]])

  # Naively backtransform with the inverse of the transformation
  bt_naive$bt_naive_mu_ij <- (sin(bt_naive[[eblup]]))^2

  return(bt_naive)
}

log_backtransform_bc <- function(data, eblup, subdomains, varu, varv){
  #point_backtransformed <- exp(eblup)
  bt_in_smp <- data[data$ObsSub == "yes", ]
  # exp(eblup_transformed_scale + 0.5 * g1)
  g1_in <- (1 - bt_in_smp$gamma_ij)^2 * (varu + (1 - bt_in_smp$gamma_i) * varv) +
    bt_in_smp$gamma_ij^2 * bt_in_smp$vare
  bt_in_smp$bt_mu_ij <- exp(bt_in_smp[[eblup]] + 0.5 *  g1_in)

  bt_out_smp <- data[data$ObsSub == "no", ]
  g1_out <- varu * (1 + (bt_out_smp$gamma_i/bt_out_smp$gamma_i0) *
                      (1 - 2 * bt_out_smp$gamma_ij))
  bt_out_smp$bt_mu_ij <- exp(bt_out_smp[[eblup]] + 0.5 * g1_out)

  # naive back transformation for the out of sample grids in out of sample domains
  bt_out_smp$bt_mu_ij <- ifelse(bt_out_smp$gamma_ij == 0 &
                                  bt_out_smp$gamma_i0 == 0,
                                exp(bt_out_smp[[eblup]]), bt_out_smp$bt_mu_ij)

  bt_eblup <- rbind(bt_in_smp[, c("bt_mu_ij", subdomains)],
                    bt_out_smp[, c("bt_mu_ij", subdomains)])
  return(bt_eblup)
}

log_backtransform_bc_area <- function(data, eblup, domains, varu){
  g1 <- varu * (data$gamma_i/data$gamma_i0)
  data$bt_eblup <- exp(data[[eblup]] + 0.5 * g1)
  return(data)
}

arcsin_backtransform_bc <- function(data, eblup, subdomains, varu, varv) {
  bt_in_smp <- data[data$ObsSub == "yes", ]

  for (k in 1:length(bt_in_smp$mu_ij)) {
    g1_in <- (1 - bt_in_smp$gamma_ij)^2 * (varu + (1 - bt_in_smp$gamma_i) * varv) +
      bt_in_smp$gamma_ij^2 * bt_in_smp$vare
    bt_in_smp$bt_mu_ij[k] <- integrate(integ_function, lower = 0, upper = pi / 2,
                                       mean = bt_in_smp$mu_ij[k],
                                       sd = sqrt(g1_in[k]))$value
  }
  bt_out_smp <- data[data$ObsSub == "no", ]
  for (k in 1:length(bt_out_smp$mu_ij)) {
    g1_out <- varu * (1 + (bt_out_smp$gamma_i/bt_out_smp$gamma_i0) *
                        (1 - 2 * bt_out_smp$gamma_ij))
    #-------------------------------------------------------------------------------
    # Naive back transformation for out of sample grids in out of sample domains
    if(bt_out_smp$gamma_ij[k] == 0 & bt_out_smp$gamma_i0[k] == 0){
      if(bt_out_smp$mu_ij[k] < 0) {bt_out_smp$bt_mu_ij[k] <- (sin(0))^2}
      else if(bt_out_smp$mu_ij[k] > pi/2) {bt_out_smp$bt_mu_ij[k] <- (sin(pi/2))^2}
      else if(bt_out_smp$mu_ij[k] < pi/2 & bt_out_smp$mu_ij[k] > 0) {bt_out_smp$bt_mu_ij[k] <- (sin(bt_out_smp$mu_ij[k]))^2}
    } else {
      # bias corrected back transformation for out of sample grids but in sample domains
      bt_out_smp$bt_mu_ij[k] <- integrate(integ_function,
                                          lower = 0,
                                          upper = pi / 2,
                                          bt_out_smp$mu_ij[k],
                                          sqrt(g1_out[k]))$value
    }
    #-------------------------------------------------------------------------------
  }
  bt_eblup <- rbind(bt_in_smp[, c("bt_mu_ij", subdomains)],
                    bt_out_smp[, c("bt_mu_ij", subdomains)])
  return(bt_eblup)
}

arcsin_backtransform_bc_area <- function(data, eblup, domains, varu) {
  g1 <- varu * (data$gamma_i/data$gamma_i0)

  for (k in 1:length(data$gamma_i)) {
    data$bt_eblup[k] <- integrate(integ_function, lower = 0, upper = pi / 2,
                                  mean = data[[eblup]][k],
                                  sd = sqrt(g1[k]))$value
  }

  return(data)
}

integ_function <- function(x, mean, sd) {
  sin(x)^2 * dnorm(x, mean = mean, sd = sd)
}
