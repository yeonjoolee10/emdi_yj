library(devtools)
load_all("C:/Users/lynye/Desktop/Diss/Code_Integration/emdi_yj")
#document("C:/Users/lynye/Desktop/Diss/Code_Integration/emdi_yj")
# Loading data - population and sample data
data("eusilcA_pop")
data("eusilcA_smp")

eusilcA_pop$dist_gen <- paste0(eusilcA_pop$district, "_",eusilcA_pop$gender)
eusilcA_smp$dist_gen <- paste0(eusilcA_smp$district, "_",eusilcA_smp$gender)

true_ind_sub <- emdi::direct(y = "eqIncome",
                         smp_data = eusilcA_pop,
                         smp_domains = "dist_gen",
                         var = F)

true_ind <- emdi::direct(y = "eqIncome",
                             smp_data = eusilcA_pop,
                             smp_domains = "district",
                             var = F)

fit_dir <- emdi::direct(y = "eqIncome",
                  smp_data = eusilcA_smp,
                  weights = "weight",
                  #threshold = 0.9 * median(eusilcA_smp$eqIncome),
                  smp_domains = "dist_gen",
                  var = T, na.rm = T, B = 50)

colnames(eusilcA_pop)
X_agg <- aggregate(eusilcA_pop[, c(3:15)], by = list(eusilcA_pop$dist_gen), FUN = mean)

X_agg <- merge(unique(eusilcA_pop[, c("district", "dist_gen")]), X_agg, by.x = "dist_gen",
               by.y = "Group.1")
comb_data <- merge(X_agg, fit_dir$ind[, c("Domain", "Mean", "Head_Count")],
                   by.x = "dist_gen", by.y = "Domain", all.x = T, all.y = F)
colnames(comb_data)[1] <- "dist_gen"

comb_data <- merge(comb_data, fit_dir$MSE[, c("Domain", "Mean", "Head_Count")],
                   by.x = "dist_gen", by.y = "Domain", all.x = T, all.y = F)
colnames(comb_data)[16:19] <- c("Dir_Mean", "Dir_HCR", "Var_Mean", "Var_HCR")

N_ik <- as.data.frame(table(eusilcA_pop$dist_gen))
colnames(N_ik) <- c("dist_gen", "N_ik")
comb_data <- merge(comb_data, N_ik, by = "dist_gen")
comb_data$smp <- ifelse(!is.na(comb_data$Dir_Mean), "in", "out")

library(sae)
eusilcA_smp$poor <- ifelse(eusilcA_smp$eqIncome <= 0.6 * median(eusilcA_smp$eqIncome), 1, 0)
dir_hcr_SRS <- sae::direct(y = poor,
                           dom = dist_gen,
                           domsize = comb_data[comb_data$smp == "in",
                                               c("dist_gen", "N_ik")],
                           data = eusilcA_smp)
dir_hcr_SRS <- merge(dir_hcr_SRS, comb_data[, c("dist_gen", "Var_HCR")],
                     by.x = "Domain", by.y = "dist_gen", all.x = T, all.y = F)
dir_hcr_SRS$DEff <- ifelse(dir_hcr_SRS$SD == dir_hcr_SRS$Var_HCR, 1, dir_hcr_SRS$SD^2/dir_hcr_SRS$Var_HCR)
dir_hcr_SRS$eff_n <- dir_hcr_SRS$SampSize/dir_hcr_SRS$DEff

comb_data <- merge(comb_data, dir_hcr_SRS[, c("Domain", "eff_n")],
                   by.x = "dist_gen", by.y = "Domain", all.x = T, all.y = F)

rm(list=setdiff(ls(), c("comb_data", "true_ind", "true_ind_sub")))
################################################################################

fit_fh_tf <- fh_tf(fixed = Dir_Mean ~ eqsize + cash + self_empl + unempl_ben,
                   vardir = "Var_Mean",
                   domains = "district",
                   subdomains = "dist_gen",
                   transformation = "no",
                   subdomsize = "N_ik",
                   MSE = T,
                   B = 10,
                   combined_data = comb_data)

fit_fh <- fh(fixed = Dir_Mean ~ eqsize + cash + self_empl + unempl_ben,
                   vardir = "Var_Mean",
                   domains = "district",
                   transformation = "log",
                   backtransformation = "bc_crude",
                   MSE = T,
                   B = 10,
                   combined_data = comb_data)

fit_fh_tf_log <- fh_tf(fixed = Dir_Mean ~ eqsize + cash + self_empl + unempl_ben,
                   vardir = "Var_Mean",
                   domains = "district",
                   subdomains = "dist_gen",
                   transformation = "log",
                   subdomsize = "N_ik",
                   MSE = T,
                   B = 10,
                   combined_data = comb_data)

fit_fh_tf_arc <- fh_tf(fixed = Dir_HCR ~ eqsize + cash + self_empl + unempl_ben,
                       vardir = "Var_HCR",
                       domains = "district",
                       subdomains = "dist_gen",
                       transformation = "arcsin",
                       subdomsize = "N_ik",
                       MSE = T,
                       B = 10,
                       combined_data = comb_data,
                       eff_smpsize = "eff_n")

summary(fit_fh)
summary(fit_fh_tf)
summary(fit_fh_tf_arc)
summary(fit_fh_tf_log)


################################################################################
res_mean <- true_ind$ind[, c("Domain", "Mean")]
res_mean <- merge(res_mean, fit_fh_tf$ind_Domain, by = "Domain")
res_mean <- merge(res_mean, fit_fh_tf_log$ind_Domain, by = "Domain")
res_mean
colnames(res_mean) <- c("Domain", "true", "TF", "TF_log")

summary((res_mean$TF - res_mean$true)/res_mean$true)
summary((res_mean$TF_log - res_mean$true)/res_mean$true)

res_mean_long <- data.frame(Domain = c(true_ind$ind$Domain, fit_fh_tf$ind_Domain$Domain,
                                      fit_fh_tf_log$ind_Domain$Domain),
                            Mean = c(true_ind$ind$Mean, fit_fh_tf$ind_Domain$EBLUP,
                                     fit_fh_tf_log$ind_Domain$EBLUP),
                            Model = rep(c("true", "TF", "TF_log"), by = length(true_ind$ind$Domain)))
library(ggplot2)

lp_mean_dom <- ggplot(data = res_mean_long) +
  geom_line(aes(x = Domain, y = Mean, group = Model, color = Model))
lp_mean_dom

res_mean_sub <- true_ind_sub$ind[, c("Domain", "Mean")]
res_mean_sub <- merge(res_mean_sub, fit_fh_tf$ind_Subdomain, by.x = "Domain",
                      by.y = "Subdomain")
res_mean_sub <- merge(res_mean_sub, fit_fh_tf_log$ind_Subdomain, by.x = "Domain",
                      by.y = "Subdomain")
res_mean_sub
colnames(res_mean_sub) <- c("Subdomain", "true", "TF", "TF_log")

summary((res_mean_sub$TF - res_mean_sub$true)/res_mean_sub$true)
summary((res_mean_sub$TF_log - res_mean_sub$true)/res_mean_sub$true)

