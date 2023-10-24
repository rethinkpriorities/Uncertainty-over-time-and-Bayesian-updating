### code for Uncertainty over time and Bayesian updating analysis

library(pacman)
p_load(readr, tidyverse, ggplot2, metafor, clubSandwich, ggpubr)

dat <- read_rds("surrogacy_data.rds")

all <- dat %>%
  filter(comparison == "same_sample", # only same sample technique
         y_wave != 2, # only one set of possible surrogates so won't contribute to identification
         ) %>% 
  mutate(yt = paste(Study, y, treatment, sep =" "),
         error = error_s_lm_est_norm, # lm as it always produces estimates
         se = error_s_lm_se_norm,
         se_weight = 1/se) %>%
  group_by(Study) %>%
  mutate(study_weight = 1/n(),
         weight = se_weight*study_weight) %>%
  ungroup()

# meta-analysis of mean error for debiasing forecasts
r_hat <- 0.6 # default parameter
v_slm_all <- impute_covariance_matrix(all$se^2, cluster = all$Study, r = r_hat)

ma_all_error <- rma.mv(yi = error, V = v_slm_all, sparse = T, W = study_weight,
                     random = ~ 1 | Study / id, slab = Study, data = all)

ct_all_error <- coef_test(ma_all_error, vcov = "CR2", cluster=all$Study)

ct_all_error 
sqrt(sum(ma_all_error$sigma2))

# correct errors by subtracting mean bias
all <- mutate(all, adj_error = error - ct_all_error$beta)

# creating squared error / variance / noise
all <- mutate(all, sqerror = adj_error^2)

# geom smooth plot by study
study_smooth <- all %>%
  ggplot(aes(x = horizon, y = sqerror, color = Study)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = F, aes(weight = se_weight)) +
  facet_wrap(~ Study, ncol = 2, scales = "free_y") +
  ylab("Squared error") + xlab("") +
  scale_y_continuous(n.breaks = 4) +
  theme_minimal() + theme(legend.position = "none") + 
  theme(text = element_text(size = 16)) + 
  theme(axis.text = element_text(size = 12))

# geom smooth plot overall
overall_smooth <- all %>%
  ggplot(aes(x = horizon, y = sqerror)) +
  geom_point(alpha = 0.5, aes(color = Study, size = weight)) +
  geom_smooth(se = F, aes(weight = weight)) +
  ylab("Squared error") + xlab("Forecast horizon (years)") +
  coorevp_cartesian(ylim = c(0,0.3)) +
  theme_minimal() + theme(text = element_text(size = 16))

ggarrange(study_smooth, overall_smooth, ncol = 1, nrow = 2)

# mse meta-analysis
## no intercepts since MSE has to equal 0 when horizon is 0
## conditional on horizon, no FEs - linear
ma_all_mse_horizon_no <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon - 1, 
                                random = ~ 1 | Study / id, data = all,
                                W = study_weight, test = "t")

ct_all_mse_horizon_no <- coef_test(ma_all_mse_horizon_no, vcov = "CR2", cluster=all$Study)

## conditional on horizon and question FEs - linear
ma_all_mse_horizon_yt <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon + y*treatment - 1, 
                                random = ~ 1 | Study / id, data = all,
                                   W = study_weight, test = "t")

ct_all_mse_horizon_yt <- coef_test(ma_all_mse_horizon_yt, vcov = "CR2", cluster=all$Study)

horizon_no_beta <- ct_all_mse_horizon_no$beta[1] # coefficient on horizon
horizon_no_se <- ct_all_mse_horizon_no$SE[1] # se on horizon
horizon_yt_beta <- ct_all_mse_horizon_yt$beta[1] # coefficient on horizon
horizon_yt_se <- ct_all_mse_horizon_yt$SE[1] # se on horizon

horizon_no_beta
horizon_no_se
horizon_no_beta - 1.96*horizon_no_se # lower bound
horizon_no_beta + 1.96*horizon_no_se # upper bound
horizon_yt_beta
horizon_yt_se
horizon_yt_beta - 1.96*horizon_yt_se # lower bound
horizon_yt_beta + 1.96*horizon_yt_se # upper bound

# mse graphs - linear
mse_no <- all %>%
  ggplot(aes(x = horizon, y = sqerror)) +
  geom_point(alpha = 0.5, aes(color = Study, size = weight)) +
  geom_abline(color = "blue", slope = horizon_no_beta, intercept = 0) +
  geom_ribbon(aes(ymin = (horizon_no_beta - 1.96 * horizon_no_se) * horizon,
                  ymax = (horizon_no_beta + 1.96 * horizon_no_se) * horizon),
              fill = "grey", alpha = 0.5) +
  coord_cartesian(ylim = c(-0.05, 0.3)) +
  ylab("Squared error") + xlab("Forecast horizon (years)") +
  theme_minimal() + theme(text = element_text(size = 16)) 

mse_yt <- all %>%
  ggplot(aes(x = horizon, y = sqerror)) +
  geom_point(alpha = 0.5, aes(color = Study, size = weight)) +
  geom_abline(color = "blue", slope = horizon_yt_beta, intercept = 0) +
  geom_ribbon(aes(ymin = (horizon_yt_beta - 1.96 * horizon_yt_se) * horizon,
                  ymax = (horizon_yt_beta + 1.96 * horizon_yt_se) * horizon),
              fill = "grey", alpha = 0.5) +
  coord_cartesian(ylim = c(-0.05, 0.3)) +
  ylab("Squared error") + xlab("Forecast horizon (years)") +
  theme_minimal() + theme(text = element_text(size = 16))

ggarrange(mse_no, mse_yt, ncol = 2, nrow = 1, 
          labels = c("No fixed effects", "Question fixed effects"),
          common.legend = T, legend = "bottom")

# computing prior from distribution of RCT treatment effects
## RCT TE meta-analyis
v_rct_all <- impute_covariance_matrix(all$rct_x_se_norm^2, cluster = all$Study, r = r_hat)
ma_all_te <- rma.mv(yi = rct_x_beta_norm, V = v_rct_all, sparse = T,
                    random = ~ 1 | Study / id, data = all,
                    W = study_weight, test = "t")
sigma_all_te <- sum(ma_all_te$sigma2)

gc() # freeing memory

# Estimating Webb model
## prior
sigma2_u <- sigma_all_te 

## error - constant linear growth
t_max <- 1000000
t <- 1:t_max
sigma2_e_no_main <- horizon_no_beta*t
sigma2_e_no_wc <- (horizon_no_beta+1.96*horizon_no_se)*t
sigma2_e_yt_main <- horizon_yt_beta*t
sigma2_e_yt_wc <- (horizon_yt_beta+1.96*horizon_yt_se)*t

## error - rate of growth increasing/decreasing 1% each year
sigma2_yt_inc <- horizon_yt_beta*(t^(3/2))
sigma2_yt_dec <- horizon_yt_beta*(t^(2/3))

# discount rate
d_wide <- tibble(t = t, d_t_no = 1/(1+sigma2_e_no_main/sigma2_u), 
                     d_t_no_wc = 1/(1+sigma2_e_no_wc/sigma2_u),
                     d_t_yt = 1/(1+sigma2_e_yt_main/sigma2_u), 
                     d_t_yt_wc = 1/(1+sigma2_e_yt_wc/sigma2_u),
                     d_t_yt_inc = 1/(1+sigma2_yt_inc/sigma2_u),
                     d_t_yt_dec = 1/(1+sigma2_yt_dec/sigma2_u)
                     )

d_long <- d_wide %>%
  pivot_longer(-t, names_to = "Estimate", values_to = "D") %>%
  mutate(
    Estimate = case_when(
      Estimate == "d_t_no"     ~ "Central, no fixed effects",
      Estimate == "d_t_no_wc"  ~ "Upper-bound, no fixed effects",
      Estimate == "d_t_yt"     ~ "Central, fixed effects",
      Estimate == "d_t_yt_wc"  ~ "Upper-bound, fixed effects",
      Estimate == "d_t_yt_inc" ~ "Increasing rate of increase",
      Estimate == "d_t_yt_dec" ~ "Decreasing rate of increase",
      TRUE                     ~ Estimate
    )
  )

# graphs
# linear
## 1000: level-level
d_long %>%
  filter(t<=1000 & 
           Estimate != "Increasing rate of increase" &
           Estimate != "Decreasing rate of increase") %>%
  ggplot(aes(x=t, y = D, color = Estimate)) + geom_line() + 
    ylab("D(t)") + xlab("t (years)") + coord_cartesian(ylim=c(0,1)) +
  theme_minimal() + theme(text = element_text(size = 16)) 

gc() # freeing memory

## 1000000: log-log
d_long %>%
  filter(t<=1000000 & 
           Estimate != "Increasing rate of increase" & Estimate != "Decreasing rate of increase") %>%
  ggplot(aes(x=t, y = D, color = Estimate)) + geom_line() +
  ylab("D(t)") + xlab("t (years)") +
  scale_x_log10(breaks = 10^(0:6), 
                labels = c("1", "10", "100", "1,000", "10,000", "100,000", 
                           "1,000,000")) +
  scale_y_log10(breaks = 10^(0:-6), 
                labels = c("1", "0.1", "0.01", "0.001", "0.0001", "0.00001", 
                           "0.000001")) +
  theme_minimal() + theme(text = element_text(size = 16)) 

gc() # freeing memory

# non-linear
# sigma
tibble(t = t, `Constant rate of increase` = sigma2_e_yt_main, 
       `Increasing rate of increase`= sigma2_yt_inc, 
       `Decreasing rate of increase` = sigma2_yt_dec) %>%
  filter(t < 1000) %>%
  pivot_longer(-t, names_to = "Assumption", values_to = "Forecast noise") %>%
  ggplot(aes(x = t, y = `Forecast noise`, color = Assumption)) + geom_line() + 
  theme_minimal() + theme(text = element_text(size = 16)) 

## d
#d_long %>%
#  filter(t<=1000 & (Estimate == "Central, fixed effects" | 
#                    Estimate == "Increasing rate of increase" | 
#                    Estimate == "Decreasing rate of increase")) %>%
#  mutate(Assumption = if_else(Estimate == "Central, fixed effects", 
#                              "Constant rate of increase", Estimate)) %>%
#  ggplot(aes(x=t, y = D, color = Assumption)) + geom_line() +
#  ylab("D(t)") + xlab("t (years)") +
#  theme_minimal() + theme(text = element_text(size = 16)) 

d_long %>%
  filter(t<=1000000 & (Estimate == "Central, fixed effects" | 
                         Estimate == "Increasing rate of increase" | 
                         Estimate == "Decreasing rate of increase"))  %>%
  mutate(Assumption = if_else(Estimate == "Central, fixed effects", 
                              "Constant rate of increase", Estimate)) %>%
  ggplot(aes(x=t, y = D, color = Assumption)) + geom_line() +
  ylab("D(t)") + xlab("t (years)") +
  scale_x_log10(breaks = 10^(0:6), 
                labels = c("1", "10", "100", "1,000", "10,000", "100,000", 
                           "1,000,000")) +
  scale_y_log10(breaks = 10^(0:-6), 
                labels = c("1", "0.1", "0.01", "0.001", "0.0001", "0.00001", 
                           "0.000001")) +
  theme_minimal() + theme(text = element_text(size = 16)) 

## key values
### have to increase t_max to 100000000 to estimate some of these
index_10_yt <- which(d_wide$d_t_yt < 0.1)[1]
index_1_yt <- which(d_wide$d_t_yt < 0.01)[1]
index_0.1_yt <- which(d_wide$d_t_yt < 0.001)[1]
index_0.01_yt <- which(d_wide$d_t_yt < 0.0001)[1]
index_0.001_yt <- which(d_wide$d_t_yt < 0.00001)[1]

index_10_yt_wc <- which(d_wide$d_t_yt_wc < 0.1)[1]
index_1_yt_wc <- which(d_wide$d_t_yt_wc < 0.01)[1]
index_0.1_yt_wc <- which(d_wide$d_t_yt_wc < 0.001)[1]
index_0.01_yt_wc <- which(d_wide$d_t_yt_wc < 0.0001)[1]
index_0.001_yt_wc <- which(d_wide$d_t_yt_wc < 0.00001)[1]

index_10_no <- which(d_wide$d_t_no < 0.1)[1]
index_1_no <- which(d_wide$d_t_no < 0.01)[1]
index_0.1_no <- which(d_wide$d_t_no < 0.001)[1]
index_0.01_no <- which(d_wide$d_t_no < 0.0001)[1]
index_0.001_no <- which(d_wide$d_t_no < 0.00001)[1]

index_10_no_wc <- which(d_wide$d_t_no_wc < 0.1)[1]
index_1_no_wc <- which(d_wide$d_t_no_wc < 0.01)[1]
index_0.1_no_wc <- which(d_wide$d_t_no_wc < 0.001)[1]
index_0.01_no_wc <- which(d_wide$d_t_no_wc < 0.0001)[1]
index_0.001_no_wc <- which(d_wide$d_t_no_wc < 0.00001)[1]

index_10_inc <- which(d_wide$d_t_yt_inc < 0.1)[1]
index_1_inc <- which(d_wide$d_t_yt_inc < 0.01)[1]
index_0.1_inc <- which(d_wide$d_t_yt_inc < 0.001)[1]
index_0.01_inc <- which(d_wide$d_t_yt_inc < 0.0001)[1]
index_0.001_inc <- which(d_wide$d_t_yt_inc < 0.00001)[1]

index_10_dec <- which(d_wide$d_t_yt_dec < 0.1)[1]
index_1_dec <- which(d_wide$d_t_yt_dec < 0.01)[1]
index_0.1_dec <- which(d_wide$d_t_yt_dec < 0.001)[1]
index_0.01_dec <- which(d_wide$d_t_yt_dec < 0.0001)[1]
index_0.001_dec <- which(d_wide$d_t_yt_dec < 0.00001)[1]

gc()

# appendix
## MSE rather than noise

app <- all %>%
  mutate(sqerror = error^2) # no debiasing

# mse meta-analysis
## no intercepts since MSE has to equal 0 when horizon is 0
## conditional on horizon, no FEs - linear
ma_app_mse_horizon_no <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon - 1, 
                                random = ~ 1 | Study / id, data = app,
                                W = study_weight, test = "t")

ct_app_mse_horizon_no <- coef_test(ma_app_mse_horizon_no, vcov = "CR2", cluster=app$Study)

## conditional on horizon and question FEs - linear
ma_app_mse_horizon_yt <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon + y*treatment - 1, 
                                random = ~ 1 | Study / id, data = app,
                                W = study_weight, test = "t")

ct_app_mse_horizon_yt <- coef_test(ma_app_mse_horizon_yt, vcov = "CR2", cluster=app$Study)

horizon_yt_beta1 <- ct_app_mse_horizon_yt$beta[1] # coefficient on horizon
horizon_yt_se1 <- ct_app_mse_horizon_yt$SE[1] # se on horizon
horizon_no_beta1 <- ct_app_mse_horizon_no$beta[1] # coefficient on horizon
horizon_no_se1 <- ct_app_mse_horizon_no$SE[1] # se on horizon

horizon_yt_beta1
horizon_yt_se1
horizon_yt_beta1 - 1.96*horizon_yt_se1 # lower bound
horizon_yt_beta1 + 1.96*horizon_yt_se1 # upper bound
horizon_no_beta1
horizon_no_se1
horizon_no_beta1 - 1.96*horizon_no_se1 # lower bound
horizon_no_beta1 + 1.96*horizon_no_se1 # upper bound


## without study weights
# mse meta-analysis
## no intercepts since MSE has to equal 0 when horizon is 0
## conditional on horizon, no FEs - linear
ma_nsw_mse_horizon_no <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon - 1, 
                                random = ~ 1 | Study / id, data = all,
                                test = "t")

ct_nsw_mse_horizon_no <- coef_test(ma_nsw_mse_horizon_no, vcov = "CR2", cluster=all$Study)

## conditional on horizon and question FEs - linear
ma_nsw_mse_horizon_yt <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon + y*treatment - 1, 
                                random = ~ 1 | Study / id, data = all,
                                test = "t")

ct_nsw_mse_horizon_yt <- coef_test(ma_nsw_mse_horizon_yt, vcov = "CR2", cluster=all$Study)

horizon_yt_beta2 <- ct_nsw_mse_horizon_yt$beta[1] # coefficient on horizon
horizon_yt_se2 <- ct_nsw_mse_horizon_yt$SE[1] # se on horizon
horizon_no_beta2 <- ct_nsw_mse_horizon_no$beta[1] # coefficient on horizon
horizon_no_se2 <- ct_nsw_mse_horizon_no$SE[1] # se on horizon

horizon_yt_beta2
horizon_yt_se2
horizon_yt_beta2 - 1.96*horizon_yt_se2 # lower bound
horizon_yt_beta2 + 1.96*horizon_yt_se2 # upper bound
horizon_no_beta2
horizon_no_se2
horizon_no_beta2 - 1.96*horizon_no_se2 # lower bound
horizon_no_beta2 + 1.96*horizon_no_se2 # upper bound


## bias correction at question level
ma_per_q <- function(data){
  ma_result <- rma(yi = data$error, sei = data$se, data = data)
  return(ma_result$beta[1])
}

qlv <- group_by(all, yt) %>%
  mutate(q_bias = ma_per_q(cur_data_all()),
         adj_error = error - q_bias,
         sqerror = adj_error^2) %>%
  ungroup()

# mse meta-analysis
## no intercepts since MSE has to equal 0 when horizon is 0
## conditional on horizon, no FEs - linear
ma_qlv_mse_horizon_no <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon - 1, 
                                random = ~ 1 | Study / id, data = qlv,
                                W = study_weight, test = "t")

ct_qlv_mse_horizon_no <- coef_test(ma_qlv_mse_horizon_no, vcov = "CR2", cluster=qlv$Study)

## conditional on horizon and question FEs - linear
ma_qlv_mse_horizon_yt <- rma.mv(yi = sqerror, V = v_slm_all, sparse = T,
                                mods = ~ horizon + y*treatment - 1, 
                                random = ~ 1 | Study / id, data = qlv,
                                W = study_weight, test = "t")

ct_qlv_mse_horizon_yt <- coef_test(ma_qlv_mse_horizon_yt, vcov = "CR2", cluster=qlv$Study)

horizon_yt_beta3 <- ct_qlv_mse_horizon_yt$beta[1] # coefficient on horizon
horizon_yt_se3 <- ct_qlv_mse_horizon_yt$SE[1] # se on horizon
horizon_no_beta3 <- ct_qlv_mse_horizon_no$beta[1] # coefficient on horizon
horizon_no_se3 <- ct_qlv_mse_horizon_no$SE[1] # se on horizon

horizon_yt_beta3
horizon_yt_se3
horizon_yt_beta3 - 1.96*horizon_yt_se3 # lower bound
horizon_yt_beta3 + 1.96*horizon_yt_se3 # upper bound
horizon_no_beta3
horizon_no_se3
horizon_no_beta3 - 1.96*horizon_no_se3 # lower bound
horizon_no_beta3 + 1.96*horizon_no_se3 # upper bound

# With non-zero mean prior
ct_all_te <- coef_test(ma_all_te, vcov = "CR2", cluster=all$Study)
sqrt(sum(ma_all_te$sigma2))
mu_u <- ma_all_te$beta[1]

evp_wide <- d_wide %>%
  mutate(evp_t_no = mu_u*(1-D_t_no) + D_t_no,
         evp_t_no_wc = mu_u*(1-D_t_no_wc) + D_t_no_wc,
         evp_t_yt = mu_u*(1-D_t_yt) + D_t_yt,
         evp_t_yt_wc = mu_u*(1-D_t_yt_wc) + D_t_yt_wc,
         evp_t_yt_inc = mu_u*(1-D_t_yt_inc) + D_t_yt_inc,
         evp_t_yt_dec = mu_u*(1-D_t_yt_dec) + D_t_yt_dec) %>%
  select(-contains("D_"))

evp_long <- evp_wide %>%
  pivot_longer(-t, names_to = "Estimate", values_to = "EVP") %>%
  mutate(
    Estimate = case_when(
      Estimate == "evp_t_no"     ~ "Central, no fixed effects",
      Estimate == "evp_t_no_wc"  ~ "Upper-bound, no fixed effects",
      Estimate == "evp_t_yt"     ~ "Central, fixed effects",
      Estimate == "evp_t_yt_wc"  ~ "Upper-bound, fixed effects",
      Estimate == "evp_t_yt_inc" ~ "Increasing rate of increase",
      Estimate == "evp_t_yt_dec" ~ "Decreasing rate of increase",
      TRUE                     ~ Estimate
    )
  )

# graphs
## 1000000: log-log
evp_long %>%
  filter(t<=1000000) %>%
  ggplot(aes(x=t, y = EVP, color = Estimate)) + geom_line() +
  ylab("Posterior EV") + xlab("t (years)") +
  scale_x_log10(breaks = 10^(0:6), 
                labels = c("1", "10", "100", "1,000", "10,000", "100,000", 
                           "1,000,000")) +
  coord_cartesian(ylim = c(0,1)) +
  theme_minimal() + theme(text = element_text(size = 16)) 

gc()
