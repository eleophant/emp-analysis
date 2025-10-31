########################################
#
# EMP analysis - alpha diversity checks
#
# Eleonore Lebeuf-Taylor
#
########################################

################ PACKAGES ################
library(phyr)
library(ggplot2)
library(patchwork) # to arrange plots

################ RAW MODELS ################
mod_raw_obs_poisson <- pglmm(
  adiv_observed_otus ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "poisson", # because this is count data
  cov_ranef = list(host_scientific_name = host_phylo)
)
mod_raw_obs_gaussian <- pglmm(
  adiv_observed_otus ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

mod_raw_chao <- pglmm(
  adiv_chao1 ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian", # because this is continuous data
  cov_ranef = list(host_scientific_name = host_phylo)
)

mod_raw_shannon <- pglmm(
  adiv_shannon ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

mod_raw_faith <- pglmm(
  adiv_faith_pd ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

################ LOG MODELS ################

mod_log_obs_poisson <- pglmm(
  log(adiv_observed_otus) ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "poisson",
  cov_ranef = list(host_scientific_name = host_phylo)
)

mod_log_obs_gaussian <- pglmm(
  log(adiv_observed_otus) ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)
mod_log_chao <- pglmm(
  log(adiv_chao1) ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

mod_log_shannon <- pglmm(
  log(adiv_shannon) ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

mod_log_faith <- pglmm(
  log(adiv_faith_pd) ~ basic_sociality + basic_diet + (1|host_scientific_name),
  data = df_metadata_sub,
  family = "gaussian",
  cov_ranef = list(host_scientific_name = host_phylo)
)

################ GET RESIDS ################
# get residuals & fitted - raw
residuals_raw_obs_poisson <- residuals(mod_raw_obs_poisson)
fitted_raw_obs_poisson <- fitted(mod_raw_obs_poisson)

residuals_raw_obs_gaussian <- residuals(mod_raw_obs_gaussian)
fitted_raw_obs_gaussian <- fitted(mod_raw_obs_gaussian)

residuals_raw_chao <- residuals(mod_raw_chao)
fitted_raw_chao <- fitted(mod_raw_chao)

residuals_raw_shannon <- residuals(mod_raw_shannon)
fitted_raw_shannon <- fitted(mod_raw_shannon)

residuals_raw_faith <- residuals(mod_raw_faith)
fitted_raw_faith <- fitted(mod_raw_faith)

# get residuals & fitted - log
residuals_log_obs_poisson <- residuals(mod_log_obs_poisson)
fitted_log_obs_poisson <- fitted(mod_log_obs_poisson)

residuals_log_obs_gaussian <- residuals(mod_log_obs_gaussian)
fitted_log_obs_gaussian <- fitted(mod_log_obs_gaussian)

residuals_log_chao <- residuals(mod_log_chao)
fitted_log_chao <- fitted(mod_log_chao)

residuals_log_shannon <- residuals(mod_log_shannon)
fitted_log_shannon <- fitted(mod_log_shannon)

residuals_log_faith <- residuals(mod_log_faith)
fitted_log_faith <- fitted(mod_log_faith)

################ NORMALITY ################
#### shapiro-wilk ####
shapiro.test(residuals_raw_obs_poisson) # W = 0.72994, p-value < 2.2e-16
shapiro.test(residuals_raw_obs_gaussian) # W = 0.98358, p-value = 0.0002883
shapiro.test(residuals_raw_chao) # W = 0.98046, p-value = 5.795e-05
shapiro.test(residuals_raw_shannon) # W = 0.99164, p-value = 0.03292
shapiro.test(residuals_raw_faith) # W = 0.98763, p-value = 0.002795

shapiro.test(residuals_log_obs_poisson) # W = 0.91274, p-value = 6.574e-14
shapiro.test(residuals_log_obs_gaussian) # W = 0.96943, p-value = 4.404e-07
shapiro.test(residuals_log_chao) # W = W = 0.96604, p-value = 1.198e-07
shapiro.test(residuals_log_shannon) # W = 0.91668, p-value = 1.494e-13
shapiro.test(residuals_log_faith) # W = 0.96723, p-value = 1.871e-07

#### qq plots ####
# raw models
qq_raw_obs_poisson <- ggplot(data.frame(resid = residuals(mod_raw_obs_poisson)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Observed OTUs \n(Poisson)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_raw_obs_gaussian <- ggplot(data.frame(resid = residuals(mod_raw_obs_gaussian)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Observed OTUs \n(Gaussian)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_raw_chao <- ggplot(data.frame(resid = residuals(mod_raw_chao)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Chao1") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_raw_shannon <- ggplot(data.frame(resid = residuals(mod_raw_shannon)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_raw_faith <- ggplot(data.frame(resid = residuals(mod_raw_faith)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Faith PD") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

# log models
qq_log_obs_poisson <- ggplot(data.frame(resid = residuals(mod_log_obs_poisson)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Observed OTUs \n(Poisson)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_log_obs_gaussian <- ggplot(data.frame(resid = residuals(mod_log_obs_gaussian)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Observed OTUs \n(Gaussian)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_log_chao <- ggplot(data.frame(resid = residuals(mod_log_chao)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Chao1") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_log_shannon <- ggplot(data.frame(resid = residuals(mod_log_shannon)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Shannon") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

qq_log_faith <- ggplot(data.frame(resid = residuals(mod_log_faith)), aes(sample = resid)) +
  stat_qq() + stat_qq_line(color = "red") +
  ggtitle("Faith PD") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12))

#### combine plots ####
qq_raw_obs_poisson + qq_raw_obs_gaussian + qq_raw_chao +  qq_raw_shannon + qq_raw_faith +
  plot_layout(ncol = 5) +
  plot_annotation(title = "A. Raw")

qq_log_obs_poisson + qq_raw_obs_gaussian + qq_log_chao + qq_log_shannon + qq_log_faith +
  plot_layout(ncol = 5) +
  plot_annotation(title = "B. Log-transformed")


################ HOMOSCEDASTICITY ################

# residuals vs fitted

#### raw ####

plot(fitted_raw_obs_poisson, residuals_raw_obs_poisson,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nObs OTUs (Poisson)")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_raw_obs_poisson, residuals_raw_obs_poisson), col = "blue", lwd = 2) # looks terrible

plot(fitted_raw_obs_gaussian, residuals_raw_obs_gaussian,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nObs OTUs (Gaussian)")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_raw_obs_gaussian, residuals_raw_obs_gaussian), col = "blue", lwd = 2) # looks much better

plot(fitted_raw_chao, residuals_raw_chao,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nChao1")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_raw_chao, residuals_raw_chao), col = "blue", lwd = 2) # looks good

plot(fitted_raw_shannon, residuals_raw_shannon,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nShannon")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_raw_shannon, residuals_raw_shannon), col = "blue", lwd = 2) # looks good

plot(fitted_raw_faith, residuals_raw_faith,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nFaith PD")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_raw_faith, residuals_raw_faith), col = "blue", lwd = 2) # looks good

#### log ####

plot(fitted_log_obs_poisson, residuals_log_obs_poisson,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nLog Obs OTUs (Poisson)")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_log_obs_poisson, residuals_log_obs_poisson), col = "blue", lwd = 2) # looks not great

plot(fitted_log_obs_gaussian, residuals_log_obs_gaussian,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nLog Obs OTUs (Gaussian)")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_log_obs_gaussian, residuals_log_obs_gaussian), col = "blue", lwd = 2) # looks a bit better

plot(fitted_log_chao, residuals_log_chao,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nLog Chao1")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_log_chao, residuals_log_chao), col = "blue", lwd = 2)

plot(fitted_log_shannon, residuals_log_shannon,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nLog Shannon")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_log_shannon, residuals_log_shannon), col = "blue", lwd = 2)

plot(fitted_log_faith, residuals_log_faith,
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted \nLog Faith PD")
abline(h = 0, col = "red", lwd = 2)
lines(lowess(fitted_log_faith, residuals_log_faith), col = "blue", lwd = 2)
