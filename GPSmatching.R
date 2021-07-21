# A guideline on the use of CausalGPS on Medicare data with individual-level strata

library(devtools)
#try(detach("package:CausalGPS", unload = TRUE), silent = TRUE)
#install_github("fasrc/CausalGPS", ref="develop", force = TRUE)
require(CausalGPS)
library("gnm")

load("/nfs/home/X/xwu/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_qd.RData")
covariates_qd$year <- as.factor(covariates_qd$year)
covariates_qd$region <- as.factor(covariates_qd$region)
a.vals <- seq(min(covariates_qd$pm25_ensemble), max(covariates_qd$pm25_ensemble), length.out = 50)
delta_n <- (a.vals[2] - a.vals[1])

match_pop_all <- generate_pseudo_pop(Y = covariates_qd$zip,
                                     w = covariates_qd$pm25_ensemble,
                                     c = covariates_qd[, c(4:19)],
                                     ci_appr = "matching",
                                     pred_model = "sl",
                                     gps_model = "parametric",
                                     use_cov_transform = TRUE,
                                     transformers = list("pow2", "pow3"),
                                     sl_lib = c("m_xgboost"),
                                     params = list(xgb_nrounds=c(50)),
                                     nthread = 8, # number of cores, you can change,
                                     covar_bl_method = "absolute",
                                     covar_bl_trs = 0.1,
                                     trim_quantiles = c(0.05,0.95), # trimed, you can change,
                                     optimized_compile = TRUE, #created a column counter for how many times matched,
                                     max_attempt = 5,
                                     matching_fun = "matching_l1",
                                     delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                     scale = 1.0)
match_pop_data <- match_pop_all$pseudo_pop
covariates_qd_trim <- subset(covariates_qd,
                             pm25_ensemble < quantile(covariates_qd$pm25_ensemble,0.95)&
                             pm25_ensemble > quantile(covariates_qd$pm25_ensemble,0.05))
match_pop_data <- cbind(match_pop_data, covariates_qd_trim[, 1:2])

dir_data = '/nfs/home/X/xwu/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/X/xwu/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
load(paste0(dir_data, "aggregate_data_qd.RData"))
aggregate_data_qd$year <- as.factor(aggregate_data_qd$year)
aggregate_data_qd$region <- as.factor(aggregate_data_qd$region)
aggregate_data_qd2 <- merge(aggregate_data_qd, match_pop_data[, c("year", "zip", "counter")], by = c("year", "zip"), all.y = TRUE)

matchingqd_gnm <- summary(gnm(dead ~ pm25_ensemble + offset(log(time_count)), 
                              eliminate = (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)),
                              data = aggregate_data_qd2,
                              family = poisson(link = "log"),
                              weights = counter))
exp(10*matchingqd_gnm$coefficients[1])
