### reproduce the PM2.5 ERC analysis in HEI final reports 
### using the newly developed CausalGPS packages

library(devtools)
try(detach("package:CausalGPS", unload = TRUE), silent = TRUE)
install_github("fasrc/CausalGPS", ref="af7a58f2d8b8715b515b0554082d8781ff0323dd")

library("CausalGPS")
library("dplyr")
library("SuperLearner")
library("xgboost")

set.seed(1)

dir_data = '/nfs/home/X/xwu/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'
dir_out = '/nfs/home/X/xwu/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/'

# Import data from Priyanka's folder
load("/nfs/home/X/xwu/shared_space/ci3_analysis/pdez_measurementerror/National_Causal-master/balance_qd/covariates_qd.RData")
covariates_qd$year<-as.factor(covariates_qd$year)
covariates_qd$region<-as.factor(covariates_qd$region)
a.vals <- seq(min(covariates_qd$pm25_ensemble), max(covariates_qd$pm25_ensemble), length.out = 100)
delta_n <- (a.vals[2] - a.vals[1])

load(paste0(dir_data,"aggregate_data_qd.RData"))
aggregate_data_qd$year<-as.factor(aggregate_data_qd$year)
aggregate_data_qd$region<-as.factor(aggregate_data_qd$region)

dead_personyear<-aggregate(cbind(aggregate_data_qd$dead,
                                 aggregate_data_qd$time_count),
                           by=list( aggregate_data_qd$zip,
                                    aggregate_data_qd$year),
                           FUN=sum)
colnames(dead_personyear)[1:4]<-c("zip","year","dead","time_count")
dead_personyear[,"mortality"] <- dead_personyear[,"dead"]/dead_personyear[,"time_count"]

prematch_data <- merge(dead_personyear, covariates_qd,
                       by=c("zip", "year"))

treat = prematch_data["pm25_ensemble"]$pm25_ensemble
c = prematch_data[,c("mean_bmi", "smoke_rate", "hispanic", "pct_blk", "medhouseholdincome" ,"medianhousevalue",
                     "poverty", "education", "popdensity", "pct_owner_occ",
                     "summer_tmmx", "winter_tmmx", "summer_rmax", "winter_rmax",
                     "region", "year")]
Y = prematch_data[,"mortality"]

# quantile (0.01, 0.99) of exposures for triming
q1 <- quantile(prematch_data$pm25_ensemble, 0.01)
q2 <- quantile(prematch_data$pm25_ensemble, 0.99)
trimed_index <- which(a.vals >= q1 & a.vals <= q2)

# This implementation reproduces analyses in HEI reports)
# it is different than the default setting in the following way:
# 1. it uses a numeric matrix as the input of Xgboost model
# 2. it matches on the levels of a.vals
# 3. it fits a kernel smoother on untrimmed matched set and trims the ERC

# numeric matrix of covariates (matched the previous implementation that directly calls xgboost)
c2 = as.data.frame(data.matrix(c))

match_pop_all_noncompile_notrim <- generate_pseudo_pop(Y = Y,
                                                       w = treat,
                                                       c = c2,
                                                       ci_appr = "matching",
                                                       pred_model = "sl",
                                                       gps_model = "parametric",
                                                       use_cov_transform = FALSE,
                                                       sl_lib = c("m_xgboost"),
                                                       params = list("xgb_nrounds" = 50,
                                                                     "xgb_max_depth" = 6,
                                                                     "xgb_eta" = 0.3,
                                                                     "xgb_min_child_weight" = 1),
                                                       nthread = 8, 
                                                       covar_bl_method = "absolute",
                                                       covar_bl_trs = 0.5,
                                                       trim_quantiles = c(0.0,1.0), # no trimed
                                                       optimized_compile = FALSE, #created a column counter for how many times matched,
                                                       max_attempt = 1,
                                                       matching_fun = "matching_l1",
                                                       bin_seq = a.vals, # matched on a.vals
                                                       delta_n = delta_n,
                                                       scale = 1.0)
match_pop_data_notrim <- match_pop_all_noncompile_notrim$pseudo_pop

erf_notrim <- estimate_npmetric_erf(matched_Y = match_pop_data_notrim$Y,
                                    matched_w = match_pop_data_notrim$w,
                                    bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                    w_vals = a.vals,
                                    nthread = 10)
plot(erf_notrim)

# trimmed ERC
plot(erf_notrim$erf[trimed_index])
plot(a.vals[trimed_index], erf_notrim$erf[trimed_index]/erf_notrim$erf[trimed_index][1])

# After the development of CausalGPS, we recommend 
# the following implementation is advantegious
# Caveat: this recommendation only applies to a stable version of CausalGPS 
# If users prefer a numerical transformation of categorical variables in the GPS model
match_pop_all_noncompile <- generate_pseudo_pop(Y = Y,
                                                w = treat,
                                                c = c2,
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transform = FALSE,
                                                sl_lib = c("m_xgboost"),
                                                params = list("xgb_nrounds" = 50,
                                                              "xgb_max_depth" = 6,
                                                              "xgb_eta" = 0.3,
                                                              "xgb_min_child_weight" = 1),
                                                nthread = 8, # number of cores, you can change,
                                                covar_bl_method = "absolute",
                                                covar_bl_trs = 0.5,
                                                trim_quantiles = c(0.01,0.99), # trimed, you can change
                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                max_attempt = 1,
                                                matching_fun = "matching_l1",
                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                scale = 1.0)
match_pop_data_trim <- match_pop_all_noncompile$pseudo_pop

erf_trim <- estimate_npmetric_erf(matched_Y = match_pop_data_trim$Y,
                                  matched_w = match_pop_data_trim$w,
                                  bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                  w_vals = a.vals,
                                  nthread = 10)
plot(erf_trim)

plot(erf_trim$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim$erf[trimed_index]/erf_trim$erf[trimed_index][1])

# If users prefer a one-hot encoding of categorical variables in the GPS model
match_pop_all_noncompile_onehot <- generate_pseudo_pop(Y = Y,
                                                       w = treat,
                                                       c = c,
                                                ci_appr = "matching",
                                                pred_model = "sl",
                                                gps_model = "parametric",
                                                use_cov_transform = FALSE,
                                                sl_lib = c("m_xgboost"),
                                                params = list("xgb_nrounds" = 50,
                                                              "xgb_max_depth" = 6,
                                                              "xgb_eta" = 0.3,
                                                              "xgb_min_child_weight" = 1),
                                                nthread = 8, # number of cores, you can change,
                                                covar_bl_method = "absolute",
                                                covar_bl_trs = 0.5,
                                                trim_quantiles = c(0.01,0.99), # trimed, you can change
                                                optimized_compile = FALSE, #created a column counter for how many times matched,
                                                max_attempt = 1,
                                                matching_fun = "matching_l1",
                                                delta_n = delta_n, # you can change this to the one you used in previous analysis,
                                                scale = 1.0)
match_pop_data_trim_onehot <- match_pop_all_noncompile_onehot$pseudo_pop

erf_trim_onehot <- estimate_npmetric_erf(matched_Y = match_pop_data_trim_onehot$Y,
                                         matched_w = match_pop_data_trim_onehot$w,
                                         bw_seq = seq(8*delta_n, 40*delta_n, 2*delta_n),
                                         w_vals = a.vals,
                                         nthread = 10)
plot(erf_trim_onehot)

plot(erf_trim_onehot$erf[trimed_index])
plot(a.vals[trimed_index], erf_trim_onehot$erf[trimed_index]/erf_trim_onehot$erf[trimed_index][1])

