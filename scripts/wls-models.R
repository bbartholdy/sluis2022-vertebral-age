# Weighted least squares regression


# Load dependencies -------------------------------------------------------

library(here)
library(MASS)
library(equatiomatic)
#source(here("scripts", "script-1.R"))
#source(here("scripts", "lin-models.R"))

# Accuracy calculator function --------------------------------------------

acc_calc <- function(object, interval = "pred", level = 0.68){
  #if(class(object) == "earth") interval = "pint"
  pred <- do.call("predict", list(object, "interval" = interval, "level" = level))
  #pred <- predict.lm(object, interval, level)
  known <- model.frame(object)[,1] 
  fit <- pred[,1]
  lwr <- pred[,2]
  upr <- pred[,3]
  acc <- known >= lwr & known <= upr
  acc <- mean(acc) * 100
  # Accuracy (%) as the number of cases where the predicted value is within 5,
    # ...10, and 15 years of the known age
  acc_5 <- mean(abs(known - fit) < 5) * 100
  acc_10 <- mean(abs(known - fit) < 10) * 100
  acc_15 <- mean(abs(known - fit) < 15) * 100
  
  rss <- sum(object$residuals**2)
  press <- sum(rstandard(object, type = "predictive")**2)
  res_sd <- sqrt(rss / object$df.residual)
  rmse <- sqrt((1 / nrow(object$model) ) * rss)
  aic <- AIC(object)
  see <- sqrt(rss / object$df.residual)
  r.sq <- summary(object)$r.squared 
  out <- data.frame("R2" = r.sq, "SEE" = see, "RMSE" = rmse, "AIC" = aic,
                    "accuracy" = acc, "accuracy_5" = acc_5, 
                    "accuracy_10" = acc_10, "accuracy_15" = acc_15)
  
  # out <- list("accuracy" = acc, "accuracy_5" = acc_5, "accuracy_10" = acc_10, 
  #             "accuracy_15" = acc_15, "PRESS" = press, "AIC" = aic, "SEE" = see, 
  #             "SEE2" = see2)
  return(out)
}

# https://online.stat.psu.edu/stat501/lesson/13/13.1
# weight calculation for WLS regression
wt_calc <- function(data){
  #res <- residuals(model)
  res <- lm(known_age ~ region_mean, data)$residuals
  res <- do.call("abs", list(x = res))
  predictor <- model.frame(model)[,1]
  wt <- 1 / lm(abs(res) ~ predictor)$fitted.values^2
  return(wt)
}

wls_reg <- function(data){
  data <- na.omit(data)
  model <- lm(known_age ~ region_mean, data)
  abs_res <- abs(model$residuals)
  pred <- model.frame(model)[[1]]
  wt <- 1 / lm(abs_res ~ pred)$fitted.values^2
  wls <- lm(known_age ~ region_mean, data, weights = wt)
  return(wls)
}


sep_methods <- split(region_means_all, region_means_all$method) # separate by methods

sep_regions <- lapply(sep_methods, 
                      function(x) split(x, x$region)) # separate by regions

# Snodgrass models --------------------------------------------------------

# snodgrass_all <- sep_regions$snodgrass$all
# snodgrass_all_lm <- lm(known_age ~ region_mean, snodgrass_all)
# abs_res <- abs(snodgrass_all_lm$residuals)
# pred <- snodgrass_all$known_age
# snodgrass_all_wt <- 1 / lm(abs_res ~ pred)$fitted.values^2
# snodgrass_all_wls <- lm(known_age ~ region_mean, snodgrass_all, weights = wt)
# 
# snodgrass_cerv <- sep_regions$snodgrass$cervical
# res <- lm(known_age ~ region_mean, snodgrass_cerv)$residuals
# abs_res <- abs(res)
# pred <- na.omit(snodgrass_cerv)$known_age
# snodgrass_cerv_wt <- 1 / lm(abs_res ~ pred)$fitted.values^2
# snodgrass_cerv_wls <- lm(known_age ~ region_mean, snodgrass_cerv, weights = wt)
# 
# snodgrass_thor <- sep_regions$snodgrass$thoracic
# res <- lm(known_age ~ region_mean, snodgrass_thor)$residuals
# abs_res <- abs(res)
# pred <- snodgrass_thor$known_age
# snodgrass_thor_wt <- 1 / lm(abs_res ~ pred)$fitted.values^2
# snodgrass_thor_wls <- lm(known_age ~ region_mean, snodgrass_thor, weights = wt)
# 
# snodgrass_thor <- sep_regions$snodgrass$thoracic
# res <- lm(known_age ~ region_mean, snodgrass_thor)$residuals
# abs_res <- abs(res)
# pred <- snodgrass_thor$known_age
# snodgrass_thor_wt <- 1 / lm(abs_res ~ pred)$fitted.values^2
# snodgrass_thor_wls <- lm(known_age ~ region_mean, snodgrass_thor, weights = wt)
# 
# snodgrass_lumb <- sep_regions$snodgrass$lumbar
# res <- lm(known_age ~ region_mean, snodgrass_lumb)$residuals
# abs_res <- abs(res)
# pred <- snodgrass_lumb$known_age
# snodgrass_lumb_wt <- 1 / lm(abs_res ~ pred)$fitted.values^2
# snodgrass_lumb_wls <- lm(known_age ~ region_mean, snodgrass_lumb, weights = wt)


# snodgrass_wls <- list("all" = snodgrass_all_wls, "cerv" = snodgrass_cerv_wls, 
#                      "thor" = snodgrass_thor_wls, "lumb" = snodgrass_lumb_wls)

snodgrass_wls <- lapply(sep_regions$snodgrass, wls_reg)
watanabe_wls <- lapply(sep_regions$watanabe, wls_reg)
praneatpolgrang_wls <- lapply(sep_regions$praneatpolgrang, wls_reg)


# Extract equations -------------------------------------------------------

snodgrass_wls_eq <- lapply(snodgrass_wls, extract_eq, use_coefs = T)
watanabe_wls_eq <- lapply(watanabe_wls, extract_eq, use_coefs = T)
praneatpolgrang_wls_eq <- lapply(praneatpolgrang_wls, extract_eq, use_coefs = T)


# Model evaluation --------------------------------------------------------

snodgrass_eval <- lapply(snodgrass_wls, acc_calc)
watanabe_eval <- lapply(watanabe_wls, acc_calc)
praneatpolgrang_eval <- lapply(praneatpolgrang_wls, acc_calc)
all_eval_tbl <- data.frame(rbind(snodgrass_eval$all, watanabe_eval$all, 
                                 praneatpolgrang_eval$all))
all_eval_tbl$method <- c("snodgrass", "watanabe", "praneatpolgrang") 
cerv_eval_tbl <- data.frame(rbind(snodgrass_eval$cervical, watanabe_eval$cervical, 
                                 praneatpolgrang_eval$cervical))
cerv_eval_tbl$method <- c("snodgrass", "watanabe", "praneatpolgrang") 
thor_eval_tbl <- data.frame(rbind(snodgrass_eval$thoracic, watanabe_eval$thoracic, 
                                 praneatpolgrang_eval$thoracic))
thor_eval_tbl$method <- c("snodgrass", "watanabe", "praneatpolgrang") 
lumb_eval_tbl <- data.frame(rbind(snodgrass_eval$lumbar, watanabe_eval$lumbar, 
                                 praneatpolgrang_eval$lumbar))
lumb_eval_tbl$method <- c("snodgrass", "watanabe", "praneatpolgrang") 


acc_calc(snodgrass_wls$all)



# Iteratively re-weighted least squares -----------------------------------

snodgrass_iwls <- lapply(sep_regions$snodgrass, 
                            MASS:::rlm.formula, formula = known_age ~ region_mean)

watanabe_iwls <- lapply(sep_regions$watanabe, 
                           MASS:::rlm.formula, formula = known_age ~ region_mean)

praneatpolgrang_iwls <- lapply(sep_regions$praneatpolgrang, 
                                  MASS:::rlm.formula, formula = known_age ~ region_mean)
acc_calc(snodgrass_iwls$all,interval = "prediction", level = 0.68)
# acc_eval(snodgrass_models$all,interval = "prediction", level = 0.68)
# 
# pred_snodgrass_all <- predict(snodgrass_iwls$all, interval ="pred", level = 0.68)
# 
# pred_snodgrass_all
acc_calc(snodgrass_iwls$all)
# acc_calc(snodgrass_models$all)


# Age bias ----------------------------------------------------------------

snodgrass_all_res <- sep_regions$snodgrass$all %>%
  mutate(residual = snodgrass_wls$all$residuals)  # create column with residuals

snodgrass_cerv_res <- sep_regions$snodgrass$cervical %>%
  na.omit() %>%
  mutate(residual = snodgrass_wls$cervical$residuals)

snodgrass_thor_res <- sep_regions$snodgrass$thoracic %>%
  na.omit() %>%
  mutate(residual = snodgrass_wls$thoracic$residuals)

snodgrass_lumb_res <- sep_regions$snodgrass$lumbar %>%
  na.omit() %>%
  mutate(residual = snodgrass_wls$lumbar$residuals)

watanabe_all_res <- sep_regions$watanabe$all %>%
  mutate(residual = watanabe_wls$all$residuals) 

watanabe_cerv_res <- sep_regions$watanabe$cervical %>%
  na.omit() %>%
  mutate(residual = watanabe_wls$cervical$residuals)

watanabe_thor_res <- sep_regions$watanabe$thoracic %>%
  na.omit() %>%
  mutate(residual = watanabe_wls$thoracic$residuals)

watanabe_lumb_res <- sep_regions$watanabe$lumbar %>%
  na.omit() %>%
  mutate(residual = watanabe_wls$lumbar$residuals)

praneatpolgrang_all_res <- sep_regions$praneatpolgrang$all %>%
  mutate(residual = praneatpolgrang_wls$all$residuals) 

praneatpolgrang_cerv_res <- sep_regions$praneatpolgrang$cervical %>%
  na.omit() %>%
  mutate(residual = praneatpolgrang_wls$cervical$residuals)

praneatpolgrang_thor_res <- sep_regions$praneatpolgrang$thoracic %>%
  na.omit() %>%
  mutate(residual = praneatpolgrang_wls$thoracic$residuals)

praneatpolgrang_lumb_res <- sep_regions$praneatpolgrang$lumbar %>%
  na.omit() %>%
  mutate(residual = praneatpolgrang_wls$lumbar$residuals)

# Elements vs. residuals --------------------------------------------------

n_elements <- element_data_long %>% 
  group_by(id,method,region) %>%
  summarise(n = sum(!is.na(score))) %>%
  ungroup()

n_elements_all <- n_elements %>%
  group_by(id, method) %>%
  summarise(n = sum(n)) %>%
  mutate(region = "all") %>%
  ungroup()

n_all_snodgrass <- n_elements_all %>%
  filter(method == "snodgrass") %>%
  mutate(residual = snodgrass_all_res$residual,
         age = snodgrass_all_res$known_age,
         sex = snodgrass_all_res$sex)

n_all_snodgrass %>%
  ggplot(aes(x = n, y = abs(residual))) +
    geom_point(aes(col = sex)) +
    geom_smooth(method = "lm", se = F)
lm(residual ~ n + age, n_all_snodgrass)

n_elements %>%
  add_row(n_elements_all) %>%
  filter(region == "all")
