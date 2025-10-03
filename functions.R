#-------------------------------------------------------------------------------
# @Organization - Climate Change Epidemiology Lab, Scripps Institution of Oceanography
# @Project - wildfire_transportability
# @Author - Michael Cheung, Candus Shi
# @Description - Functions for main analysis file
#-------------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Documentation: read_wildfire_data
# Usage: read_wildfire_data()
# @description: Read and format the 2018 California Wildfire data, 
#               merge the wildfire data with Census data

# Arguments/Options:

# @return: the data as described above
#-----------------------------------------------------------------------------

read_wildfire_data <- function() {
    
    ### from benmarhnia-lab/Wildfires_social_vulnerability/1_Wildfires_social_vulnerability_GSC.R
    
    # load and prepare data
    CA_hosp_County <- read.csv("https://raw.githubusercontent.com/benmarhnia-lab/Wildfires_social_vulnerability/refs/heads/main/CA_hosp_County_2018.csv")
    
    # define study period: from Sept 13 to Dec 5, 2018
    CA_hosp_County_week37_49 <- CA_hosp_County[which(CA_hosp_County$week>36 & CA_hosp_County$week<49),]
    
    CA_hosp_County_week45_49<-CA_hosp_County[which(CA_hosp_County$week>44 & CA_hosp_County$week<49),]
    
    # average smoke exposure during exposure period
    Counties_exp <- CA_hosp_County_week45_49 |>
        group_by(COUNTY_1) |>
        summarize(smoke_week45_49=mean(smoke))
    
    CA_hosp_County_week37_49 <- left_join(CA_hosp_County_week37_49, Counties_exp,
                                          by = c("COUNTY_1"))
    
    # define smoke exposure level
    CA_hosp_County_week37_49$Smokebin<-ifelse(CA_hosp_County_week37_49$smoke_week45_49>=10, 1,
                                              ifelse(CA_hosp_County_week37_49$smoke_week45_49<10, 0, NA))
    
    CA_hosp_County_week37_49$week<- as.numeric(CA_hosp_County_week37_49$week)
    
    ### census population DP5
    dp5 <- read.csv("ACSDP5Y2020.DP05-2025-02-19T205730.csv")
    dp5 <- dp5[2, colnames(dp5)[grepl("Estimate", colnames(dp5))]]
    colnames(dp5) <- str_extract(colnames(dp5), regex("(^.+)County"))
    colnames(dp5) <- str_replace(colnames(dp5), ".County", "")
    colnames(dp5) <- str_replace_all(colnames(dp5), "\\.", " ")
    dp5 <- data.frame(t(dp5))
    dp5 <- cbind(rownames(dp5), dp5)
    rownames(dp5) <- NULL
    colnames(dp5) <- c("countyname", "total.pop")
    dp5$total.pop <- as.numeric(gsub(",", "", dp5$total.pop))
    
    ### combine
    data <- left_join(CA_hosp_County_week37_49, dp5, by = c("countyname"))
    data$time_unit <- rep(1:84, 58)
    data$exp_day <- ifelse(data$Smokebin==1 & data$week >= 45, 1, 0)
    data$post_time <- ifelse(data$exp_day==1, data$time_unit - 56, 0)
    return(data)
}

#-----------------------------------------------------------------------------
# Documentation: estimate_effect
# Usage: estimate_effect(county_ind)
# @description: Estimate the effect of wildfire exposure using an interrupted
#               time series analysis. 

# Arguments/Options:
# @param county_ind: indices of the counties to be included in the analysis
# @param count: a boolean specifying whether the outcome is a count or continuous variable (DEFAULT: TRUE)

# @return: a vector of the effect estimate, 2.5% and 97.5% confidence interval
#-----------------------------------------------------------------------------

estimate_effect <- function(county_ind, count = TRUE) {
    
    # Interrupted time series regression analysis
    if (count == TRUE) {
        
        county_model <- glm(resp ~ offset(log(total.pop)) + 
                                exp_day + 
                                time_unit +
                                tmax + prec + hum + shrtwv_rad + wind,
                            family=quasipoisson,
                            data=exposed_data[county_ind, ])
    } else {
        
        county_model <- glm(resp_prop ~ exp_day + 
                                time_unit +
                                tmax + prec + hum + shrtwv_rad + wind,
                            family=gaussian,
                            data=exposed_data[county_ind, ])
    }
    
    estimate <- county_model$coefficients[2][1]
    estimate_ci <- confint(county_model)[2, ]
    
    if (count == TRUE) {
        estimate <- exp(estimate)
        estimate_ci <- exp(estimate_ci)
    }
    
    return(c(estimate, estimate_ci))
}

#-----------------------------------------------------------------------------
# Documentation: transport_iosw
# Usage: transport_iosw(county_1, county_0)
# @description: Estimate the transported effect using inverse odds of weighting.
#               Logistic regression and SuperLearner is used to estimate the 
#               probability of county inclusion.

# Arguments/Options:
# @param county_1: the county or set of counties from which the effect is transported
# @param county_0: the county or set of counties to which the effect is transported
# @param cps_SL_library: the set of prediction algorithms used in the SuperLearner

# @return: a tibble with the estimate type (logistic regression or SuperLearner),
#          the effect estimate, and the 95% confidence interval
#-----------------------------------------------------------------------------

transport_iosw <- function(county_1, county_0, 
                           zsite.SL.library = c('SL.glm', 'SL.ranger', 'SL.glmnet')) {
    
    # Define indicator for counties to transport to/from
    county_1_ind <- which(exposed_data$countyname %in% county_1)
    county_0_ind <- which(exposed_data$countyname %in% county_0)
    
    county_1_data <- exposed_data |>
        mutate(COUNTY_1 = case_when(countyname %in% county_1 ~ 1,
                                    countyname %in% county_0 ~ 0,
                                    T ~ NA)) |>
        filter(!is.na(COUNTY_1))
    
    # Estimate the probability of county inclusion using logistic regression
    cps_lr <- predict(glm(COUNTY_1 ~ tmax + prec + hum + shrtwv_rad + wind,
                          family = binomial,
                          data = county_1_data),
                      type = "response")
    
    # Estimate the probability of county inclusion using SuperLearner
    cps_sl <- SuperLearner(Y=county_1_data$COUNTY_1,
                           X=county_1_data[, c("tmax", "prec", "hum", "shrtwv_rad", "wind")],
                           SL.library=zsite.SL.library,
                           family='binomial')$SL.predict[,1]
    
    p_county_1 <- sum(county_1_data$COUNTY_1 == 1)/nrow(county_1_data)
    p_county_0 <- sum(county_1_data$COUNTY_1 == 0)/nrow(county_1_data)
    
    county_1_data <- county_1_data |>
        mutate(iosw_lr = ifelse(county_1_data$COUNTY_1 == 1, (1 - cps_lr)/cps_lr * p_county_1/p_county_0, 0),
               iosw_sl = ifelse(county_1_data$COUNTY_1 == 1, (1 - cps_sl)/cps_sl * p_county_1/p_county_0, 0))
    
    # county_1_data <- county_1_data |>
    #     group_by(COUNTY_1) |>
    #     mutate(iosw = pmin(pmax(iosw, quantile(iosw, 0.05)), quantile(iosw, 0.95)))
    
    outcome_mod_lr <- glm(resp ~ offset(log(total.pop)) + 
                              exp_day + 
                              time_unit + 
                              tmax + prec + hum + shrtwv_rad + wind, 
                          family = quasipoisson, 
                          data = county_1_data, 
                          weights = county_1_data$iosw_lr,
                          subset = county_1_data$iosw_lr != 0)
    
    outcome_mod_sl <- glm(resp ~ offset(log(total.pop)) + 
                              exp_day + 
                              time_unit + 
                              tmax + prec + hum + shrtwv_rad + wind, 
                          family = quasipoisson, 
                          data = county_1_data, 
                          weights = county_1_data$iosw_sl,
                          subset = county_1_data$iosw_sl != 0)
    
    estimate_lr <- exp(outcome_mod_lr$coefficients[2][1])
    ci_lr <- exp(confint(outcome_mod_lr)[2,])
    estimate_sl <- exp(outcome_mod_sl$coefficients[2][1])
    ci_sl <- exp(confint(outcome_mod_sl)[2,])
    
    transport_dat <- tibble(type = c("Transported (logistic regression)", "Transported (SuperLearner)"),
                            estimate = c(estimate_lr, estimate_sl),
                            lower_ci = c(ci_lr[1], ci_sl[1]),
                            upper_ci = c(ci_lr[2], ci_sl[2]))
    
    return(transport_dat)
}
