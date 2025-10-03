#-------------------------------------------------------------------------------
# @Organization - Climate Change Epidemiology Lab, Scripps Institution of Oceanography
# @Project - wildfire_transportability
# @Author - Michael Cheung, Candus Shi
# @Description - Main analysis file
#-------------------------------------------------------------------------------


# Preamble ------------------
library(tidyverse)
library(table1)
library(SuperLearner)
library(ggpubr)
library(grid)
# library(forecast)
# library(zoo)

# Set seed
set.seed(10149101)

# Source functions
source("functions.R")


# 1. Read data ------------------

# Read data
data <- read_wildfire_data()

# Create exposed county data set
exposed_data <- data |>
    filter(Smokebin == 1) |>
    mutate(resp_prop = resp/total.pop,
           date = ymd(paste0(year, "/", month, "/", day)))
    

# 2. Descriptive analysis ------------------

# Table 2: descriptive statistics
descriptive_data <- exposed_data |>
    mutate(tmax = structure(tmax, label = "Maximum temperature"),
           prec = structure(prec, label = "Precipitation"),
           hum = structure(hum, label = "Relative humidity"),
           shrtwv_rad = structure(shrtwv_rad, label = "Shortwave radiation"),
           wind = structure(wind, label = "Wind velocity"),
           exp_day = factor(exp_day, levels = c(0,1), labels = c("Before smoke exposure", "After smoke exposure")))

table1(~ tmax + prec + hum + shrtwv_rad + wind | exp_day,
       overall = c(left = "Total"),
       data = descriptive_data)

# Figure 1A: generate scatter plots by county with simple level change model prediction overlay (without time-varying confounders)
plot_data <- exposed_data |>
    group_by(countyname) |>
    mutate(resp_fitted = predict(glm(resp ~ offset(log(total.pop)) + exp_day + time_unit,
                                     family = quasipoisson),
                                 type = "response")) |>
    ungroup()

plot <- plot_data |>
    ggplot(aes(x = date, y = resp)) +
    geom_point(alpha = 5/10, size = 0.5) +
    geom_line(aes(x = date, y = resp_fitted), color = "red") +
    annotate('rect', 
             xmin=ymd("2018-11-08"), 
             xmax=ymd("2018-12-05"), 
             ymin=min(exposed_data$resp), 
             ymax=max(exposed_data$resp), 
             alpha=.2) +
    theme_light()

plot + 
    facet_wrap(vars(countyname)) + 
    labs(x = "Date (2018)", 
         y = "Total Respiratory Hospitalizations")
    

# 3. County-level interrupted time series analysis ------------------

# Estimate the effect of wildfire exposure for each county
exposed_county_effects <- t(sapply(unique(exposed_data$countyname), function(s) {
    county_ind <- which(exposed_data$countyname == s)
    s_estimate <- estimate_effect(county_ind, count = TRUE)
    return(s_estimate)
}))

# Figure 1: plot the effect estimates
exposed_county_effects <- tibble::rownames_to_column(data.frame(exposed_county_effects), "county")
colnames(exposed_county_effects) <- c("county", "estimate", "lower_ci", "upper_ci")
exposed_site_ates <- exposed_county_effects |>
    mutate(county = as.factor(county)) |>
    arrange(desc(estimate))

exposed_county_effects |>
    ggplot(aes(x = reorder(county, estimate),
               y = estimate)) +
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci)) + 
    geom_hline(yintercept = 1) + 
    coord_flip() +
    labs(x = "County", y = "Effect Estimate (RR)", linetype = "") + 
    theme_minimal()


# 4. Transport county-level effect estimates ------------------

observed <- exposed_county_effects |>
    mutate(type = "Observed")

# Estimate the transported effect for each county using all remaining counties
transport_all <- t(sapply(unique(exposed_data$countyname), function(s) {
    county_1 <- unique(exposed_data$countyname)[!(unique(exposed_data$countyname) == s)]
    effect <- transport_iosw(county_1 = county_1, county_0 = s)
    return(effect)
}))

transport_all <- as.data.frame(apply(transport_all, 2, unlist)) |>
    mutate(county = rep(rownames(transport_all), each = 2),
           type = case_when(type == "Transported (logistic regression)" ~ "Transported - all (logistic regression)",
                            type == "Transported (SuperLearner)" ~ "Transported - all (SuperLearner)"),
           across(c(estimate, lower_ci, upper_ci), as.numeric)) 

# Estimate the transported effect for each county using 50% of remaining counties (randomly selected)
transport_half <- t(sapply(unique(exposed_data$countyname), function(s) {
    other_counties <- unique(exposed_data$countyname)[!(unique(exposed_data$countyname) == s)]
    county_1 <- sample(other_counties, 0.5*length(other_counties), replace = F)
    effect <- transport_iosw(county_1 = county_1, county_0 = s)
    return(effect)
}))

transport_half <- as.data.frame(apply(transport_half, 2, unlist)) |>
    mutate(county = rep(rownames(transport_half), each = 2),
           type = case_when(type == "Transported (logistic regression)" ~ "Transported - 50% (logistic regression)",
                            type == "Transported (SuperLearner)" ~ "Transported - 50% (SuperLearner)"),
           across(c(estimate, lower_ci, upper_ci), as.numeric))

# Figure A2: plot the observed and logistic regression transported effect estimates (from all remaining counties and 50% remaning counties)
bind_rows(transport_all, transport_half, observed) |>
    filter(type == "Observed" 
           | type == "Transported - all (logistic regression)"
           | type == "Transported - 50% (logistic regression)") |>
    ggplot(aes(x = reorder(county, estimate),
               y = estimate,
               color = type,
               linetype = type)) +
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci),
                    position = position_dodge(0.7)) + 
    geom_hline(yintercept = 1) + 
    labs(x = "County", y = "Effect Estimate (RR)", linetype = "Estimate", color = "Estimate") + 
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    
top_5_counties <- c("Colusa", "Marin", "Mendocino", "San Mateo", "Merced")
sj_observed <- filter(exposed_county_effects, county == "San Joaquin") |>
    mutate(type = "Observed")

# Estimate the transported effects for each of Colusa, Marin, Mendocino, San Mateo, and Merced to San Joaquin
transport_pairwise <- t(sapply(top_5_counties, function(s) transport_iosw(county_1 = s, county_0 = "San Joaquin")))

transport_pairwise <- as.data.frame(apply(transport_pairwise, 2, unlist)) |>
    mutate(county = paste0(rep(rownames(transport_pairwise), each = 2), " to ", "San Joaquin"),
           across(c(estimate, lower_ci, upper_ci), as.numeric)) |>
    bind_rows(sj_observed)
    
fig_2A <- transport_pairwise %>%
    ggplot(aes(x = factor(county, level = unique(county)),
               y = estimate,
               color = type,
               linetype = type)) + 
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci),
                    position = position_dodge(0.25)) +
    geom_hline(yintercept = 1) +
    labs(x = "", y = "Effect Estimate (RR)", linetype = "") + 
    theme_minimal() +
    labs(color = "Estimate", linetype = "Estimate")

# Estimate the transported effect for Colusa, Marin, Mendocino, San Mateo, and Merced to San Joaquin
transport_jointly <- transport_iosw(top_5_counties, "San Joaquin")

transport_jointly <- sj_observed |>
    select(-county) |>
    bind_rows(transport_jointly)
    
fig_2B <- transport_jointly %>%
    ggplot(aes(x = type,
               y = estimate,
               color = type, linetype = type)) + 
    geom_pointrange(aes(ymin = lower_ci, ymax = upper_ci),
                    position = position_dodge(0.25)) +
    geom_hline(yintercept = 1) +
    labs(x = "San Joaquin", y = "Effect Estimate (RR)") + 
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    labs(color = "Estimate", linetype = "Estimate")

# Figure 2: plot the pairwise transported effects and jointly transported effect from Colusa, Marin, Mendocino, San Mateo, and Merced to San Joaquin
fig_2 <- ggarrange(fig_2A, fig_2B, 
                   nrow = 2, 
                   labels = c("A", "B"),
                   common.legend = TRUE, 
                   legend = "right")
