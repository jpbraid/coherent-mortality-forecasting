library(tidyverse)
library(demography)

# load in full mortality data; also convert qx's to mx's
mortality_data <- read_csv("S:/Agencies/ALT/ALT/ALT2015-17/jesse/4_mortality_improvement/output/qx_all_interpolated_M.csv") %>% 
  mutate(gender = "M") %>% 
  rbind(read_csv("S:/Agencies/ALT/ALT/ALT2015-17/jesse/4_mortality_improvement/output/qx_all_interpolated_F.csv") %>% 
  mutate(gender = "F")) %>% 
  gather(key = "year", value = "q_x", -c(age, gender)) %>%
  mutate(q_x_shift = lag(q_x, default = 0),
         m_x = ifelse(age == 0, q_x, q_x/(1 - (1/12)*(q_x_shift/(1 - q_x_shift)) - (5/12)*q_x))) %>% 
  select(-q_x_shift)

# this m_x calculation works because the age 0 m_x's are calculated first
# if things were calculated in a different order this totally wouldn't work

# we want a general "convert m_x to q_x" function as well
convert_to_qx <- function(mx) {
  qx <- mx
  for (i in 2:length(qx)) {
    qx[i] <- mx[i]*(1 - (1/12)*(qx[i - 1]/(1 - qx[i - 1])))/(1 + (5/12)*mx[i])
  }
  qx
}

# set up params
max_age <- 120
min_analysis_year <- 1961
max_analysis_year <- 2016

# select the relevant subset of the mortality data
mortality_data_subset <- mortality_data %>% filter(age %in% 0:max_age, year %in% min_analysis_year:max_year)

# the demogdata object wants matrices, not data frames, so...
mortality_matrix_M <- mortality_data_subset %>% 
  filter(gender == "M") %>%
  rename(Male = m_x, Year = year, Age = age) %>%
  select(c(Year, Age, Male)) %>% 
  spread(Year, Male) %>% 
  select(-Age) %>% 
  as.matrix()

mortality_matrix_F <- mortality_data_subset %>% 
  filter(gender == "F") %>%
  rename(Female = m_x, Year = year, Age = age) %>%
  select(c(Year, Age, Female)) %>% 
  spread(Year, Female) %>% 
  select(-Age) %>% 
  as.matrix()

# create the demogdata object from smoothed mortality rates
aus_mx <- read.demogdata("S:/Agencies/ALT/ALT/ALT2015-17/jesse/4_mortality_improvement/Mx_1x1.txt", 
                         "S:/Agencies/ALT/ALT/ALT2015-17/jesse/4_mortality_improvement/Exposures_1x1.txt", 
                         type = "mortality", label = sprintf("Australia_%s_%s", min_analysis_year, max_year))
names(aus_mx$rate) <- c("Female", "Male")
aus_mx$rate <- list(Female = mortality_matrix_F, Male = mortality_matrix_M)
aus_mx$pop <- list(Female = matrix(10000, nrow = max_age + 1, ncol = length(unique(mortality_data_subset$year))), 
                   Male = matrix(10000, nrow = max_age + 1, ncol = length(unique(mortality_data_subset$year))))
aus_mx$year <- min_analysis_year:max_year
aus_mx$age <- 0:max_age

# now fit fPCA / time series models to the smoothed mortality rates, and project forward
HBY_fit <- coherentfdm(aus_mx)
HBY_forecast <- forecast(object = HBY_fit, h = 2170 - max_year, level = 80, 
                         K = 100, drange = c(0, 0.5))





MI <- data.frame(gender = NA, start_year = NA, end_year = NA, age = NA, improvement_range = NA, MI = NA)

# project qx's based on different subsets of historical data (eg 1961 to 2001)
for (max_year in seq(from = max_analysis_year - 15, to = max_analysis_year, by = 5)) {
  # select the relevant subset of the mortality data
  mortality_data_subset <- mortality_data %>% filter(age %in% 0:max_age, year %in% min_analysis_year:max_year)
  
  # the demogdata object wants matrices, not data frames, so...
  mortality_matrix_M <- mortality_data_subset %>% 
    filter(gender == "M") %>%
    rename(Male = m_x, Year = year, Age = age) %>%
    select(c(Year, Age, Male)) %>% 
    spread(Year, Male) %>% 
    select(-Age) %>% 
    as.matrix()
  
  mortality_matrix_F <- mortality_data_subset %>% 
    filter(gender == "F") %>%
    rename(Female = m_x, Year = year, Age = age) %>%
    select(c(Year, Age, Female)) %>% 
    spread(Year, Female) %>% 
    select(-Age) %>% 
    as.matrix()
  
  # create the demogdata object from smoothed mortality rates
  aus_mx <- read.demogdata("S:/Agencies/ALT/ALT/ALT2015-17/jesse/4_mortality_improvement/Mx_1x1.txt", 
                           "S:/Agencies/ALT/ALT/ALT2015-17/jesse/4_mortality_improvement/Exposures_1x1.txt", 
                           type = "mortality", label = sprintf("Australia_%s_%s", min_analysis_year, max_year))
  names(aus_mx$rate) <- c("Female", "Male")
  aus_mx$rate <- list(Female = mortality_matrix_F, Male = mortality_matrix_M)
  aus_mx$pop <- list(Female = matrix(10000, nrow = max_age + 1, ncol = length(unique(mortality_data_subset$year))), 
                     Male = matrix(10000, nrow = max_age + 1, ncol = length(unique(mortality_data_subset$year))))
  aus_mx$year <- min_analysis_year:max_year
  aus_mx$age <- 0:max_age
  
  # now fit fPCA / time series models to the smoothed mortality rates, and project forward
  HBY_fit <- coherentfdm(aus_mx)
  HBY_forecast <- forecast(object = HBY_fit, h = 2170 - max_year, level = 80, 
                           K = 100, drange = c(0, 0.5))
  
  # for gender in c("M", "F") [i'm ok with doing two separate mortality data matrices]
    qx_forecast <- HBY_forecast$Male$rate$Male %>% apply(2, convert_to_qx) #%>% write.csv("male_forecasts.csv")

    # calculate year-on-year improvement factors
    indices <- 2:ncol(qx_forecast)
    qx_forecast_shifted <- cbind(qx_forecast[, indices], 0)
    yearly_improvement <- qx_forecast_shifted/qx_forecast - 1
  
    # calculate 25 and 125 year improvement factors
    # NOTE: this assumes columns are sorted by increasing year, gap is 1, etc.
    MI_25 <- (qx_forecast[, 25]/mortality_matrix_M[, ncol(mortality_matrix_M)])^(1/25) - 1
    MI_25_df <- as.data.frame(cbind(gender = "M", start_year = min_analysis_year, end_year = max_year, 
                               age = 0:max_age, improvement_range = 25, MI = MI_25))
    MI_125 <- (qx_forecast[, 125]/mortality_matrix_M[, ncol(mortality_matrix_M)])^(1/125) - 1
    MI_125_df <- as.data.frame(cbind(gender = "M", start_year = min_analysis_year, end_year = max_year, 
                                   age = 0:max_age, improvement_range = 125, MI = MI_125))
    
    MI <- rbind(MI, MI_25_df, MI_125_df)
}

MI <- MI %>% filter(complete.cases(MI))
rownames(MI) <- NULL


