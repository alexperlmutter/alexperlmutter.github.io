## Suneet's analysis

# Load libraries
packages <- c("dplyr", "ggpubr", "geepack", "ggplot2", "car", "data.table", "reshape2")
lapply(packages, require, character.only = TRUE)

# Importing data
Data <- as.data.frame(read.csv("/Users/alexanderperlmutter/Google Drive (perlphd@gmail.com)/Suneet/Data.csv"))

# Creating binary strain variables at cutoffs -0.19 and -0.20 to 0.
Data <- Data %>% mutate(bin_strain = ifelse(LV.Global.Strain.. > -.20, 1, 0))

# Creating binary dilation variables at cutoffs 2 or greater
Data <- Data %>% mutate(LVED_Vol_Z_gt_2 = ifelse(LVED.Vol..Z.score. > 2, 1, 0))

# Creating visit number variable for each subject
Data <-Data %>%
  group_by(ID, group = cumsum(ID != lag(ID, default = first(ID)))) %>%
  mutate(Visit = row_number()) %>%
  ungroup() %>%
  dplyr::select(-group)

# Removing post-intervention observation for each subject
Data <- Data %>% 
  group_by(ID) %>% 
  slice(1:(n()-1))

# Organizing data by ID and Visit
Data$ID <- as.integer(gsub("S","",as.character(Data$ID)))
Data <- Data[order(Data$ID, Data$Visit),]

# Create time variable based dates of study
Data <- Data %>% 
  group_by(ID) %>%
  mutate(age_days = Age.at.time.of.study..days.) %>%
  dplyr::select(-Age.at.time.of.study..days.)

# Creating survival function - once patient has the outcome, they have it regardless of if it reverses
modify.vars <- function(x){
  x$LVED_Vol_Z_gt_2 <- ifelse(x$LVED_Vol_Z_gt_2 > 2, 1, x$LVED_Vol_Z_gt_2)
  x$bin_strain <- ifelse(x$bin_strain < -.20, 1, x$bin_strain)
  x
}
Data <- modify.vars(Data)

# Strain Z-scores
strain_Z <- Data$LV.Global.Strain.. - mean(na.omit(Data$LV.Global.Strain..))
Data$strain_Z <- strain_Z/sd(na.omit(Data$LV.Global.Strain..)) # positive z-scores will be strain
summary(Data$strain_Z)

# Strain centered at cutoff
Data$strain_center <- Data$LV.Global.Strain..+0.20

# Time to dilation and strain
Data <- Data %>%
  group_by(ID) %>%
  mutate(time_to_dilation = ifelse(LVED_Vol_Z_gt_2 == 1, age_days, 0)) %>%
  mutate(time_to_strain = ifelse(bin_strain == 1, age_days, 0))

# Create counter starting from 1 when positive for dilation or strain or both
Data <- na.omit(Data) %>%
  group_by(ID,LVED_Vol_Z_gt_2) %>%
  mutate(dilation_times = ifelse(LVED_Vol_Z_gt_2 == 1, row_number(), 0))
    
Data <- na.omit(Data) %>%
  group_by(ID,bin_strain) %>%
  mutate(strain_times = ifelse(bin_strain == 1, row_number(), 0))

Data <- na.omit(Data) %>%
  group_by(ID,LVED_Vol_Z_gt_2,bin_strain) %>%
  mutate(dil_strain_times = ifelse(LVED_Vol_Z_gt_2 == 1 & bin_strain == 1, row_number(), 0))

# What average number of follow-up days corresponds to each strain and dilation time? pretty much corresponds to the density plots.
age_times <- data.frame(matrix(NA, 2, 4))
age_times[,1:2] <- c(rbind(
  cbind(mean(na.omit(Data$age_days)[na.omit(Data$dilation_times == 1 & Data$LVED_Vol_Z_gt_2 == 1)]), last(addmargins(table(na.omit(Data$age_days)[na.omit(Data$dilation_times == 1 & Data$LVED_Vol_Z_gt_2 == 1)])))), 
  cbind(mean(na.omit(Data$age_days)[na.omit(Data$dilation_times == 2 & Data$LVED_Vol_Z_gt_2 == 1)]), last(addmargins(table(na.omit(Data$age_days)[na.omit(Data$dilation_times == 2 & Data$LVED_Vol_Z_gt_2 == 1)]))))))
age_times[,3:4] <- c(rbind(
  cbind(mean(na.omit(Data$age_days)[na.omit(Data$strain_times == 1 & Data$bin_strain == 1)]), last(addmargins(table(na.omit(Data$age_days)[na.omit(Data$strain_times == 1 & Data$bin_strain == 1)])))),
  cbind(mean(na.omit(Data$age_days)[na.omit(Data$strain_times == 2 & Data$bin_strain == 1)]), last(addmargins(table(na.omit(Data$age_days)[na.omit(Data$strain_times == 2 & Data$bin_strain == 1)]))))))
colnames(age_times) <- c("dilation age", "dilation N", "strain age", "strain N")
rownames(age_times) <- c("T1","T2")
age_times  # not sure how meaningful this is considering neonates were not systematically examined at the same age and N's are small, especially at T2. Regardless looks like strain is detectable a full week before dilation.

# Creating survival function - keep first time to event regardless of if patient maintains event
Data <- Data %>%
  group_by(ID) %>%
  mutate(time_to_dilation = ifelse(dilation_times == 0, 0, ifelse(dilation_times >= 1, time_to_dilation[dilation_times == 1], 0)),
         time_to_strain = ifelse(strain_times == 0, 0, ifelse(strain_times >= 1, time_to_strain[strain_times == 1], 0)))

# On average, how much earlier does strain occur than dilation at each time point?
Data <- Data %>%
  mutate(dilation_b4_strain = time_to_dilation - time_to_strain) # if negative, strain occurs before dilation

# Capture just the first time that strain or dilation is evident. If both occur at same Visit, use strain time.
Data <- Data %>%
  group_by(ID) %>%
  mutate(strain_diff_surv = ifelse(strain_times == 1 & dilation_times == 1, dilation_b4_strain[strain_times == 1],
                                      ifelse(strain_times == 1 & dilation_times == 0, dilation_b4_strain[strain_times == 1],
                                      ifelse(strain_times == 0 & dilation_times == 1, dilation_b4_strain[dilation_times == 1], 0))))

# Create event variable for when strain, dilation, or both occur
Data <- Data %>%
  group_by(ID) %>%
  mutate(event = ifelse(bin_strain == 1 & LVED_Vol_Z_gt_2 == 1, 3,
                        ifelse(bin_strain == 0 & LVED_Vol_Z_gt_2 == 1, 2,
                               ifelse(bin_strain == 1 & LVED_Vol_Z_gt_2 == 0, 1, 0))))

# Create first event variable for keeping only the event above that first occurs
Data <- Data %>%
  group_by(ID) %>%
  mutate(first_event = ifelse(strain_times == 1 & dilation_times == 1, first(event[strain_times == dilation_times]),
                              ifelse(strain_times == 1 & dilation_times == 0, event[strain_times == 1],
                                     ifelse(strain_times == 0 & dilation_times == 1, event[dilation_times == 1], 0))))
Data$first_event[Data$first_event == 0] <-"No event"
Data$first_event[Data$first_event == 1] <-"Strain"
Data$first_event[Data$first_event == 2] <-"Dilation"
Data$first_event[Data$first_event == 3] <-"Both"

# Create plot of age and time to each event by which event (dilation = positive, strain = negative) occurs first or both = 0
first_event_plot <- ggplot(na.omit(Data[Data$first_event != "No event",]), aes(x = age_days, y = time_to_strain)) +
  geom_line(aes(x = age_days, y = strain_diff_surv, group = factor(first_event), color = factor(first_event)), alpha = 0.5) +
  geom_point(aes(x = age_days, y = strain_diff_surv, group = factor(first_event), color = factor(first_event)), alpha = 0.5) +
  xlim(min(Data$age_days),max(Data$age_days)) +
  labs(x = "Age",  # this is population z-score
       y = "Time to strain (negative) or dilation (positive)",
       color = "First event") +
  theme_minimal() # Most of the dilation values are likely to remain evident on the second visit. If 0 is the population mean, a lot values are detectable at baseline dilation_times = 0.
first_event_plot

# based on event_plot create list of tables for each age day up to 40 that events occur
# more strain occurred from day 14 up to day 23
# more dilation occurred from day 25 on up to day 23
event_list <- vector("list", 65)
for(i in 1:length(event_list)){
  event_list[[i]] <- table(Data[Data$age_days<i,]$first_event)
  event_table <- plyr::ldply(event_list, rbind)[,-1]
}
event_table <- cbind(event_table,Days = as.integer(1:65))
for(j in 1:length(event_table)){
  event_table[,j] <- ifelse(is.na(event_table[,j]), 0, event_table[,j])
}
colnames(event_table) <- c("Strain", "Dilation", "Both", "Days")
library(reshape2)
event_table <- melt(event_table,id="Days")
event_plot <- ggplot(event_table,aes(x=Days,y=value,colour=variable,group=variable)) + 
  geom_line() +
  labs(x = "Days",  # this is population z-score
       y = "Number of events") +
  theme_minimal() # Most of the dilation values are likely to remain evident on the second visit. If 0 is the population mean, a lot values are detectable at baseline dilation_times = 0.
event_plot

# Extracting the first time a difference is evident for each patient
mean(Data$strain_diff_surv[Data$strain_diff_surv > 0]) # this is the average number of days that dilation is evident earlier than strain
mean(Data$strain_diff_surv[Data$strain_diff_surv < 0]) # this is the average number of days that strain is evident earlier than dilation

# To what extent is strain associated with dilation, accounting for within patient correlation and adjusting for height to weight ratio?
gee <- geeglm(LVED_Vol_Z_gt_2 ~ bin_strain, data=na.omit(Data), id=ID, waves = Visit, family=binomial, corstr="ex")
summary(gee) # No evidence of association

# Maybe your strain cutoffs are wrong because using Z-scores we see an association
gee.linear <- geeglm(LVED.Vol..Z.score. ~ strain_Z, data=na.omit(Data), id=ID, waves = Visit, family=gaussian, corstr="ex")
summary(gee.linear) # No evidence of association

# Plotting strain and Dilation Z-scores by visit number, so each participant who has a measure is plotted once on each plot
z_biv_plot_list <- vector("list", max(Data$Visit))
for(i in 1:max(Data$Visit)){
  z_biv_plot_list[[i]] <- 
    ggplot(na.omit(Data[Data$Visit == i,]), aes(x = strain_Z, y = LVED.Vol..Z.score.)) +
    geom_line(aes(y = LVED.Vol..Z.score.)) +
    geom_point(aes(y = LVED.Vol..Z.score.)) +
    geom_smooth(method = lm, se = FALSE) +
  labs(x = "Strain Z-score",  # this is sample z-score with mean 0 = -0.19 and variance 1
       y = "Dilation Z-score") # this is sample z-score based on population z-score
  }
z_biv_plot_list

