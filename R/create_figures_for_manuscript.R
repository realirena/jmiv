rm(list=ls())
library(ggplot2)
library(lubridate)
library(dplyr)
library(scales) 
library(data.table)
# 
fat_change_data <- read.csv("~/research_project_data/SWAN data/prepped_model_data/fat_data_cov_03252022.csv")

lean_change_data <- read.csv("~/research_project_data/SWAN data/prepped_model_data/lean_data_cov_04202022.csv")


## before detrending 
hormone_data <- read.csv("/Users/irena/Documents/repos/swan_models/R/data/SWAN Longitudinal Hormone Data 100220.csv")

### after detrending
#hormone_data <- read.csv("/Users/irena/Documents/repos/swan_models/R/data/hormone_data_fat_03252022.csv")

### plot the hormone trajectories 

g1 <- ggplot(hormone_data, 
             aes(x=time_fmp, y=log(ESTRADIOL))) + 
  geom_point(color="#8d8df2")+ 
  geom_smooth(method = "loess", se = TRUE, size=1.25 , color="#d66418") +
  labs(x="", y = "E2 Measurement",  title="E2 Measurements")+
  scale_x_continuous(breaks = pretty_breaks(n=10)) +
  theme(plot.title=element_text(size=30, face="bold", hjust = 0.5),
        plot.subtitle=element_text(size=25, face="bold", hjust = 0.5),
        legend.position = "none",
        axis.title =  element_text(size=20),
        axis.text.x = element_text(size=20,angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=20)) 

g2 <- ggplot(hormone_data, 
             aes(x=time_fmp, y=log(FSH))) + 
  geom_point(color="#49bab6")+ 
  geom_smooth(method = "loess", se = TRUE, size=1.25 , color="#d66418") +
  scale_x_continuous(breaks = pretty_breaks(n=10)) +
  labs(x="Time to FMP (Final Menstrual Period)", y = "FSH Measurement",  title="FSH Measurements")+
  theme(plot.title=element_text(size=30, face="bold", hjust = 0.5),
        plot.subtitle=element_text(size=25, face="bold", hjust = 0.5),
        legend.position = "none",
        axis.title =  element_text(size=20),
        axis.text.x = element_text(size=20,angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=20)) 

## save as 1800 x 900
gridExtra::grid.arrange(g1, g2)


## compute summary statistics for hormones and bm outcomes 

mean(fat_hormone$e2_loess_resid)
sd(fat_hormone$e2_loess_resid)
mean(fat_hormone$fsh_loess_resid)
sd(fat_hormone$fsh_loess_resid)


(table(fat_change_data$re)/841) * 100

table(fat_change_data$SPORTSGROUP)
(table(fat_change_data$SPORTSGROUP)/841) * 100

mean(fat_change_data$fat_change_rate)
mean(lean_change_data$lean_change_rate)

mean(fat_change_data$first_fat_perc)
sd(fat_change_data$first_fat_perc)

mean(lean_change_data$first_lean_perc)
sd(lean_change_data$first_lean_perc)


## create histograms of the rates of change 
g1 <- ggplot(fat_change_data, aes(x=fat_change_rate)) +  
  geom_histogram(fill="#4a95b2", colour="#56b1d4")+ 
  geom_vline(aes(xintercept=median(na.omit(fat_change_rate))),
             color="#e27abd", linetype="dashed", size=1.5) + 
  labs(title="Fat Mass Rate Change",
       subtitle = "Median Change: 0.001 rate increase",
       x="Change from First to Last Visit") + 
  theme(plot.title=element_text(size=30, face="bold", hjust=0.5),
        plot.subtitle=element_text(size=26, face="bold", hjust=0.5),
        axis.title= element_text(size=26),
        axis.text.x = element_text(size=20,angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=20))


g2 <- ggplot(lean_change_data, aes(x=lean_change_rate)) +  
  geom_histogram(fill="#378582", colour="#BFE4DE")+ 
  geom_vline(aes(xintercept=median(na.omit(lean_change_rate))),
             color="#e27abd", linetype="dashed", size=1.5) + 
  labs(title="Lean Mass Rate Change",
       subtitle = "Median Change: 0.002 rate decrease",
       x="Change from First to Last Visit", y="") + 
  theme(plot.title=element_text(size=30, face="bold", hjust=0.5),
        plot.subtitle=element_text(size=26, face="bold", hjust=0.5),
        axis.title= element_text(size=26),
        axis.text.x = element_text(size=20,angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=20))

## save as 1600 x 1200
gridExtra::grid.arrange(g1, g2, ncol=2)

### create ppd distributions 
e2_ppd_values <- read.csv("/Users/irena/repos/swan_models/R/figures/revisions_1013/fat_model_e2_ppd_values_0921.csv")

e2_ppd_values$obs_less_ind <- ifelse(e2_ppd_values$observed_value<e2_ppd_values$sim_value, 1,0)

e2_p_values <- e2_ppd_values %>%
  group_by(id) %>%
  summarise(sum_inds = sum(obs_less_ind))

e2_p_values$ppd_p_value <- round(e2_p_values$sum_inds/1000,2 )
#e2_p_values$e2_ppd_quantiles <- cut(e2_p_values$ppd_p_value,breaks =seq(0,1,0.05))

g1 <- ggplot(e2_p_values, aes(x=ppd_p_value)) +
  geom_histogram(bins=50,fill="#BFBFEC", color="#8d8df2") +
  labs(title="E2 posterior predictive p-values (Fat Mass Model)",
       subtitle = expression(paste("Based on 1,000 draws of model posterior samples")), x="") + 
  theme(plot.title = element_text(hjust = 0.5,size=30, face="bold"),
        plot.subtitle = element_text(hjust = 0.5,size=30, face="bold"),
        axis.title = element_text(size=25, face="bold"),
        axis.text.x = element_text(size=20,vjust = 1, hjust = 1),
        axis.text = element_text(size=20, face="bold"))

fsh_ppd_values <- read.csv("/Users/irena/repos/swan_models/R/figures/revisions_1013/fat_model_fsh_ppd_values_0921.csv")

fsh_p_values <- fsh_ppd_values %>%
  group_by(id) %>%
  summarise(sum_inds = sum(obs_less_ind))
fsh_p_values$ppd_p_value <- round(fsh_p_values$sum_inds/1000, 2)

g2 <- ggplot(fsh_p_values, aes(x=ppd_p_value)) +
  geom_histogram(bins=50,fill="#A2D7D5", color="#49bab6") +
  labs(title="FSH posterior predictive p-values (Fat Mass Model)",
       subtitle = expression(paste("Based on 1,000 draws of model posterior samples")), x="", y="") + 
  theme(plot.title = element_text(hjust = 0.5,size=30, face="bold"),
        plot.subtitle = element_text(hjust = 0.5,size=30, face="bold"),
        axis.title = element_text(size=25, face="bold"),
        axis.text.x = element_text(size=20,vjust = 1, hjust = 1),
        axis.text = element_text(size=20, face="bold"))


fat_mass_ppd <- grid.arrange(g1, g2, ncol=2)

### waist circumference ppd 

### create ppd distributions 
e2_waist <- read.csv("/Users/irena/repos/swan_models/R/figures/revisions_1013/waist_e2_ppd_values_1013.csv")

waist_e2_p_values <- e2_waist %>%
  group_by(id) %>%
  summarise(sum_inds = sum(obs_less_ind))

waist_e2_p_values$ppd_p_value <- round(waist_e2_p_values$sum_inds/1000,2 )
#e2_p_values$e2_ppd_quantiles <- cut(e2_p_values$ppd_p_value,breaks =seq(0,1,0.05))

g3 <- ggplot(waist_e2_p_values, aes(x=ppd_p_value)) +
  geom_histogram(bins=50,fill="#BFBFEC", color="#8d8df2") +
  labs(title="E2 posterior predictive p-values (Waist Circumference Model)",
       subtitle = expression(paste("Based on 1,000 draws of model posterior samples")), x="p-value") + 
  theme(plot.title = element_text(hjust = 0.5,size=30, face="bold"),
        plot.subtitle = element_text(hjust = 0.5,size=30, face="bold"),
        axis.title = element_text(size=25, face="bold"),
        axis.text.x = element_text(size=20,vjust = 1, hjust = 1),
        axis.text = element_text(size=20, face="bold"))

fsh_waist <- read.csv("/Users/irena/repos/swan_models/R/figures/revisions_1013/waist_fsh_ppd_values_1013.csv")

waist_fsh_p_values <- fsh_waist  %>%
  group_by(id) %>%
  summarise(sum_inds = sum(obs_less_ind))
waist_fsh_p_values$ppd_p_value <- round(waist_fsh_p_values$sum_inds/1000, 2)

g4 <- ggplot(waist_fsh_p_values, aes(x=ppd_p_value)) +
  geom_histogram(bins=50,fill="#A2D7D5", color="#49bab6") +
  labs(title="FSH posterior predictive p-values (Waist Circumference Model)",
       subtitle = expression(paste("Based on 1,000 draws of model posterior samples")), x="p-value", y="") + 
  theme(plot.title = element_text(hjust = 0.5,size=30, face="bold"),
        plot.subtitle = element_text(hjust = 0.5,size=30, face="bold"),
        axis.title = element_text(size=25, face="bold"),
        axis.text.x = element_text(size=20,vjust = 1, hjust = 1),
        axis.text = element_text(size=20, face="bold"))

grid.arrange(g3, g4, ncol=2)

## save as 2800 x 1400
grid.arrange(g1, g2, g3, g4, ncol=2)

fat_mass_ppd <- grid.arrange(g1, g2, ncol=2)

## create plots of 2 individuals 

e2_cis <- read.csv("/Users/irena/repos/swan_models/R/figures/plot_4_individuals_fat_0520.csv")
fsh_cis <- read.csv("/Users/irena/repos/swan_models/R/figures/plot_4_individuals_fsh_cis_0520.csv")
B_summary <- read.csv("/Users/irena/repos/swan_models/R/figures/B_summary_fat_mass_0520.csv")


plot_4_ids <- c(15, 82)
B_e2_subset  <- B_summary[B_summary$new_id%in%plot_4_ids&B_summary$hormone=="e2",c("mean", "beta_par", "new_id")]
B_fsh_subset  <- B_summary[B_summary$new_id%in%plot_4_ids&B_summary$hormone=="fsh",c("mean", "beta_par", "new_id")]

B_e2_wide <- reshape(B_e2_subset , idvar = "new_id", timevar = c("beta_par"), direction = "wide")
B_fsh_wide <- reshape(B_fsh_subset , idvar = "new_id", timevar = c("beta_par"), direction = "wide")

merged_ind_traj_e2 <- merge(B_e2_wide, e2_cis, by.x="new_id", by.y="id")

merged_ind_traj_e2$predicted_traj <- merged_ind_traj_e2$mean.b0 + (merged_ind_traj_e2$mean.b1*merged_ind_traj_e2$time_fmp)

merged_ind_traj_e2 <- merge(merged_ind_traj_e2, fat_hormone[,c("new_id","time_fmp" , "e2_loess_resid")], by = c("new_id", "time_fmp"))

merged_ind_traj_fsh <- merge(B_fsh_wide,fsh_cis,  by.x="new_id", by.y="id")
merged_ind_traj_fsh$predicted_traj <- merged_ind_traj_fsh$mean.b0 + (merged_ind_traj_fsh$mean.b1*merged_ind_traj_fsh$time_fmp)
merged_ind_traj_fsh <- merge(merged_ind_traj_fsh, fat_hormone[,c("new_id","time_fmp" , "fsh_loess_resid")], by=c("new_id", "time_fmp"))




bi_vars <- read.csv("/Users/irena/repos/swan_models/R/figures/bis_v_cov_2_individuals_0520.csv")
estimated_vars <- read.csv("/Users/irena/repos/swan_models/R/figures/posterior_means_vcov_2_individuals_0520.csv")

colnames(bi_vars)[1] <- "new_id"
colnames(estimated_vars)[5:6] <-c("perc_5", "perc_95")

merged_ind_traj_e2 <- merge(bi_vars[,c("new_id", "e2_intercept", "e2_slope", "e2_cov")], merged_ind_traj_e2, by="new_id")
merged_ind_traj_fsh <- merge(bi_vars[,c("new_id", "fsh_intercept", "fsh_slope", "fsh_cov")], merged_ind_traj_fsh, by="new_id")


## compute the variance of the mean trajectory (var(b_0i + b_1i*t))
merged_ind_traj_e2$var_e2_mean_traj <- merged_ind_traj_e2$e2_intercept + (((merged_ind_traj_e2$time_fmp)^2)*merged_ind_traj_e2$e2_slope) + 2*merged_ind_traj_e2$time_fmp*merged_ind_traj_e2$e2_cov

merged_ind_traj_fsh$var_fsh_mean_traj <- merged_ind_traj_fsh$fsh_intercept + (((merged_ind_traj_fsh$time_fmp)^2)*merged_ind_traj_fsh$fsh_slope) + 2*merged_ind_traj_fsh$time_fmp*merged_ind_traj_fsh$fsh_cov

merged_ind_traj_e2 <- merge(merged_ind_traj_e2, estimated_vars[estimated_vars$variable=="e2_var",c("new_id", "variable", "mean", "perc_5", "perc_95")], by="new_id")
merged_ind_traj_fsh <- merge(merged_ind_traj_fsh, estimated_vars[estimated_vars$variable=="fsh_var",c("new_id", "variable", "mean", "perc_5", "perc_95")], by="new_id")

## add the posterior means of the variances 
merged_ind_traj_e2$estimated_sigma <- sqrt(merged_ind_traj_e2$mean)
merged_ind_traj_fsh$estimated_sigma <- sqrt(merged_ind_traj_fsh$mean)

## the 5th and 95th order statistics for the variance estimates
merged_ind_traj_e2$estimated_5_sigma <- sqrt(merged_ind_traj_e2$perc_5)
merged_ind_traj_e2$estimated_95_sigma <- sqrt(merged_ind_traj_e2$perc_95)
merged_ind_traj_fsh$estimated_5_sigma <- sqrt(merged_ind_traj_fsh$perc_5)
merged_ind_traj_fsh$estimated_95_sigma <- sqrt(merged_ind_traj_fsh$perc_95)

## plot the prediction intervals 
g3 <- ggplot(data=merged_ind_traj_e2, aes(x=time_fmp, y=predicted_traj)) + 
  geom_line(color="#42427d", size=3.5) +
  #  scale_color_manual(values = c("e2" = "", "fsh"="#4b9c98")) + 
  geom_ribbon(aes(ymin=(predicted_traj+1.64*sqrt(var_e2_mean_traj)),ymax=(predicted_traj-1.64*sqrt(var_e2_mean_traj))), fill='#490880', alpha=0.35) + 
  geom_ribbon(aes(ymin=(predicted_traj+1.64*estimated_sigma),ymax=(predicted_traj-1.64*estimated_sigma)), fill='#c061ff', alpha=0.25) + 
  geom_ribbon(aes(ymin=predicted_traj-1.64*estimated_5_sigma, ymax=predicted_traj+1.64*estimated_5_sigma), size=1.5, alpha=0, linetype="dotted", color="#570b7a" ##b351f0
  ) + 
  geom_ribbon(aes(ymin=predicted_traj-1.64*estimated_95_sigma, ymax=predicted_traj+1.64*estimated_95_sigma), size=1.5, alpha=0, linetype="dotted", color="#570b7a") + 
  # geom_line(aes( aes(x=time_fmp, y=predicted_traj+1.64*estimated_5_sigma)), size=1.5, linetype="dashed", color="black") + 
  geom_vline(xintercept = 0, colour = "#FF4500", linetype="dashed", size=1.5) + 
  facet_wrap(~new_id) +
  geom_point(data=merged_ind_traj_e2, aes(x=time_fmp, y=e2_loess_resid,
                                          fill=as.factor(new_id)),shape=24, size=3.5) +
  scale_fill_manual(values=c("#d66418", "#f2168f")) +
  labs(
    #title="Estimated and Observed Hormone Trajectories for 2 Individuals", subtitle=paste0("Plotted with estimated mean trajectory standard deviations and estimated individual residual standard deviations"),
    x="", y="E2 Hormone Residuals") + 
  theme(plot.title=element_text(size=30, face="bold", hjust = 0.5),
        plot.subtitle = element_text(size=26, hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(size=20,angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=30))

g4 <- ggplot(data=merged_ind_traj_fsh, aes(x=time_fmp, y=predicted_traj)) + 
  geom_line(color="#4b9c98", size=3.5) +
  geom_ribbon(aes(ymin=(predicted_traj+1.64*sqrt(var_fsh_mean_traj)),ymax=(predicted_traj-1.64*sqrt(var_fsh_mean_traj))), fill='#49d17d', alpha=0.35) + 
  geom_ribbon(aes(ymin=(predicted_traj+1.64*estimated_sigma),ymax=(predicted_traj-1.64*estimated_sigma)), fill="#4b9c98", alpha=0.25) + 
  geom_ribbon(aes(ymin=predicted_traj-1.64*estimated_5_sigma, ymax=predicted_traj+1.64*estimated_5_sigma), size=1.5, alpha=0, linetype="dotted", color="#057869"##45b59f
  ) + 
  geom_ribbon(aes(ymin=predicted_traj-1.64*estimated_95_sigma, ymax=predicted_traj+1.64*estimated_95_sigma), size=1.5, alpha=0, linetype="dotted", color="#057869") + 
  geom_vline(xintercept = 0, colour = "#FF4500", linetype="dashed", size=1.5) + 
  facet_wrap(~new_id) +
  geom_point(data=merged_ind_traj_fsh, aes(x=time_fmp, y=fsh_loess_resid,
                                           fill=as.factor(new_id)),shape=22, size=3.5) +
  scale_fill_manual(values=c("#d66418", "#f2168f"))+
  labs(x="Time (in years) to Final Menstrual Period (FMP)", y="FSH Hormone Residuals") + 
  theme(legend.position = "none",
        axis.text.x = element_text(size=30,angle =45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=30))

# save as 1800 x 1600 
gridExtra::grid.arrange(g3,g4)
