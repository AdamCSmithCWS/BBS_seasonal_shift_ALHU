### modeling the seasonal pattern in BBS observations through time


# packages ----------------------------------------------------------------

library(bbsBayes)
library(tidyverse)
library(mgcv)


source("Functions/GAM_basis_function_mgcv.R")

# prep the BBS data -------------------------------------------------------

st_dat = stratify(by = "bbs_cws")

species = "Allen's Hummingbird"


spsData = prepare_data(st_dat,species_to_run = species,
                            model = "gamye",
                            heavy_tailed = TRUE)


#require(lubridate)
spsData$doy = lubridate::yday(as.Date(paste(spsData$r_year,
                            spsData$month,
                            spsData$day,
                            sep = "-")))

spsData$decadeF = cut(spsData$r_year,breaks = c(1965,
                                                1979.5,
                                                1989.5,
                                               1999.5,
                                               2009.5,
                                               2020),
                     labels = c("70s",
                                "80s",
                                "90s",
                                "00s",
                                "10s"),
                     ordered_result = TRUE)

spsData$decade = as.integer(spsData$decadeF)



# Explore the raw data of counts by day of year for each decade ---------------------------


tmp = data.frame(year = spsData$year,
                 route = factor(spsData$route),
                 obser = factor(spsData$obser),
                 count = spsData$count,
                 lcount = log(spsData$count+1,base = 10),
                 doy = spsData$doy,
                 decadeF = spsData$decadeF,
                 yr_d = as.integer(str_sub(spsData$r_year,4,4)),
                 week = cut(spsData$doy,breaks = seq(134,195,by = 7)))



bp = ggplot(data = tmp,aes(y = lcount,x = doy,colour = yr_d))+
  #geom_boxplot(varwidth = TRUE)+
  #geom_violin()+
  geom_point(alpha = 0.3,position = position_jitter(width = 0.25))+
  #scale_y_log10()+
  scale_colour_viridis_c(begin = 0.1,end = 0.7)+
  geom_smooth()+
  facet_wrap(~decadeF,nrow = 2)
print(bp)




# generating the GAM basis function data ----------------------------------


spsData$season <- spsData$doy-(min(spsData$doy)-1)

nseason= max(spsData$season)

season_gam = gam_basis(orig.preds = 1:nseason,
                    nknots = floor(nseason/7),
                    npredpoints = nseason,
                    sm_name = "season")


# export the base model from bbsBayes -------------------------------------


bbsBayes::model_to_file(model = "gamye",
                        filename = "models/gamye_season_blank.R",
                        heavy_tailed = TRUE)



# now, there's a manual step required to modify the JAGS code for the base model.
# once modified the "gamye_season_blank.R" is re-saved as
# "gamye_season.R" in the "models" directory
# the model estimates a non-linear smooth for the effect of season (day of the season)
# the smooths vary by decade, and the parameters of the smooth
# are fit with a first-difference time-series 
# beta_season[year,knot] ~ dnorm(beta_season[year-1,knot],taubeta_season)




# add the new data to the bbsBayes data list ------------------------------


spsData$season_basis <- season_gam$season_basis
spsData$nknots_season <- season_gam$nknots_season
spsData$decade <- as.integer(spsData$decadeF)
spsData$ndecades <- max(spsData$decade)
spsData$nseason = nseason

# names(spsData)
spsData$decadeF <- NULL

params = c("n",
          "n3",
          "seasoneffect",
          "beta_season",
          "B_season",
          "tau_B_season",
          "taubeta_season")



fit <- bbsBayes::run_model(jags_data = spsData,
                           model_file_path = "models/gamye_season.R",
                           parameters_to_save = params,
                           parallel = TRUE)

inds = generate_indices(fit,jags_data = spsData,alternate_n = "n3")
trends = generate_trends(indices = inds,Min_year = 1970)
inds2 = generate_indices(fit,jags_data = spsData,alternate_n = "n")
ip = plot_indices(inds2,min_year = 1970)


seas <- tidybayes::gather_draws(fit$samples,seasoneffect[decade,day])
seas_sum <- seas %>% group_by(.variable,decade,day) %>% 
  summarise(mean = mean(exp(.value)),
            lci = quantile(exp(.value),0.025),
            uci = quantile(exp(.value),0.975)) %>% 
  mutate(decadeF = factor(decade,
                          levels = 1:5,
                          labels = c("1966-1979",
                             "1980-1989",
                             "1990-1999",
                             "2000-2009",
                             "2010-2019"),
         ordered = TRUE))
seas_p = ggplot(data = seas_sum,aes(x = day,y = mean))+
  #geom_ribbon(aes(ymin = lci,ymax = uci,fill = decadeF),alpha = 0.05) + 
  geom_line(aes(colour = decadeF))+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(begin = 0.2,end = 0.8,aesthetics = c("colour","fill"))+
  ylab("Effect of season on counts (mean additional birds due to season)")+
  xlab("day of BBS season")+
  labs(title = "Seasonal pattern in counts by decade for Allen's Hummingbird")

seas_p_f = ggplot(data = seas_sum,aes(x = day,y = mean))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = decadeF),alpha = 0.2) + 
  geom_line(aes(colour = decadeF))+
  scale_y_continuous(limits = c(0,NA))+
  scale_colour_viridis_d(begin = 0.2,end = 0.8,aesthetics = c("colour","fill"))+
  ylab("Effect of season on counts (mean additional birds due to season)")+
  xlab("day of BBS season")+
  labs(title = "Seasonal pattern in counts by decade for Allen's Hummingbird")+
  theme(legend.position = "none")+
  facet_wrap(~decadeF,nrow = 2,scales = "fixed")


pdf(file = paste0("figures/",species,"_seasonal_effect_by_decade.pdf"))
print(seas_p)
print(seas_p_f)
for(i in length(ip):1){
  rr = unique(inds2$data_summary$Region)[i]
tmp1 = inds$data_summary
wtt = which(trends$Region == rr)
ttr = round(trends[wtt,"Trend"],1)
ttr1 = round(trends[wtt,"Trend_Q0.025"],1)
ttr2 = round(trends[wtt,"Trend_Q0.975"],1)

lbl = paste0(ttr," [",ttr1," : ",ttr2,"]")
tmp = tmp1 %>% filter(Region == rr,Year >= 1970)
yup = 0.9*max(tmp1$Index_q_0.975)
ggadd <- ip[[i]] + geom_ribbon(data = tmp,aes(x = Year, ymin = Index_q_0.025,ymax = Index_q_0.975),
                               alpha = 0.2,fill = "darkorange")+
  geom_line(data = tmp,aes(x = Year, y = Index),colour = "darkorange")+
  annotate("text",x = 2000,y = yup*0.7,label = lbl)
print(ggadd)
}
dev.off()


