####install libraries####
# { 
#   install.packages("nlme")
#   install.packages("data.table")
#   install.packages("emmeans")
#   install.packages("lmerTest")
#   install.packages("lattice")
#   install.packages("plyr")
#   install.packages("ggplot2")
#   install.packages("Renext")
#   install.packages("ggsn")
#   install.packages("doBy")
#   install.packages("Rmisc")
#   install.packages("dplyr")
#   install.packages("gmodels")
#   install.packages("googlesheets")
#   install.packages("httpuv")
#   install.packages("RCurl")
#   install.packages("gsheet")
#   install.packages("standardize")
#   install.packages("lme4")
#   install.packages("Matrix")
#   }

####used libraries####
{
  library(nlme)
  library(data.table)
  library(emmeans)
  library(lmerTest)
  library(lattice)
  library(plyr)
  library(ggplot2)
  library(Renext)
  library(ggsn)
  library(doBy)
  library(Rmisc)
  library(dplyr)
  library(gmodels)
  library(googlesheets)
  library(httpuv)
  library(RCurl)
  library(gsheet)
  library(standardize)
  library(lme4)
  library(Matrix)
}

#Dropbox (Neurasmus BV) Balance Beam test

# set wd - working directory
rm(wd)
wd <- setwd("~/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/balancebeam")
wd <- setwd("~/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/balancebeam")
print(wd)


# read csv file, select data, convert to dataframe, recode data
#csvfile <- "BB Total Data.csv"
csvfile <- "BB Total Data.csv"

BBRaw <- read.csv(csvfile)
View(BBRaw)


# convert CSV to data.frame
BBRaw <- data.frame(BBRaw)

write.csv(BBRaw, file = "BBRaw.csv")

# new data frame without NA, NA removed
BBRaw_na_omit <- na.omit(BBRaw)
View(BBRaw_na_omit)

write.csv(BBRaw_na_omit, file = "BBRaw_na.csv")


# Plot1, ggplot geom_boxplot + geom_dotplot ALL BBRaw data points
plot1 <- ggplot(BBRaw_na_omit, aes(x=factor(Genotype, level=c("Control", "CE")),y=Time)) +
  facet_wrap(~ Beam) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1.1, dotsize=0.7, width = 3, aes(fill = Genotype)) +
  scale_x_discrete (name = "Genotype") +
  scale_y_continuous (name = "Time to cross beam (s.)") +
  expand_limits(y=c(0,30)) +
  theme_bw() +
  ggtitle ("Balance Beam Test") +
  theme(legend.text = element_text(size = 12)) +
  guides(col = guide_legend(ncol = 1, byrow = TRUE, keywidth = 0.1, keyheight = 0.1, default.unit = "inch")) +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.y  = element_text(size=14),
        axis.line = element_line(colour = 'black', size = 0.3), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

plot(plot1)

# Linear Mixed-Effects Models
LME1<-lme(Time ~ Genotype*Beam,
          data = BBRaw, 
          correlation = NULL,
          random = ~1 | AnimalID, 
          method = "REML", 
          na.action = na.exclude)

summary(LME1)
anova(LME1)
VarCorr(LME1)
plot(LME1)
qqnorm(LME1)



# Plot2, Sex differences ggplot geom_boxplot + geom_dotplot ALL BBRaw data points
plot2 <- ggplot(BBRaw_na_omit, aes(x=factor(Genotype, level=c("Control", "CE")),y=Time)) +
  facet_grid(~ Beam + Sex) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1.1, dotsize=0.5, width = 3, aes(fill = Genotype)) +
  scale_x_discrete (name = "") +
  scale_y_continuous (name = "Time to cross beam (s.)") +
  expand_limits(y=c(0,30)) +
  theme_bw() +
  ggtitle ("") +
  theme(legend.text = element_text(size = 12)) +
  guides(col = guide_legend(ncol = 1, byrow = TRUE, keywidth = 0.1, keyheight = 0.1, default.unit = "inch")) +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.y  = element_text(size=14),
        axis.line = element_line(colour = 'black', size = 0.3), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

plot(plot2)

# Linear Mixed-Effects Models sex difference 
LME2<-lme(Time ~ Genotype*Beam*Sex,
          data = BBRaw, 
          correlation = NULL,
          random = ~1 | AnimalID, 
          method = "REML", 
          na.action = na.exclude)

summary(LME2)
anova(LME2)
VarCorr(LME2)
plot(LME2)
qqnorm(LME2)


#-------------------------------------------------------------------------------------------------------------------------------------------
# Data frame with the variables of only the 6mm beam
beam6 <- BBRaw_na_omit %>% 
  filter(Beam == 6)

# LME for the 6mm beam Time ~ Group
LME1<-lme(Time ~ Genotype,
          data = beam6, 
          correlation = NULL,
          random = ~1 | AnimalID, 
          method = "REML", 
          na.action = na.exclude)

summary(LME1)
anova(LME1)
VarCorr(LME1)
plot(LME1)
qqnorm(LME1)


# Data frame with the variables of only the 12mm beam
beam12 <- BBRaw_na_omit %>% 
  filter(Beam == 12)

# LME for the 12mm beam Time ~ Group
LME1<-lme(Time ~ Genotype,
          data = beam12, 
          correlation = NULL,
          random = ~1 | AnimalID, 
          method = "REML", 
          na.action = na.exclude)

summary(LME1)
anova(LME1)
VarCorr(LME1)
plot(LME1)
qqnorm(LME1)




# make BBRawAggr file 1 which is aggregate of aggregate
rm(BBRaw)
{BBRawAggr <- BBRaw_na_omit %>% 
    group_by(Genotype, Construct, AnimalID,Beam) %>% 
    summarise_at(c("Time"), funs( 
  mean (., na.rm=T),
  median (., na.rm=T),
  #n = sum(!is.na(.)),
  #min = min (.,na.rm=T),
  #max = max (.,na.rm=T),
  sd = sd (.,na.rm=T),
  se = sd(., na.rm=T)/sqrt(sum(!is.na(.))), 
  #cv = sd(., na.rm=T)/mean (., na.rm=T),
  # meanCI = ci (., na.rm=T)[1],
  #lowCI = ci (., na.rm=T)[2],
  #highCI = ci (., na.rm=T)[3],
  # sdCI = ci (., na.rm=T)[4]
  ci95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2))
}

View(BBRawAggr)

write.csv(BBRawAggr, file = "BBRawAggr.csv")


#plot3, ggplot geom_boxplot + geom_dotplot mean Time of BBRawaggr
plot3 <- ggplot(BBRawAggr, aes(x=factor(Genotype, level=c("Control", "CE")),y=mean)) +
  facet_wrap(~ Beam) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1.1, dotsize=0.8, width = 3, aes(fill = Genotype)) +
  scale_x_discrete (name = "") +
  scale_y_continuous (name = "Time to cross beam (s.)") +
  #expand_limits(y=c(0,20)) +
  theme_bw() +
  ggtitle ("Balance Beam Test") +
  theme(legend.text = element_text(size = 12)) +
  guides(col = guide_legend(ncol = 1, byrow = TRUE, keywidth = 0.1, keyheight = 0.1, default.unit = "inch")) +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.y  = element_text(size=14),
        axis.line = element_line(colour = 'black', size = 0.3), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

plot(plot3)

# Plot3, ggplot geom_bar
#ggplot(BBRawAggr, aes(x = factor(Genotype, level=c("Control", "CE")), fill = Genotype)) + 
  theme_bw() +
  facet_wrap(~ Beam) +
  geom_bar() +
  labs(y = " Time to cross beam(s.)",
       title = "Balance Beam Test",
       x = "")

# Linear Mixed-Effects Models
LME1<-lme(mean ~ Genotype*Beam,
          data = BBRawAggr, 
          correlation = NULL,
          random = ~1 | AnimalID, 
          method = "REML", 
          na.action = na.exclude)

summary(LME1)
anova(LME1)
VarCorr(LME1)
plot(LME1)
qqnorm(LME1)


# make BBRawAggr1 file 2 which is aggregate of aggregate
rm(BBRawAggr)
{BBAggr1 <- BBRaw_na_omit %>% 
    group_by(Genotype, Beam) %>% 
    summarize_at (c("Time"), funs( 
  mean (., na.rm=T),
  median (., na.rm=T),
  #n = sum(!is.na(.)),
  #min = min (.,na.rm=T),
  #max = max (.,na.rm=T),
  sd = sd (.,na.rm=T),
  se = sd(., na.rm=T)/sqrt(sum(!is.na(.))), 
  #cv = sd(., na.rm=T)/mean (., na.rm=T),
  # meanCI = ci (., na.rm=T)[1],
  #lowCI = ci (., na.rm=T)[2],
  #highCI = ci (., na.rm=T)[3],
  # sdCI = ci (., na.rm=T)[4]
  ci95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2))
}

View(BBAggr1)

write.csv(BBAggr1, file = "BBAggr1.csv")





 
